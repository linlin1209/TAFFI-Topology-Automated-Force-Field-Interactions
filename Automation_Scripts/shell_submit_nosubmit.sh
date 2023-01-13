#!/bin/bash

# Collect a bash-safe pathway to the current folder
# Inserts escape characters before parentheses so that the script can
# handle filenames/folders that involve these
function esc_str () {
    a=${1}
    a=$( echo ${a} | sed 's/(/\\(/' )
    a=$( echo ${a} | sed 's/)/\\)/' )
    echo ${a}
}

if [ $# -lt 2 -o "$1" == "?" ]; then
    echo -e "\nDescription: supply a minimum of an input string that holds the shell call (e.g., 'python foo.py -arg1 1 -arg2 2 ...')"
    echo -e "and an output filename."
    echo -e "\nOptional arguments:\n"
    echo -e "\t-p\t: number of processors (integer; nodes are automatically scaled to match desired cores)"
    echo -e "\t-t\t: walltime (integer; units of hours)"
    echo -e "\t-q\t: queue name (debug;serial;CLUSTER)"
    echo -e "\t-ppn\t: processors per node"
    echo -e "\t-d\t: duplicate job (integer; the job is submitted in duplicate with sequential dependecies)"
    echo -e "\t-h\t: submit job in hold state (no argument)"
    echo -e "\t-shell\t: any arguments sandwiched between -shell flags (e.g., -shell mv foo folder/foo -shell) "
    echo -e "\t\t  will be added to submission script and executed in the run before the python call."
    echo -e ""
    exit
fi

# Parse required positional arguments
command=$1
filename_prefix=$(echo "$1" | awk -F '[.]' '{print $1}')
filename_tail=$(echo "$1" | awk -F '[.]' '{print $NF}'); shift
keep_opt='off';
outputname=$1; shift
name_sub=$( echo ${outputname} | sed 's/(//' )
name_sub=$( echo ${name_sub} | sed 's/)//' )
outputname=$( esc_str ${outputname} )

# List defaults
where_is_this_being_run='tachus'
num_procs=1
queue='CLUSTER'
walltime=4
hold_opt=0
ppn=24;  # processors per node
duplicates=1;
shell_command="";
id="";
min_flag=0;

# Parse command-line inputs
while [ $# -gt 0 ]; do

    if [ $1 = "-w" ]; then
	shift;
	where_is_this_being_run=$1; shift;
    elif [ $1 = "-p" ]; then
	shift;
	num_procs=$1; shift;	
    elif [ $1 = "-k" ]; then
	shift;
	keep_opt='on'
    elif [ $1 = "-q" ]; then
	shift;
	queue=$1; shift;
    elif [ $1 = "-t" ]; then
	shift;
	walltime=$1; shift;
    elif [ $1 = "-h" ]; then
	shift;
	hold_opt=1; shift;
    elif [ $1 = "-d" ]; then
	shift;
	duplicates=$1; shift;
    elif [ $1 = "-id" ]; then
	shift;
	id=$1; shift;
    elif [ $1 = "-ppn" ]; then
	shift;
	ppn=$1; shift;
    elif [ $1 = "-shell" ]; then
	shift;
	while [ $1 != "-shell" ]; do
	    shell_command="${shell_command} $1"; shift;
	done
	shift;
    else
	echo -e "ERROR: command-line argument \"${1}\" is not recognized. Exiting..."
	exit
    fi
done

# Check for minute argument 
if [ ! -z $( echo ${walltime} | grep "min" ) ]; then
    walltime=$( echo ${walltime} | sed 's/min//' )
    min_flag=1
fi

# Remove old submission scripts
if [ -f ${outputname}.submit ]; then
    echo ${outputname}
    rm ${outputname}.submit
fi

# Save the working directory to variable
current_dir=$(pwd)
current_dir=$( esc_str "${current_dir}" )

# calculate the number of nodes required 4*num_procs on somasim
# since jobs are always run in fours. 
num_nodes=$(( ( ${num_procs} + ${ppn} - 1 ) / ${ppn} ))

# Calculate the tasks per node (i.e., number of cores per node)
if [ ${num_procs} -gt ${ppn} ]; then
    tasks_per_node=${ppn}
else
    tasks_per_node=${num_procs}
fi

# Write the submission script
echo -e "#!/bin/bash" > ${name_sub}.submit
echo -e "#" >> ${name_sub}.submit
echo -e "#PBS -N ${outputname}" >> ${name_sub}.submit 
echo -e "#PBS -A chm114" >> ${name_sub}.submit 
if [ ${ppn} -eq 0 ]; then
    echo -e "#PBS -l nodes=${num_nodes}" >> ${name_sub}.submit 
else
    echo -e "#PBS -l nodes=${num_nodes}:ppn=${ppn}" >> ${name_sub}.submit 
fi
if [ ${min_flag} -eq 1 ]; then
    echo -e "#PBS -l walltime=00:${walltime}:00" >> ${name_sub}.submit 
else
    echo -e "#PBS -l walltime=${walltime}:00:00" >> ${name_sub}.submit 
fi
echo -e "#PBS -q ${queue}" >> ${name_sub}.submit
echo -e "#PBS -S /bin/sh" >> ${name_sub}.submit 
echo -e "#PBS -o ${outputname}.out" >> ${name_sub}.submit 
echo -e "#PBS -e ${outputname}.err" >> ${name_sub}.submit 

echo -e "\n# cd into the submission directory" >> ${name_sub}.submit
echo -e "cd -- ${current_dir}" >> ${name_sub}.submit
echo -e "echo Working directory is ${current_dir}" >> ${name_sub}.submit
echo -e "echo Running on host `hostname`" >> ${name_sub}.submit
echo -e "echo Time is `date`\n" >> ${name_sub}.submit

echo -e "# Copy the local path" >> ${name_sub}.submit
echo -e "PATH=${PATH}\n" >> ${name_sub}.submit

# If a shell command is supplied it is embedded before the lammps run
if [ "${shell_command}" != ""  ]; then
    echo -e "${shell_command}\n" >> ${name_sub}.submit
fi

echo -e "# Run python script" >> ${name_sub}.submit 
echo -e "sleep 2" >> ${name_sub}.submit
echo -e "${command}" >> ${name_sub}.submit
echo -e "sleep 2\n" >> ${name_sub}.submit

if [ ${duplicates} -eq 1 ]; then
    if [ ${hold_opt} -eq 1 ]; then      # Submit with a hold
	#qsub ${name_sub}.submit -h    
    elif [ "${id}" != ""  ]; then       # Submit with job dependency	
	#qsub -W depend=afterany:${id} ${name_sub}.submit
    else                                # Submit with no dependence
	#qsub ${name_sub}.submit
    fi

elif [ ${duplicates} -gt 1 ]; then  

    if [ ${hold_opt} -eq 1 ]; then      # Submit with a hold
    JOB=1
	#JOB=$(qsub ${name_sub}.submit -h)
    elif [ "${id}" != ""  ]; then       # Submit with job dependency
    JOB=1
	#JOB=$(qsub -W depend=afterany:${id} ${name_sub}.submit)
    else                                # Submit with no dependence
    JOB=1
	#JOB=$(qsub ${name_sub}.submit)
    fi
    printf "Submitting job:\t%s\n" "${JOB}"
    
    for i in $(seq $((${duplicates}-1))); do
	#JOB=$(qsub -W depend=afterany:${JOB} ${name_sub}.submit )
	printf "Submitting job:\t%s\n" "${JOB}"
    done
fi
