#!/bin/bash

# IMPORT SYSTEM SPECIFIC PATHS AND BATCH ARGUMENTS
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # Saves the location of this script to a variable


# Write command line parser
# List defaults
config_file="${DIR}/config.txt"
FF=""
solvents=()
run_folder="VDW_Params"
base="temp"
# Parse command-line inputs
while [ $# -gt 0 ]; do

    if [ $1 = "-mol" ]; then
	shift;
	solvents=( $1 ); 
	base=${1};
   echo $solvents
   echo $base
	shift;
    elif [ $1 = "-FF" ]; then
	shift;
	FF=$1; shift;	
    elif [ $1 = "-f" ]; then
	shift;
	run_folder="${1}"; shift;
    elif [ $1 = "-c" ]; then
	shift;
	config_file="${1}"; shift;
    else
	echo -e "ERROR: command-line argument \"${1}\" is not recognized. Exiting..."
	exit
    fi
done

# If the force-field is empty then it is read from the config file
if [ -z "$FF" ]; then
    FF=$( awk -v v="FF" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
fi

# Check that the configuration file exists
if [ ! -f ${config_file} ]; then
    echo "ERROR: the file ${config_file} does not exist."
    exit
fi

TAFFI_PATH=${DIR}/../                                     # It is easier to define a lot of the commands relative to the taffi root
PYTHON_ABS=$(which python)                                # Saves the absolute path to the favored python dist to variable
ORCA_SUBMIT=${DIR}/orca_submit.py
SHELL_SUBMIT=${DIR}/shell_submit.sh
LAMMPS_SUBMIT=${DIR}/lammps_submit.py
LAMMPS_EXE=$( awk -v v="LAMMPS_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
ORCA_EXE=$( awk -v v="ORCA_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )

VDW_MD_PROCS=$( awk -v v="VDW_MD_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )   
VDW_MD_WT=$( awk -v v="VDW_MD_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )      
VDW_MD_Q=$( awk -v v="VDW_MD_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )       
VDW_MD_NPP=$( awk -v v="VDW_MD_NPP" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     
VDW_MD_SCHED=$( awk -v v="VDW_MD_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )   
VDW_MD_PPN=$( awk -v v="VDW_MD_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )     
VDW_MD_SIZE=$( awk -v v="VDW_MD_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 

VDW_QC_PROCS=$( awk -v v="VDW_QC_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
VDW_QC_WT=$( awk -v v="VDW_QC_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
VDW_QC_Q=$( awk -v v="VDW_QC_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
VDW_QC_SCHED=$( awk -v v="VDW_QC_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
VDW_QC_PPN=$( awk -v v="VDW_QC_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
VDW_QC_SIZE=$( awk -v v="VDW_QC_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 
################

solvent_path="."
d_N_thresh="0.01"   # number density threshold for rerunning MD at NVT
d_N_fix="0.05"      # jobs whose number density falls below d_N_thresh are resubmitted with a number density of d_N_fix
L2_s="0.1"
L2_e="1.0"
max_cycles=20
max_attempts=10     # Maximum number of md resubmissions to attempt. 
sleep 5

# Waits for completion of jobs in the queue, expects an 
# array of job ids to be supplied as the first argument
function check_queue () {
    arr=("$@")
    result=('1')
    count=1
    while [ ${#result[@]} -gt 0 ]; do
	tmp=$( python ${DIR}/qstat.py | grep ${USER} | awk '{ if ( NF > 4 && $5 != "C" ) { print $1 }}')
	readarray -t ids_curr <<< "${tmp}"
	l2=" ${ids_curr[*]} "                    # add framing blanks
	result=()
	for item in ${arr[@]}; do
	    if [[ $l2 =~ " $item " ]] ; then    # use $item as regexp
		result+=($item)
	    fi
	done	
	if [ ${#result[@]} -gt 0 ]; then
	    if [ ${count} == 1 ]; then
		echo "waiting for jobs (${result[*]}) to complete...."
		count=$(( ${count} +  1 ))
	    fi	    
	    sleep 100
	fi
    done
} 

# If needed create the vdw directory
#if [ ! -d ${run_folder} ]; then
#   echo "Error ${run_folder} doesn't exist"
#    mkdir ${run_folder}
#fi
#cd ${run_folder}
pwd


# Determine the starting cycle if jobs are already present
if [ -d vdw ]; then

    # Find the number of cycle folders
    counter=0
    for i in vdw/cycle*; do
	if [ ${i} == "vdw/cycle*" ]; then
	    break
	else
	    counter=$(( ${counter} + 1 ))
	fi
    done

    # Determine the starting cycle
    cycle=1
    rem_flag=0

    # If no cycle folders are present, start on cycle 1 and delete old folders
    if [ ${counter} -eq 0 ]; then
	cycle=1

	# remove first md folder if incomplete
	if [ -d vdw/md-1 ]; then
	    if [ ! -f vdw/md-1/md-1.end.data ]; then
		rm -r vdw/md-1
	    fi
	fi

	# remove first md folder if incomplete
	if [ -d vdw/md-2 ]; then
	    if [ ! -f vdw/md-2/md-2.end.data ]; then
		rm -r vdw/md-2
	    fi
	fi

	# remove third md folder if incomplete
	if [ -d vdw/md-3 ]; then
	    if [ ! -f vdw/md-3/md-3.end.data ]; then
		rm -r vdw/md-3
	    fi
	fi

	# remove configs folder if present
	if [ -d vdw/configs ]; then
	    rm -r vdw/configs
	fi

    # Else, start after the first complete cycle
    else
	
	for i in `seq ${max_cycles}`; do

	    # Check the status of existing cycles
	    # If incomplete, restart here and delete the rest
	    md=$(( ${i} + 2 ))
	    if [ ! -f vdw/cycle-${i}/DFT-AA.db -o ! -f vdw/md-${md}/md-${md}.end.data -o ! -d vdw/configs/${md}-0 ]; then

		# loop over md folders and removed any that are greater than or equal to the current md
		for j in vdw/md-*/; do		    
		    num=$( python -c "print (int('${j}'.split('/')[-2].split('-')[-1]))" )
		    if [ ${num} -ge ${md} ]; then
			echo "removing md-${num}... (rm -r ${j})"
			rm -r ${j}

			# Remove MD configs for this cycle if present
			for k in vdw/configs/${num}-* ; do
			    if [ k == "${base}_vdw/configs/${num}-*" ]; then
				break
			    else
				rm -r vdw/configs/${num}-*
				break
			    fi
			done
#			if [ -d ${base}_vdw/configs/${num}-0 ]; then
#			    echo "Deleting ${base}_vdw/configs/${num}-*"
#			    rm -r ${base}_vdw/configs/${num}-*
#			fi
		    fi
		done

		# loop over cycle folders and remove any that are greater than or equal to the current cycle
		for j in vdw/cycle-*/; do		    
		    num=$( python -c "print int('${j}'.split('/')[-2].split('-')[-1])" )
		    if [ ${num} -ge ${cycle} ]; then
			echo "removing cycle-${num}... (rm -r ${j})"
		        rm -r ${j}
		    fi
		done

		# Break out of checking folders
		break
	    fi

	    # If complete up to this point, increment starting cycle
	    cycle=$(( ${cycle} + 1 ))

	done
    fi

# If no job data is present then start with 1
else
    cycle=1
fi

# Outer loop controls cycles
fix_d=0
while [ ${cycle} -le ${max_cycles} ]; do 

    # Step 1 - Run LAMMPS jobs to sample configurations
    # while loop handles resubmission conditions if the number density of the md simulation drops below d_N_thresh
    cond=0
    attempts=1
    while [ ${cond} -eq 0 ]; do

	# Check if the maximum number of resubmission attempts has been reached
        if [ "${attempts}" -gt "${max_attempts}" ]; then
            echo "ERROR: md failed. Maximum number of resubmission attempts reached."
	    exit
        fi
	
	ids_sub=()
	if [ ${cycle} -eq 1 ]; then

	    # Collate the preliminary intramolecular modes 
	    if [ ! -f intermediate_params.db ]; then
		python ${TAFFI_PATH}FF_functions/merge_FF.py intermediate_params.db ./Intra/after_charges/Intra-DFT.db -only "atom bond angle torsion charge"
	    fi

	    # STEP 1: RUN MD (cycle == 1)
	    #for i in ${solvents[*]}; do

		N_atoms=$(cat "./out_inter.xyz" | wc -l | bc)
		N_atoms=$(echo "${N_atoms} - 2" | bc)
		N_mol=$(echo "1000/${N_atoms} + 1" | bc)

		# Make runfolder
		if [ ! -d vdw ]; then
		    mkdir vdw
		fi

		
		# For resubmissions use a lower density
		if [ "${attempts}" -gt 1 ]; then
		    ts=1.0
		    d_N=0.08
		    t_A="1E5"
		else
		    ts=1.0
		    d_N=0.09
		    t_A="1E5"
		fi

		# Make runs (if/then are in place in case of resubmission)
		if [ ! -d vdw/md-1 ]; then
		    python ${TAFFI_PATH}FF_functions/gen_md_for_sampling.py "./out_inter.xyz" -FF "intermediate_params.db ./charges/CHELPG_calcs/charges/fit_charges.db ${FF}" \
			   -N ${N_mol} -o vdw/md-1 -d_N "${d_N}" -sigma_scale 0.6 -charge_scale 0.0 -T_A 10 -T 298 -t_A "${t_A}" -t 1E5 --UFF -q 0 -ts "${ts}" --molecule > /dev/null
		fi

		if [ ! -d vdw/md-2 ]; then
		    python ${TAFFI_PATH}FF_functions/gen_md_for_sampling.py "./out_inter.xyz" -FF "intermediate_params.db ./charges/CHELPG_calcs/charges/fit_charges.db ${FF}" \
		       -N ${N_mol} -o vdw/md-2 -d_N "${d_N}" -sigma_scale 0.8 -charge_scale 0.0 -T_A 10 -T 298 -t_A "${t_A}" -t 1E5 --UFF -q 0 -ts "${ts}" --molecule > /dev/null
		fi

		if [ ! -d vdw/md-3 ]; then
		    python ${TAFFI_PATH}FF_functions/gen_md_for_sampling.py "./out_inter.xyz" -FF "intermediate_params.db ./charges/CHELPG_calcs/charges/fit_charges.db ${FF}" \
		       -N ${N_mol} -o vdw/md-3 -d_N "${d_N}" -sigma_scale 1.0 -charge_scale 0.0 -T_A 10 -T 298 -t_A "${t_A}" -t 1E5 --UFF -q 0 -ts "${ts}" --molecule > /dev/null
		fi
		
		# If/else loops bracketing the directory switches have been added throughout the script for graceful failure. 
		if [ -d vdw ]; then
		    cd vdw
		else
		    echo "ERROR in vdw_loop: failed to create initial md jobs. Exiting..."
		    exit
		fi
		
		# Submit the jobs
		for j in md*/; do
		    cd ${j}
		    name=$( echo $j | awk -F '[/]' '{ print $1 }' )
		    echo "Submitting ${name}..."
          ids_sub+=( $(python ${LAMMPS_SUBMIT} -f ${name}.in.init -p ${VDW_MD_PROCS} -t ${VDW_MD_WT} -q ${VDW_MD_Q} -o ${name} -sched ${VDW_MD_SCHED} -ppn ${VDW_MD_PPN} -size ${VDW_MD_SIZE} -npp ${VDW_MD_NPP} -path_to_exe ${LAMMPS_EXE} --silent) )
		    cd ..
		done
		cd .. #Back to inchi dir
	   #done

	# STEP 1: RUN MD (cycle>1)
	else

	    md=$(($cycle + 2))
	    #for i in ${solvents[*]}; do    

		N_atoms=$(cat "./out_inter.xyz" | wc -l | bc)
		N_atoms=$(echo "${N_atoms} - 2" | bc)
		N_mol=$(echo "1000/${N_atoms} + 1" | bc)

		# If the constant number density flag has been triggered, run NVT
		if [ ${fix_d} -eq 1 ]; then
		    python ${TAFFI_PATH}FF_functions/gen_md_for_sampling.py "vdw" -N ${N_mol} -o vdw -T 298 -t_A 1E5 -t 1E5 -T_A 100 -d_N ${d_N_fix} --molecule > /dev/null

		# Else, run NPT
		else
		    python ${TAFFI_PATH}FF_functions/gen_md_for_sampling.py "vdw" -N ${N_mol} -o vdw -T 298 -t_A 1E5 -t 1E5 -T_A 100 --molecule  > /dev/null
		fi

		# Submit the MD simulation to the cluster (ids_sub is used to keep track of jobs awaiting completion)
		# If/else loops bracketing the directory switches have been added throughout the script for graceful failure.  
		if [ -d vdw/md-${md} ]; then
		    cd vdw/md-${md}
		else
		    echo "ERROR in vdw_loop: failed to create vdw/md-${md}. Exiting..."
		    exit
		fi

		# Check if input file was actually generated
		if [ ! -f md-${md}.in.init ]; then
		    echo "ERROR in vdw_loop: failed to create vdw/md-${md}.in.init.  Exiting..."    
		    exit
		fi

		name=md-${md}
		echo "Submitting ${name}..."
		ids_sub+=( $(python ${LAMMPS_SUBMIT} -f ${name}.in.init -p ${VDW_MD_PROCS} -t ${VDW_MD_WT} -q ${VDW_MD_Q} -o ${name} -sched ${VDW_MD_SCHED} -ppn ${VDW_MD_PPN} -size ${VDW_MD_SIZE} -npp ${VDW_MD_NPP} -path_to_exe ${LAMMPS_EXE} --silent) )

		cd ../.. #Back to inchikey dir
	    #done
	fi

	# Wait until jobs complete
	check_queue "${ids_sub[@]}"

	# Check if the jobs completed / need to be resubmitted
	#for i in ${solvents[*]}; do    	

	    # First cycle, check that the job completed
	    if [ ${cycle} -eq 1 ]; then
		cond=1
		folders=("md-1" "md-2" "md-3")
		for j in ${folders[*]}; do
		    if [ ! -f vdw/${j}/${j}.end.data ]; then
			echo "WARNING: vdw/${j} didn't complete. Attempting resubmission..."
			cond=0
			rm -r vdw/${j}
		    fi
		done	   

	    # Other cycles, check that the system hasn't evaporated and that the job completed
	    else

		# Check the number density in the latest MD simulation
		# Since in the first cycle the number density is held constant, the break condition is guarranteed to be met
		md=$((${cycle} + 2))
		d_N=$(python ${TAFFI_PATH}FF_functions/check_density.py vdw/md-${md}/md-${md}.end.data)
		if (( $(echo "0.01 > ${d_N}" |bc -l) )); then 	    

		    # Print diagnostic
		    echo "Switching over to NVT simulations at a constant number density of ${d_N_fix} atoms/A^3..." 

		    # Toggle NVT flag, remove the most recent run, backup the penultimate *end.data file, 
		    # and rescale the box in the penultimate *end.data file for use in the resubmitted job
		    fix_d=1
		    rm -r vdw/md-${md}
		    md_old=$((${cycle} +1 ))
		    cp vdw/md-${md_old}/md-${md_old}.end.data vdw/md-${md_old}/md-${md_old}.end.data.unscaled
		    python ${TAFFI_PATH}FF_functions/rescale_data.py vdw/md-${md_old}/md-${md_old}.end.data -d_N ${d_N_fix}

		elif [ ! -f vdw/md-${md}/thermo.avg -o ! -f vdw/md-${md}/equil.lammpstrj ]; then
		    echo "WARNING: md-${md} didn't complete. Attempting resubmission..."
		    cond=0
		    rm -r vdw/md-${md}	    
		else
		    cond=1
		fi
	    fi
	#done
	    
	# Increment attempts
	attempts=$(( ${attempts} + 1 ))

    done
	
    # STEP 2: RUN ORCA JOBS
    ids_sub=()
    #for i in ${solvents[*]}; do
    	python ${TAFFI_PATH}FF_functions/gen_jobs_for_vdw.py vdw -r_min_scale 1.5 -every 100 -N 50 -p ${VDW_QC_PROCS} -QC_types dft -FF "${FF}"

	# If/else loops bracketing the directory switches have been added throughout the script for graceful failure.
	if [ -d vdw ]; then
    	    cd vdw
	else
	    echo "ERROR in vdw_loop: failed to create vdw. Exiting..."
	    exit
	fi

    	tmp=$(echo $(python ${ORCA_SUBMIT} -p ${VDW_QC_PROCS} -t ${VDW_QC_WT} -ppn ${VDW_QC_PPN} -q ${VDW_QC_Q} -sched ${VDW_QC_SCHED} -size ${VDW_QC_SIZE} -path_to_exe "${ORCA_EXE}" --silent) | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
    	readarray -t ids_sub <<< "${tmp}"

    	cd .. # Bakc to inchikey folder
    #done

    # Wait until jobs complete
    echo "start_check"
    check_queue "${ids_sub[@]}"
    echo "over check"

    # remove *.gbw files (take up a LOT of memory)
    #for i in ${solvents[*]}; do
	cd vdw
	N_gbw=$( find -name *.gbw | wc -l )
	if [ ${N_gbw} != 0 ]; then
	    find -name *.gbw | xargs rm
	fi
	N_prop=$( find -name *.prop | wc -l )
	if [ ${N_prop} != 0 ]; then
	    find -name *prop | xargs rm
	fi    
	cd ..
    #done

    # STEP 3: FIT VDW PARAMETERS
    ids_sub=()
    #for i in ${solvents[*]}; do
	
	# Use previous cycle as an initial guess
	ids_sub+=( $( ${SHELL_SUBMIT} "${PYTHON_ABS} ${TAFFI_PATH}FF_functions/extract_vdw.py -f vdw -FF_DFT './charges/CHELPG_calcs/charges/fit_charges.db ${FF}' -o cycle-${cycle} -E_max 0.0 -xhi2_thresh 1E-8 -q_a 0.0 -q_b 0.0 -mixing_rule wh -L2_sigma ${L2_s} -L2_eps ${L2_e}" vdw_parse -p 1 -t 4 -q ${VDW_QC_Q} ) )

    #done
    cycle=$(( ${cycle} + 1 ))

    # Wait until jobs complete
    check_queue "${ids_sub[@]}"

    # Update convergence condition
    if [ ${cycle} -le ${max_cycles} ]; then 
	echo "incomplete" > vdw/status.txt
    else
	echo "complete" > vdw/status.txt
    fi

done

# Plot parameter convergence
if [ ! -f vdw.eps.pdf -o ! -f vdw.sigmas.pdf ]; then
    ${PYTHON_ABS} ${TAFFI_PATH}Parsers/plot_vdw_convergence.py vdw
fi

# Merge final vdw parameters and fit params
#cd ..

# Merge modes and charges into intermediate param dictionary
${PYTHON_ABS} ${TAFFI_PATH}FF_functions/merge_FF.py final_vdw.db vdw/cycle-${max_cycles}/DFT-AA.db -R
${PYTHON_ABS} ${TAFFI_PATH}FF_functions/merge_FF.py final_vdw.db vdw/cycle-${max_cycles}/DFT-UA.db -R

