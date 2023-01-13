#!/bin/bash

# SYSTEM SPECIFIC PATHS
######################
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"   # Saves the location of this script to a variable

# Determine which configuration file to use
if [ "$#" -lt 1 ]; then
    config_file="${DIR}/config.txt"
else
    config_file=$1
fi

# Check that the configuration file exists
if [ ! -f ${config_file} ]; then
    echo "ERROR: the file ${config_file} does not exist."
    exit
fi

# IMPORT SYSTEM SPECIFIC PATHS AND BATCH ARGUMENTS
TAFFI_PATH=${DIR}/../                                     # It is easier to define a lot of the commands relative to the taffi root
PYTHON_ABS=$(which python)                                # Saves the absolute path to the favored python dist to variable
ORCA_SUBMIT=${DIR}/orca_submit.py
SHELL_SUBMIT=${DIR}/shell_submit.sh
LAMMPS_SUBMIT=${DIR}/lammps_submit.py
LAMMPS_EXE=$( awk -v v="LAMMPS_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
ORCA_EXE=$( awk -v v="ORCA_EXE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FF=$( awk -v v="FF" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
MODULE_STRING=$( awk -v v="MODULE_STRING" -F "*" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
CHARGE=$( awk -v v="CHARGE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
GENS=$( awk -v v="GENS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
BASIS=$( awk -v v="BASIS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
FUNCTIONAL=$( awk -v v="FUNCTIONAL" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )

# FOM SPECIFIC ARGUMENTS
PARAM_GEOOPT_PROCS=$( awk -v v="PARAM_GEOOPT_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_GEOOPT_WT=$( awk -v v="PARAM_GEOOPT_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_GEOOPT_Q=$( awk -v v="PARAM_GEOOPT_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_GEOOPT_SCHED=$( awk -v v="PARAM_GEOOPT_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_GEOOPT_PPN=$( awk -v v="PARAM_GEOOPT_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
PARAM_GEOOPT_SIZE=$( awk -v v="PARAM_GEOOPT_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 

PARAM_BA_PROCS=$( awk -v v="PARAM_BA_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_BA_WT=$( awk -v v="PARAM_BA_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_BA_Q=$( awk -v v="PARAM_BA_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_BA_SCHED=$( awk -v v="PARAM_BA_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_BA_PPN=$( awk -v v="PARAM_BA_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
PARAM_BA_SIZE=$( awk -v v="PARAM_BA_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 

PARAM_D_PROCS=$( awk -v v="PARAM_D_PROCS" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_D_WT=$( awk -v v="PARAM_D_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_D_Q=$( awk -v v="PARAM_D_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_D_SCHED=$( awk -v v="PARAM_D_SCHED" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_D_PPN=$( awk -v v="PARAM_D_PPN" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )  
PARAM_D_SIZE=$( awk -v v="PARAM_D_SIZE" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} ) 

PARAM_FIT_WT=$( awk -v v="PARAM_FIT_WT" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
PARAM_FIT_Q=$( awk -v v="PARAM_FIT_Q" '{ if ( NF > 1 && $1 == v  ) { print $2; } }' ${config_file} )
######################

# Waits for completion of jobs in the queue, expects an 
# array of job ids to be supplied as the first argument
function check_queue () {
    arr=("$@")
    result=('1')
    count=1
    while [ ${#result[@]} -gt 0 ]; do
	tmp=$( python ${DIR}/qstat.py | grep $USER | awk '{ if ( NF > 4 && $5 != "C" ) { print $1 }}')
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

echo -e "\nWorking on final intramolecular parameters..."

# If local parameters are available then use them for the FF argument
if [ -f Params_for_Batch.db ]; then
    FF=$(pwd)/Params_for_Batch.db
fi

##############################
# Run the frag optimizations #
##############################

# Generate Intermolecular Modes Fragments
if [ ! -d Intramolecular_Modes ]; then

    # Assemble shell submit string for intramolecular fragment generation
    # On some systems running simple python scripts on the head node isn't allowed, so it is safer to just always run on compute nodes
    ids_sub=()
    str=""
    if [ ! -z "${MODULE_STRING}" ]; then	    
	   str=${MODULE_STRING}
    fi

    str+="${PYTHON_ABS} ${TAFFI_PATH}FF_functions/frag_gen.py '*.xyz' -FF '${FF}' -gens ${GENS} -q ${CHARGE} -o Intramolecular_Modes"
    tmp=$( echo $( ${SHELL_SUBMIT} "${str}" frag_gen_intra -p 1 -t ${PARAM_FIT_WT} -q ${PARAM_FIT_Q}  ) | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
    readarray -t ids_sub <<< "${tmp}"

    # Wait until jobs complete
    check_queue "${ids_sub[@]}"

    # QC geometry optimize starting fragments      
    cd Intramolecular_Modes
    chmod 777 mass_geoopt.sh

    # Check if there are any intramolecular modes in need of parameterization
    for i in *.xyz; do
      if [ ${i} == "*.xyz" ]; then
          echo "No intramolecular jobs in need of execution. Exiting..."
          cd ..
      else
          break
      fi
    done

    # Generate the geometry optimization files
    for i in *.xyz; do
      name=$( echo ${i} | awk -F '.' '{print $1 }')
      ${PYTHON_ABS} ${TAFFI_PATH}FF_functions/xyz_to_orca.py ${i} -p ${PARAM_GEOOPT_PROCS} -r 0.0 --no_freq -q ${CHARGE} -f ${FUNCTIONAL} -b ${BASIS} -o ${name}_geoopt.in
      mkdir ${name}_geoopt
      mv ${name}_geoopt.in ${name}_geoopt/.
    done

    ids_sub=()
    ids_sub+=$(${PYTHON_ABS} ${ORCA_SUBMIT} -p ${PARAM_GEOOPT_PROCS} -t ${PARAM_GEOOPT_WT} -ppn ${PARAM_GEOOPT_PPN} -q ${PARAM_GEOOPT_Q} -sched ${PARAM_GEOOPT_SCHED} -size ${PARAM_GEOOPT_SIZE} -o optimizations -path_to_exe "${ORCA_EXE}" --silent )$'\n'

    # Wait until jobs complete
    check_queue "${ids_sub[@]}"	

    # Remove job submission files
    rm optimizations.*

    # remove *.gbw files (take up a LOT of memory)
    N_gbw=$( find -name *.gbw | wc -l )
    if [ ${N_gbw} != 0 ]; then
	find -name *.gbw | xargs rm
    fi

    # cd back to gen-* directory
    cd ..

fi

# Check if geoopt jobs are incomplete/need resubmission
cd Intramolecular_Modes
for j in $( seq 10 ); do

    # Check for incompletes/resubmissions
    resub_flag=0
    for i in *.xyz; do
	name=$(echo ${i} | awk -F '[.xyz]' '{ print $1 }' )

	if [ ! -d ${name}_geoopt ]; then
	    echo "ERROR: ${name}_geoopt folder wasn't found in ${i}..."
	    exit
	fi

	if [ ! -f ${name}_geoopt/${name}_geoopt.out ]; then
	    echo "WARNING: ${name}_geoopt/${name}_geoopt.out wasn't found in ${i}..."
	    resub_flag=1
	fi

	flag=$(awk 'BEGIN{flag=0} /THE OPTIMIZATION HAS CONVERGED/ {flag=1}END{print flag}' ${name}_geoopt/${name}_geoopt.out)
	if [ ${flag} -eq 0 ]; then
	    resub_flag=1
	    echo "WARNING: ${name}_geoopt didn't run to completion..."
	fi
    done

    # Attempt resubmission is the resub_flag has been triggered
    if [ ${resub_flag} -eq 1 ]; then
	ids_sub=()
	ids_sub+=$(${PYTHON_ABS} ${ORCA_SUBMIT} -p ${PARAM_GEOOPT_PROCS} -t ${PARAM_GEOOPT_WT} -ppn ${PARAM_GEOOPT_PPN} -q ${PARAM_GEOOPT_Q} -sched ${PARAM_GEOOPT_SCHED} -size ${PARAM_GEOOPT_SIZE} -o optimizations -path_to_exe "${ORCA_EXE}" --silent --resubmit )$'\n'

	# Wait until jobs complete
	check_queue "${ids_sub[@]}"	

    # If nothing is in need of resubmission then break
    else

	# remove *.gbw files (take up a LOT of memory)
	N_gbw=$( find -name *.gbw | wc -l )
	if [ ${N_gbw} != 0 ]; then
	    find -name *.gbw | xargs rm
	fi
	break
    fi

done
cd ..

echo "COMPLETE: geoopt jobs"

######################
# Run the mode scans #
######################

# Check if bond/angle mode scans are needed (flag is set to 1 if so)
ba_flag=0
cd Intramolecular_Modes
for i in *.xyz; do
    name=$(echo ${i} | awk -F '[.xyz]' '{ print $1 }' )

    # Check if the bonds folders are present for all fragments
    for j in ${name}*bond*/; do found=${j}; break; done;

    # Missing folders condition
    if [ ${found} == "${name}*bond*/" ]; then
	ba_flag=1
	break
    fi

    # Check if the angles folders are present for all fragments
    for j in ${name}*angle*/; do found=${j}; break; done;

    # Missing folders condition
    if [ ${found} == "${name}*angle*/" ]; then
	ba_flag=1
	break
    fi	
done
cd ..


# Check if dihedral mode scans are needed (flag is set to 1 if so)
d_flag=0
cd Intramolecular_Modes
for i in *.xyz; do
    name=$(echo ${i} | awk -F '[.xyz]' '{ print $1 }' )

    # Check if the bonds folders are present for all fragments
    for j in ${name}*dihedral*/; do found=${j}; break; done;

    # Missing folders condition
    if [ ${found} == "${name}*dihedral*/" ]; then
	d_flag=1
	break
    fi

done
cd ..

echo "COMPLETE: check on missing modes"

# Generate bond/angle mode scan jobs if needed
ids_sub=()
if [ ${ba_flag} -eq 1 ]; then
    ${PYTHON_ABS} ${TAFFI_PATH}FF_functions/paramgen.py Intramolecular_Modes -p ${PARAM_BA_PROCS} -q ${CHARGE} -gens ${GENS} -modes "bonds angles" -theory dft -f ${FUNCTIONAL} -b ${BASIS} > /dev/null
    cd Intramolecular_Modes
    ids_sub+=$(${PYTHON_ABS} ${ORCA_SUBMIT} -f *in -d bonds_angles -p ${PARAM_BA_PROCS} -t ${PARAM_BA_WT} -q ${PARAM_BA_Q} -sched ${PARAM_BA_SCHED} -ppn ${PARAM_BA_PPN} -size ${PARAM_BA_SIZE} -path_to_exe "${ORCA_EXE}" --silent -o bonds_angles )$'\n'
    cd ..
fi

# Generate dihedral mode scan jobs if needed
if [ ${d_flag} -eq 1 ]; then
    ${PYTHON_ABS} ${TAFFI_PATH}FF_functions/paramgen.py Intramolecular_Modes -p ${PARAM_D_PROCS} -q ${CHARGE} -gens ${GENS} -modes "dihedrals" -theory dft -f ${FUNCTIONAL} -b ${BASIS} -d_step 10.0 --scan > /dev/null
    cd Intramolecular_Modes
    ids_sub+=$(${PYTHON_ABS} ${ORCA_SUBMIT} -f *in -d dihedrals -p ${PARAM_D_PROCS} -t ${PARAM_D_WT} -q ${PARAM_D_Q} -sched ${PARAM_D_SCHED} -ppn ${PARAM_D_PPN} -size ${PARAM_D_SIZE} -path_to_exe "${ORCA_EXE}" --silent -o dihedrals)$'\n'
    cd ..
fi

# Wait until jobs complete
check_queue "${ids_sub[@]}"	

echo "COMPLETE: QC jobs for missing modes"

# Remove old test and REDO folders (sometimes the result of an interuption mid-cycle)
cd Intramolecular_Modes
if [ -d test ]; then 
    rm -r test
fi
N_redo=$( find -name REDO | wc -l )
if [ ${N_redo} != 0 ]; then
    find -name REDO | xargs rm -r 
fi
N_restart=$( find -name RESTART | wc -l )
if [ ${N_restart} != 0 ]; then
    find -name RESTART | xargs rm -r 
fi

cd ..

# Check for bond/angle resubmissions
for i in $( seq 10 ); do

    # If final_params data is already in place then avoid the resubmissions
    if [ -d Intramolecular_Modes/final_params -a -f Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db ]; then
	break
    fi

    # Check for needed resubmissions using the extract_intramolecular_params.py script running in --find_restarts mode.
    ${PYTHON_ABS} ${TAFFI_PATH}FF_functions/extract_intramolecular_params.py -f Intramolecular_Modes -modes "bonds angles harm_dihedrals" -FF "${FF}" -gens ${GENS} --find_restarts -lammps_exe "${LAMMPS_EXE}" -o test >> tmp.txt
    flag=$(awk 'BEGIN{flag=0} /Generating the input files for a reoptimization based on the lowest energy/ {flag=1}END{print flag}' tmp.txt)
    rm -r Intramolecular_Modes/test
    rm tmp.txt

    # Handle flexible dihedral resubmission
    D_flag=$( ${PYTHON_ABS} ${TAFFI_PATH}FF_functions/dihedral_restart.py Intramolecular_Modes --verbose | wc -l )

    if [ ${D_flag} -ne 0 ]; then
	echo "Some flexible dihedral scans did not complete. Resuming..."

	# Submit restarted dihedral scans
	cd Intramolecular_Modes	    
	D_ids_sub=()
	tmp=$(echo $(${PYTHON_ABS} ${ORCA_SUBMIT} -f "*in" -d RESTART -p ${PARAM_D_PROCS} -t ${PARAM_D_WT} -sched ${PARAM_D_SCHED} -q ${PARAM_D_Q} -ppn ${PARAM_D_PPN} -size ${PARAM_D_SIZE} -o dihedrals_restart -path_to_exe "${ORCA_EXE}" --silent) \
	      | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
	readarray -t D_ids_sub <<< "${tmp}"
	cd ..	    
    fi

    if [ ${flag} -eq 1 ]; then
	echo "Some bond/angle jobs failed to converge. Resubmitting..."
	ids_sub=()

	# Run geometry optimizations of the unconverged modes
	cd Intramolecular_Modes
	tmp=$(echo $(${PYTHON_ABS} ${ORCA_SUBMIT} -f "*in" -d REDO -p ${PARAM_GEOOPT_PROCS} -t ${PARAM_GEOOPT_WT} -sched ${PARAM_GEOOPT_SCHED} -q ${PARAM_GEOOPT_Q} -ppn ${PARAM_GEOOPT_PPN} -size ${PARAM_GEOOPT_SIZE} -o ba_geo_resub -path_to_exe "${ORCA_EXE}" --silent) \
	      | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
	readarray -t ids_sub <<< "${tmp}"
	cd ..

	# Wait until jobs complete
	check_queue "${ids_sub[@]}"	

	# Generate and run new mode scans
	${PYTHON_ABS} ${TAFFI_PATH}FF_functions/restart_scans.py Intramolecular_Modes -f ${BASIS} -b ${FUNCTIONAL} -gens ${GENS}
	cd Intramolecular_Modes
	tmp=$(echo $(${PYTHON_ABS} ${ORCA_SUBMIT} -f "*in" -p ${PARAM_BA_PROCS} -t ${PARAM_BA_WT} -sched ${PARAM_BA_SCHED} -q ${PARAM_BA_Q} -size ${PARAM_BA_SIZE} -ppn ${PARAM_BA_PPN} -o ba_scan_resub -path_to_exe "${ORCA_EXE}" --silent) \
	      | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
	readarray -t ids_sub <<< "${tmp}"
	cd ..

	# Wait until jobs complete
	check_queue "${ids_sub[@]}"	

	# Remove REDO folder for a fresh start on the next cycle in case something went wrong
	N_redo=$( find -name REDO | wc -l )
	if [ ${N_redo} != 0 ]; then
	    find -name REDO | xargs rm -r 
	fi

    fi

    # Run the dihedral post processing if jobs were submitted.
    if [ ${D_flag} -ne 0 ]; then

	# Wait until jobs complete
	check_queue "${D_ids_sub[@]}"	

	# Run post processing
	${PYTHON_ABS} ${TAFFI_PATH}FF_functions/dihedral_restart.py Intramolecular_Modes --post

	# Remove RESTART folder for a fresh start on the next cycle in case something went wrong
	N_restart=$( find -name RESTART | wc -l )
	if [ ${N_restart} != 0 ]; then
	    find -name RESTART | xargs rm -r 
	fi	    
    fi

    # Break out of the restart loop if no resubmissions occured
    if [ ${D_flag} -eq 0 -a ${flag} -ne 1 ]; then	    
	break
    fi

done

echo "COMPLETE: QC data for intramolecular mode fit"

# If there was an aborted attempt to fit then remove the old folder
if [ -d Intramolecular_Modes/final_params -a ! -f Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db ]; then
    echo "Removing old final_params folder..."
    rm -r Intramolecular_Modes/final_params
fi

# Fit final params
if [ ! -d Intramolecular_Modes/final_params ]; then

    # remove *.gbw files (take up a LOT of memory)
    N_gbw=$( find -name *.gbw | wc -l )
    if [ ${N_gbw} != 0 ]; then
	find -name *.gbw | xargs rm
    fi

    # Assemble shell submit string
    ids_sub=()
    str=""
    if [ ! -z "${MODULE_STRING}" ]; then	    
	str=${MODULE_STRING}
    fi
    str+="${PYTHON_ABS} ${TAFFI_PATH}FF_functions/extract_intramolecular_params.py -f Intramolecular_Modes -o final_params -FF '${FF}' -min_cycles 5 -gens ${GENS} -max_cycles 10 -N_sc 3 --mixing_rule wh -lammps_exe ${LAMMPS_EXE}" 
    tmp=$( echo $( ${SHELL_SUBMIT} "${str}" final_params -p 1 -t ${PARAM_FIT_WT} -q ${PARAM_FIT_Q} ) | awk '{ for(j=1;j<=NF;j++) { if ( $j != "Submitted" && $j != "batch" && $j != "job" ) {print $j } } }')
    readarray -t ids_sub <<< "${tmp}"

    # Wait until jobs complete
    check_queue "${ids_sub[@]}"	

fi

# Check the results of the fit
if [ ! -d Intramolecular_Modes/final_params -o ! -f Intramolecular_Modes/final_params/Intramolecular_Modes-DFT.db ]; then
    echo "An error occured during the intramolecular mode fit."
    exit
fi

echo "COMPLETE: intramolecular mode fit"

exit
