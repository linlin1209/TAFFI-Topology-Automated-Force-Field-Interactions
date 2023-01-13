#!/bin/bash

# Help message
if [ $# -lt 3 ]; then
    echo ""
    echo "The script expects three arguments."
    echo ""
    echo "    1: a quoted string containing a basename with a wildcard (*) character. (e.g., 'test.*.lammpstrj'"
    echo "    2: starting index (e.g., 0)"
    echo "    3: ending index (inclusive) (e.g., 10)"
    echo ""
    echo "For example, supplying 'test.*.lammpstrj' 0 10, would combine trajectory files test.0.lammpstrj through test.10.lammpstrj"
    echo "If --thermo is supplied as an optional fourth argument, then the script will assume the files are lammps thermo files."
    exit
fi

# Check for thermo flag
if [ $# -gt 3 -a "$4" == "--thermo" ]; then
    thermo_flag=1
else
    thermo_flag=0
fi

# Save the root name to a variable
root=${1}

# Combine lammps trajectories
if [ ${thermo_flag} == 0 ]; then

    counter=0
    for i in `seq ${2} ${3}`; do
	file=$( echo "$1" | sed -r "s/\*/${i}/g")

	# Check that the file exists
	if [ ! -f ${file} ]; then
	    echo "ERROR: could not find ${file}. Exiting..."
	    exit
	fi

	if [ ${counter} -eq 0 ]; then
	    cat ${file} > combined.lammpstrj
	else
	    awk 'BEGIN{count=0}{ if ( NF == 2 && $2 == "TIMESTEP" ) { count++} if ( count > 1 ) {print $0}}' ${file} >> combined.lammpstrj
	fi
	counter=$(( ${counter} + 1 ))
    done

# Combine thermo files
elif [ ${thermo_flag} == 1 ]; then

    counter=0
    for i in `seq ${2} ${3}`; do
	file=$( echo "$1" | sed -r "s/\*/${i}/g")

	# Check that the file exists
	if [ ! -f ${file} ]; then
	    echo "ERROR: could not find ${file}. Exiting..."
	    exit
	fi

	if [ ${counter} -eq 0 ]; then
	    cat ${file} > combined.thermo
	else
	    awk '{ if ( NR > 2 ) {print $0}}' ${file} >> combined.thermo
	fi
	counter=$(( ${counter} + 1 ))
    done
fi