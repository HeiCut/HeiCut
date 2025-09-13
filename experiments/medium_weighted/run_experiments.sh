#!/bin/bash
###############################################################################
# run_experiments.sh
# #
# Runs the commands listed every 'all_commands.txt', which are all of the
# commands for the generated experiments.
# #
# Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
###############################################################################

# Paths
PATH_TO_GENERATED_EXPERIMENTS_FOLDER=`cat paths/path_to_generated_experiments_folder.txt`;
PATH_TO_TMP_FILE="./tmp_run.txt";

# Number of jobs used for the parallel execution
NUMBER_OF_JOBS_IN_PARALLEL=8;

# Timeout for each instance
TIMEOUT=8100; # 2 hours and 15 minutes

# Whether to reserve the cores exclusively for each job
RESERVE_CORES=0;

# Get CMD flags
while getopts 'n:r' flag
do
    case "$flag" in
        n)
            # Set the number of jobs
            NUMBER_OF_JOBS_IN_PARALLEL=${OPTARG};
        ;;
        r)
            # Reserve the cores exclusively for each job
            RESERVE_CORES=1
        ;;
        ?)
            echo "script usage: ./run_experiments.sh -n NUMBER_OF_JOBS_IN_PARALLEL -r" >&2;
            exit 1;
        ;;
    esac
done

# Get the number of threads defined when generating the experiments
eval $(grep '^NUM_THREADS=' generate_experiments.sh);

# Total number of logical cores
TOTAL_NUMBER_LOGICAL_CORES=$(nproc);

# Total number of physical cores 
TOTAL_NUMBER_PHYSICIAL_CORES=$(( $TOTAL_NUMBER_LOGICAL_CORES / 2 ));

# Compute the number of physical cores for each job
NUMBER_PHYSICAL_CORES_PER_JOB=$(( $TOTAL_NUMBER_PHYSICIAL_CORES / $NUMBER_OF_JOBS_IN_PARALLEL ));

if (( $NUMBER_OF_JOBS_IN_PARALLEL > $TOTAL_NUMBER_LOGICAL_CORES && $RESERVE_CORES == 1 ))
then
    echo "The number of jobs (${NUMBER_OF_JOBS_IN_PARALLEL}) is greater than the number of logical cores (${TOTAL_NUMBER_LOGICAL_CORES}).";
    exit 1;
fi

if (( NUM_THREADS > 2 * NUMBER_PHYSICAL_CORES_PER_JOB && $RESERVE_CORES == 1 )) 
then
    echo "The number of threads (${NUM_THREADS}) is greater than the number of logical cores per job ($((2 * NUMBER_PHYSICAL_CORES_PER_JOB))).";
    exit 1;
fi

if [ $RESERVE_CORES -eq 1 ]
then
    # Add a replacement to the config file for GNU Parallel to get the slot number (ranging from 0 to NUMBER_OF_JOBS_IN_PARALLEL - 1)
    # Taken from: https://stackoverflow.com/a/18625528
    replacement='--rpl "{%%} \$_=slot()-1"'

    # Check if the replacement is already in the config file
    if ! grep -Fxq -- "$replacement" "$HOME/.parallel/config"; then
        # If not found, append the replacement to the config file
        echo "$replacement" >> "$HOME/.parallel/config"
        echo "Replacement {%%} added to config file of GNU parallel."
    fi;
fi;

# Check if the 'generated/all_commands' folder exists
if [ -d "${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/all_commands" ]
then
    # Clear old temp file
    if [ -f "$PATH_TO_TMP_FILE" ]
    then
        rm ${PATH_TO_TMP_FILE};
    fi
    # Loop over all algorithms
    for algorithm_commands in `find ${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/all_commands/* -iname all_commands.txt`
    do
        # Add command to temp file
        cat "${algorithm_commands}" >> ${PATH_TO_TMP_FILE};
    done
    
    if [ -f "$PATH_TO_TMP_FILE" ]
    then
        # Execute all the commands parallel with the given number of cores 
        if [ $RESERVE_CORES -eq 1 ]
        then
            # Reserve cores exclusively for each job based on the slot number (given by {%%})
            echo "Reserving ${NUMBER_PHYSICAL_CORES_PER_JOB} (= ${TOTAL_NUMBER_PHYSICIAL_CORES} // ${NUMBER_OF_JOBS_IN_PARALLEL}) physical cores and thus $((2 * NUMBER_PHYSICAL_CORES_PER_JOB)) logical cores per job.";
            cat ${PATH_TO_TMP_FILE} | parallel -j${NUMBER_OF_JOBS_IN_PARALLEL} --timeout ${TIMEOUT} \
            'physical_core_start=$(( {%%} * '${NUMBER_PHYSICAL_CORES_PER_JOB}' )); 
            physical_core_end=$(( physical_core_start + '${NUMBER_PHYSICAL_CORES_PER_JOB}' - 1 )); 
            taskset -c $(for ((i=$physical_core_start; i<=$physical_core_end; i++)); do echo -n "$i,$((i+'${TOTAL_NUMBER_PHYSICIAL_CORES}')),"; done | sed "s/,$//") sh -c {}';
        else
            cat ${PATH_TO_TMP_FILE} | parallel -j${NUMBER_OF_JOBS_IN_PARALLEL} --timeout ${TIMEOUT} sh;
        fi;
        
        # Clear temp file
        rm ${PATH_TO_TMP_FILE};
    fi;
else
    exit 1;
fi
