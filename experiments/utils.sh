#!/bin/bash
###############################################################################
# utils.sh
# #
# Utility methods that are used by the other scripts.
# #
# Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
###############################################################################

# Utlility method to check if a value exists in a list
exists_in_list() {
    passed_list=$1
    passed_value=$2
    echo $passed_list | tr " " '\n' | grep -F -q -x "$passed_value";
}

# Utlility method to create the folders for a passed algorithm and path
create_folders_for_algorithm_and_path() {
    passed_path=$1
    passed_algorithm=$2
    
    # Create folder for the passed algorithm and path
    mkdir -p ${passed_path}/${passed_algorithm};
    # Create folder containing all commands for the passed algorithm
    mkdir -p ${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/all_commands/${passed_algorithm};
}

# Utlility method to store the command of the passed experiment
store_command() {
    passed_path=$1
    passed_algorithm=$2
    passed_command=$3
    passed_memory_limit=$4

    # Store the command running the passed algorithm under the path
    if [ -z "$passed_memory_limit" ]; then
    	 echo "/usr/bin/time -v ${passed_command} 1> ${passed_path}/${passed_algorithm}/result.txt 2> ${passed_path}/${passed_algorithm}/stats.txt" > "${passed_path}/${passed_algorithm}/command.sh";
    else
         echo "ulimit -v ${passed_memory_limit}; /usr/bin/time -v ${passed_command} 1> ${passed_path}/${passed_algorithm}/result.txt 2> ${passed_path}/${passed_algorithm}/stats.txt" > "${passed_path}/${passed_algorithm}/command.sh";
    fi;

    # Store the path to command running the passed algorithm in the file with all commands
    echo "${passed_path}/${passed_algorithm}/command.sh" >> "${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/all_commands/${passed_algorithm}/all_commands.txt";
    # Make command executable
    chmod +x "${passed_path}/${passed_algorithm}/command.sh";
}
