#!/bin/bash
###############################################################################
# extract_results.sh
# #
# Extracts the results of the experiments and writes them to the
# 'all_results' folder.
# #
# Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
###############################################################################

# Paths
PATH_TO_GENERATED_EXPERIMENTS_FOLDER=`cat paths/path_to_generated_experiments_folder.txt`;
PATH_TO_ALL_RESULTS_FOLDER="${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/all_results";

# Check if the 'generated' folder exists
if [ -d ${PATH_TO_GENERATED_EXPERIMENTS_FOLDER} ]
then
    # Clear old results folder
    if [ -d "$PATH_TO_ALL_RESULTS_FOLDER" ]
    then
        rm -r ${PATH_TO_ALL_RESULTS_FOLDER};
    fi
    
    # Find every 'result.txt' file nested inside the 'generated' folder
    for result in `find ${PATH_TO_GENERATED_EXPERIMENTS_FOLDER}/* -iwholename */result.txt`
    do
        # Get the directory of the found 'result.txt' file
        directory=`dirname $result`;
        algorithm=`basename $directory`;
        
        # Initialize values
        final_mincut_value=0;
        total_computing_time=0;
        peak_memory=0;
        
        # Get the parameters
        parameters=`dirname $directory | awk -F "/jobs/" '{print $2}' | tr / ";" | sed "s/seed_//g"`;
        hypergraph="${parameters%%;*}";
        
        # Create new results folder for the given algorithm if it does not exist yet
        if [ ! -d "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}" ]
        then
            mkdir -p ${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm};
        fi
        
        total_computing_time=`grep "^total_computing_time" ${result}  | awk '{print $NF}'`;
        peak_memory=`grep "Maximum resident set size (kbytes):" ${directory}/stats.txt | awk '{print $NF}'`;
        
        if [[ $algorithm == "hypercactus_generator"* ]]
        then
            ######################################## hypercactus_generator ########################################
            
            # Write the results to the "all_hypercactus_results.csv" file
            initial_num_edges=`grep "^initial_num_edges" ${result}  | awk '{print $NF}'`;
            initial_num_nodes=`grep "^initial_num_nodes" ${result}  | awk '{print $NF}'`;
            kernel_num_edges=`grep "^kernel_num_edges" ${result}  | awk '{print $NF}'`;
            kernel_num_nodes=`grep "^kernel_num_nodes" ${result}  | awk '{print $NF}'`;
            hypercactus_num_edges=`grep "^hypercactus_num_edges" ${result}  | awk '{print $NF}'`;
            hypercactus_num_nodes=`grep "^hypercactus_num_nodes" ${result}  | awk '{print $NF}'`;
            echo "${initial_num_edges};${initial_num_nodes};${kernel_num_edges};${kernel_num_nodes};${hypercactus_num_edges};${hypercactus_num_nodes};${total_computing_time};${peak_memory};${parameters};${algorithm}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}/all_hypercactus_results.csv";
        else
            final_mincut_value=`grep "^final_mincut_value" ${result}  | awk '{print $NF}'`;
            
            if [[ $algorithm == "trimmer"* ]]
            then
                ############################################### trimmer ###############################################
                
                is_exact=$(grep -q "^final_mincut_value" ${result} && echo 1 || echo 0);
                
            elif [[ $algorithm == "submodular"* ]]
            then
                ######################################## submodular(_parallel) ########################################
                
                is_exact=$(grep -q "^final_mincut_value" ${result} && echo 1 || echo 0);
                
                # Write the iterations to the "all_iterations.csv" file
                initial_num_nodes=`grep "^initial_num_nodes" ${result}  | awk '{print $NF}'`;
                solving_num_iterations=`grep "^solving_num_iterations" ${result}  | awk '{print $NF}'`;
                solving_contractions_per_it=`grep "^solving_contractions_per_it" ${result}  | awk '{print $NF}'`;
                echo "${initial_num_nodes};${solving_num_iterations};${solving_contractions_per_it};${parameters};${algorithm}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}/all_iterations.csv";
                
            elif [[ $algorithm == "ilp"* ]]
            then
                ############################################ ilp(_parallel) ############################################
                
                is_exact=$(grep -q "ILP: Optimal solution found." ${result} && echo 1 || echo 0);
                
            elif [[ $algorithm == "maxsat"* ]]
            then
                ############################################ maxsat(_parallel) ############################################
                
                is_exact=$(grep -q "MaxSAT: Optimal solution found." ${result} && echo 1 || echo 0);
                
                # Write the size of the clauses to the "all_clauses.csv" file
                initial_num_edges=`grep "^initial_num_edges" ${result}  | awk '{print $NF}'`;
                initial_num_nodes=`grep "^initial_num_nodes" ${result}  | awk '{print $NF}'`;
                num_hard_clauses=`grep "^num_hard_clauses" ${result}  | awk '{print $NF}'`;
                num_soft_clauses=`grep "^num_soft_clauses" ${result}  | awk '{print $NF}'`;
                echo "${initial_num_edges};${initial_num_nodes};${num_hard_clauses};${num_soft_clauses};${parameters};${algorithm}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}/all_clauses.csv";

            elif [[ $algorithm == "kernelizer"* ]]
            then
                ######################################## kernelizer(_parallel) ########################################
                
                # Write the naive mincut to the "all_naive_mincuts.csv" file (only if the algorithm is kernelizer_IT0)
                if echo "${algorithm}" | grep -q "_IT0";
                then
                    naive_mincut_value=`grep "^naive_mincut_value" ${result}  | awk '{print $NF}'`;
                    echo "${naive_mincut_value};${hypergraph}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/all_naive_mincuts.csv";
                fi
                
                base_solver=`grep "^base_solver[[:space:]]" ${result}  | awk '{print $NF}'`;
                base_solver_time=`grep "^base_solver_time" ${result}  | awk '{print $NF}'`;
                
                if [[ $base_solver == "ilp" ]]
                then
                    base_solver_exact_output="ILP: Optimal solution found.";
                else
                    base_solver_exact_output="^base_solver_mincut_value";
                fi
                
                is_exact=$(grep -q -e "${base_solver_exact_output}" -e "no base solver necessary." ${result} && \
                    (echo ${algorithm} | grep -q "_IT0") && \
                echo 1 || echo 0);
                
                # Write the reductions to the "all_reductions.csv" file
                initial_num_edges=`grep "^initial_num_edges" ${result}  | awk '{print $NF}'`;
                initial_num_nodes=`grep "^initial_num_nodes" ${result}  | awk '{print $NF}'`;
                reductions=`grep -E "^abs|coarsening_time|_pruning_time" ${result} | awk '{print $NF}' | paste -sd ';'`
                echo "${initial_num_edges};${initial_num_nodes};${reductions};${base_solver_time};${parameters};${algorithm}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}/all_reductions.csv";
            fi
            
            # Write to the 'all_results.csv' file
            echo "${final_mincut_value};${total_computing_time};${peak_memory};${is_exact};${parameters};${algorithm}" >> "${PATH_TO_ALL_RESULTS_FOLDER}/${algorithm}/all_results.csv";
        fi
    done
else
    exit 1;
fi