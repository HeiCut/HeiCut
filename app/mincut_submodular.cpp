/******************************************************************************
 * mincut_submodular.cpp
 * *
 * Solver that finds the minimum cut of a hypergraph via submodular optimization.
 * More specifically, it implements the algorithms of Klimmek and Wagner, Mak
 * and Wong and Queyranne. Note that one instance of the SubmodularMincut class
 * can only be used for one hypergraph (and its contracted versions).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/utils/definitions.h"
#include "lib/utils/random.h"
#include "lib/solvers/submodular.h"
#include "lib/io/mt_kahypar_io.h"
#include "lib/utils/output.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"
#include "mt-kahypar/utils/cast.h"
// KaHIP headers
#include "kahip/timer.h"

int main(int argc, char *argv[])
{
#ifdef NDEBUG
    std::cout << "build_mode \t\t\t\tRELEASE" << std::endl;
#else
    std::cout << "build_mode \t\t\t\tDEBUG" << std::endl;
#endif

    // Stores the config
    SubmodularConfig config;

    int parseReturnCode = ParseParams::parse_parameters_mincut_submodular(argc, argv, config);

    if (parseReturnCode)
        return 1;

    // Initialize the timer
    timer t;

    // Set the seed for the random number generator
    RandomFunctions::set_seed(config.seed);

    // Read the hypergraph from the given file
    mt_kahypar_hypergraph_t hypergraphWrapper = MtKaHyParIO::read_hypergraph_from_file(config);

    // Cast the hypergraph wrapper to a dynamic hypergraph
    DynamicHypergraph &hypergraph = mt_kahypar::utils::cast<DynamicHypergraph>(hypergraphWrapper);

    // Extract the number of nodes and edges
    EdgeIndex inputNumEdges = hypergraph.initialNumEdges() - hypergraph.numRemovedHyperedges();
    NodeIndex inputNumNodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();

    // Determine if the hypergraph has weighted edges
    bool hasWeightedEdges = false;
    for (EdgeID edgeID : hypergraph.edges())
        if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeWeight(edgeID) != 1)
        {
            hasWeightedEdges = true;
            break;
        }

    std::cout << "seed \t\t\t\t\t" << config.seed << std::endl;
    std::cout << "ordering_type \t\t\t\t" << config.orderingType << std::endl;
    std::cout << "ordering_mode \t\t\t\t" << config.orderingMode << std::endl;
    std::cout << "has_weighted_edges \t\t\t" << (hasWeightedEdges ? "true" : "false") << std::endl;
    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;
    // Initialize the total time
    double totalComputingTime = 0;

    // Restart timer for measuring the time of the preprocessing
    t.restart();

    // Initialize the submodular mincut solver
    SubmodularMincut mincutSubmodularSolver(hypergraph.initialNumNodes(), inputNumEdges, config.orderingType, config.orderingMode, hasWeightedEdges, config.numThreads);

    // Print the time of the preprocessing
    double preprocessingTime = t.elapsed();
    std::cout << "===================================================" << std::endl;
    std::cout << "preprocessing_time \t\t\t" << preprocessingTime << std::endl;
    std::cout << "===================================================" << std::endl;
    totalComputingTime += preprocessingTime;

    // Restart the timer for measuring the time of computing the minimum cut for the hypergraph
    t.restart();

    // Solve the minimum cut problem via submodular optimization
    SubmodularMincutResult result = mincutSubmodularSolver.solve(hypergraph);

    // Print the time of computing the minimum cut for the hypergraph
    double solvingTime = t.elapsed();
    std::cout << "solving_time \t\t\t\t" << solvingTime << std::endl;
    std::cout << "solving_num_iterations \t\t\t" << result.numIterations << std::endl;
    std::cout << "solving_contractions_per_it \t\t" << result.meanContractionsPerIteration << std::endl;
    std::cout << "===================================================" << std::endl;
    totalComputingTime += solvingTime;

    // Print the final minimum edge cut value
    std::cout << "final_mincut_value \t\t\t" << result.minEdgeCut << std::endl;
    // Print the total time
    std::cout << "total_computing_time \t\t\t" << totalComputingTime << std::endl;

    // Free the hypergraph wrapper
    mt_kahypar_free_hypergraph(hypergraphWrapper);

    return 0;
}