/******************************************************************************
 * mincut_trimmer.cpp
 * *
 * Exact hypergraph minimum cut algorithm using the k-trimmer of Chekuri and Xu.
 * Note that this algorithm currently only works for unweighted hypergraphs.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
#include <cassert>
#include <vector>
#include <algorithm>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/orderer/orderer.h"
#include "lib/trimmer/trimmer.h"
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
    TrimmerConfig config;

    int parseReturnCode = ParseParams::parse_parameters_mincut_trimmer(argc, argv, config);

    if (parseReturnCode)
        return 1;

    // Initialize the timer
    timer t;

    // Set the seed for the random number generator
    RandomFunctions::set_seed(config.seed);

    // Read the hypergraph from the given file
    mt_kahypar_hypergraph_t hypergraphWrapper = MtKaHyParIO::read_hypergraph_from_file(config);

    // Cast the hypergraph wrapper to a static hypergraph
    StaticHypergraph &hypergraph = mt_kahypar::utils::cast<StaticHypergraph>(hypergraphWrapper);

    // Extract the number of nodes and edges
    EdgeIndex inputNumEdges = hypergraph.initialNumEdges() - hypergraph.numRemovedHyperedges();
    NodeIndex inputNumNodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();

    // The algorithm only works for unweighted hypergraphs
    bool hasWeightedEdges = false;
    for (EdgeID edgeID : hypergraph.edges())
        if (hypergraph.edgeIsEnabled(edgeID))
        {
            hypergraph.setEdgeWeight(edgeID, 1);
        }

    std::cout << "seed \t\t\t\t\t" << config.seed << std::endl;
    std::cout << "ordering_type \t\t\t\t" << config.orderingType << std::endl;
    std::cout << "ordering_mode \t\t\t\t" << config.orderingMode << std::endl;
    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;

    // Initialize the total time
    double totalComputingTime = 0;

    // Restart timer for measuring the time of the preprocessing
    t.restart();

    // Initialize the ordering of the nodes
    std::vector<NodeID> nodeOrdering(inputNumNodes);
    // Initialize the head ordering of the edges
    std::vector<EdgeID> edgeHeadOrdering(inputNumEdges);
    // Initialize the heads of the edges
    std::vector<NodeID> edgeHead(inputNumEdges);

    // Perform the MA-ordering of the nodes on the input hypergraph
    Orderer<StaticHypergraph, EdgeWeight> inputOrderer(hypergraph.initialNumNodes(), inputNumEdges, OrderingType::MA, hasWeightedEdges, RandomFunctions::get_random_engine());
    inputOrderer.compute_ordering(hypergraph, inputNumNodes, &nodeOrdering, &edgeHeadOrdering, &edgeHead);

    // Initialize the trimmer to generate a k-trimmed hypergraph for HIGHEST_QUALITY preset (= dynamic hypergraph)
    Trimmer trimmer(hypergraph.initialNumNodes(), inputNumNodes, inputNumEdges, hypergraph, nodeOrdering, edgeHeadOrdering, edgeHead, HIGHEST_QUALITY);
    // Build the backward edges
    trimmer.build_backward_edges();
    // trimmer.print_backward_edges();

    // Start with k = 2 and perform exponential search on k
    TrimmerValue k = 2;
    // Initialize the minimum edge cut value
    CutValue minEdgeCut;
    // Initialize the submodular mincut solver
    SubmodularMincut mincutSubmodularSolver(hypergraph.initialNumNodes(), inputNumEdges, config.orderingType, config.orderingMode, hasWeightedEdges, config.numThreads);

    // Print the time of the preprocessing
    double preprocessingTime = t.elapsed();
    std::cout << "===================================================" << std::endl;
    std::cout << "preprocessing_time \t\t\t" << preprocessingTime << std::endl;
    std::cout << "===================================================" << std::endl;
    // Increase the total time
    totalComputingTime += preprocessingTime;

    // Continue doubling k until (k > minEdgeCut) is satisfied (see end of while loop)
    while (true)
    {
        // Restart the timer for measuring the time of computing the minimum cut for the current k-trimmed hypergraph
        t.restart();

        // Create the the k-trimmed certificate of the hypergraph
        mt_kahypar_hypergraph_t trimmedHypergraphWrapper = trimmer.create_k_trimmed_certificate(k);
        // Make sure that the trimmed hypergraph is a dynamic hypergraph
        assert(trimmedHypergraphWrapper.type == DYNAMIC_HYPERGRAPH);

        // Cast the trimmed hypergraph wrapper to a dynamic trimmed hypergraph
        DynamicHypergraph &trimmedHypergraph = mt_kahypar::utils::cast<DynamicHypergraph>(trimmedHypergraphWrapper);

        // Solve the minimum cut problem via submodular optimization
        SubmodularMincutResult result = mincutSubmodularSolver.solve(trimmedHypergraph);
        minEdgeCut = result.minEdgeCut;

        // Free the trimmed hypergraph wrapper
        mt_kahypar_free_hypergraph(trimmedHypergraphWrapper);

        // Print the minimum edge cut value for the current k-trimmed hypergraph
        std::cout << k << "_trimmed_mincut_value \t\t\t" << minEdgeCut << std::endl;
        std::cout << k << "_trimmed_num_iterations \t\t" << result.numIterations << std::endl;
        std::cout << k << "_trimmed_contractions_per_it \t\t" <<  result.meanContractionsPerIteration << std::endl;
        // Print the time of computing the minimum cut for the current k-trimmed hypergraph
        double kTrimmedTime = t.elapsed();
        std::cout << k << "_trimmed_time \t\t\t\t" << kTrimmedTime << std::endl;
        std::cout << "===================================================" << std::endl;
        totalComputingTime += kTrimmedTime;

        // Check if we found the minimum edge cut
        if (k > minEdgeCut)
            break;

        // Otherwise, double k
        k *= 2;
    }

    // Print the final minimum edge cut value
    std::cout << "final_mincut_value \t\t\t" << minEdgeCut << std::endl;
    // Print the total time
    std::cout << "total_computing_time \t\t\t" << totalComputingTime << std::endl;

    // Free the hypergraph wrapper
    mt_kahypar_free_hypergraph(hypergraphWrapper);

    return 0;
}