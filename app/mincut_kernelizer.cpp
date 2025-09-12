/******************************************************************************
 * mincut_kernelizer.cpp
 * *
 * Computes the (estimated) mincut of a hypergraph. For this, the kernel of the
 * hypergraph is built by coarsening the hypergraph using the label propagation
 * and pruning rules. Afterwards, the mincut of the coarsened hypergraph is
 * computed using a submodular solver or an ILP (with floating point relaxation).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/utils/definitions.h"
#include "lib/io/mt_kahypar_io.h"
#include "lib/solvers/kernelizer.h"
#include "lib/utils/random.h"
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
    KernelizerConfig config;

    int parseReturnCode = ParseParams::parse_parameters_kernelizer(argc, argv, config);

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

    std::cout << "seed \t\t\t\t\t" << config.seed << std::endl;
    std::cout << "base_solver \t\t\t\t" << config.baseSolver << std::endl;
    if (config.baseSolver == BaseSolver::ILP)
    {
        std::cout << "ilp_timeout \t\t\t\t" << config.ilpTimeout << std::endl;
        std::cout << "ilp_mode \t\t\t\t" << config.ilpMode << std::endl;
    }
    else
    {
        std::cout << "ordering_type \t\t\t\t" << config.orderingType << std::endl;
        std::cout << "ordering_mode \t\t\t\t" << config.orderingMode << std::endl;
    }
    std::cout << "pruning_mode \t\t\t\t" << config.pruningMode << std::endl;
    std::cout << "lp_num_iterations \t\t\t" << config.LPNumIterations << std::endl;
    std::cout << "lp_mode \t\t\t\t" << config.LPMode << std::endl;
    if (config.LPMode == LabelPropagationMode::PROBABILISTIC)
        std::cout << "lp_num_pins_to_sample \t\t\t" << config.LPNumPinsToSample << std::endl;
    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;

    // Initialize the kernelizer
    Kernelizer kernelizer(config);

    // Compute the minimum cut of the hypergraph using kernelization
    KernelizerResult result = kernelizer.compute_mincut(hypergraph, inputNumNodes, inputNumEdges);

    std::cout << "===================================================" << std::endl;
    // Print the final minimum edge cut value
    std::cout << "final_mincut_value \t\t\t" << result.minEdgeCut << std::endl;
    // Print the total time
    std::cout << "total_computing_time \t\t\t" << result.time << std::endl;

    // Free the hypergraph wrapper
    mt_kahypar_free_hypergraph(hypergraphWrapper);

    return 0;
}