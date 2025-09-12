/******************************************************************************
 * kcore_generator.cpp
 * *
 * Generates a (k, 2)-core subgraph of a given hypergraph. The (k, 2)-core is a
 * maximal subgraph in which all nodes have at least unweighted degree k and all
 * hyperedges have at least a size of 2. This is the direct equivalent of k-core
 * definition on normal graphs, which was implemented in VieCut. This algorithm
 * is a port of the VieCut implementation for normal graphs to hypergraphs. The
 * VieCut implementation is based on the algorithm by  Batagelj and Zaversnik.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <fstream>
#include <iostream>
#include <cassert>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/utils/random.h"
#include "lib/io/mt_kahypar_io.h"
#include "lib/decomposition/core_decomposition.h"
#include "lib/solvers/kernelizer.h"
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
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
    KCoreGeneratorConfig config;

    int parseReturnCode = ParseParams::parse_parameters_kcore_generator(argc, argv, config);

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
    NodeIndex inputNumPins = hypergraph.initialNumPins();

    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "initial_num_pins  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;

    // Initialize the core decomposition
    CoreDecomposition coreDecomposition(inputNumNodes, inputNumEdges, inputNumPins);
    // Decompose the hypergraph into its k-core components
    coreDecomposition.decompose_into_cores(hypergraph);

    std::vector<EdgeIndex> targetCores = coreDecomposition.get_target_cores();

    if (targetCores.empty())
    {
        std::cout << "No target cores found." << std::endl;
        return 0;
    }

    std::cout << "smallest_core  \t\t\t\t" << targetCores[0] << std::endl;
    std::cout << "largest_core  \t\t\t\t" << targetCores[targetCores.size() - 1] << std::endl;
    std::cout << "===================================================" << std::endl;

    // Initialize the kernelizer config
    KernelizerConfig kernelizerConfig;
    kernelizerConfig.seed = config.seed;
    kernelizerConfig.numThreads = config.numThreads;
    kernelizerConfig.baseSolver = DEFAULT_BASE_SOLVER;
    kernelizerConfig.ilpTimeout = DEFAULT_TIMEOUT_ILP;
    kernelizerConfig.ilpMode = DEFAULT_ILP_MODE;
    kernelizerConfig.orderingType = DEFAULT_ORDERING_TYPE;
    kernelizerConfig.orderingMode = DEFAULT_ORDERING_MODE;
    kernelizerConfig.pruningMode = DEFAULT_PRUNING_MODE;
    // Deactivate LP to get an exact mincut value
    kernelizerConfig.LPNumIterations = 0;
    kernelizerConfig.LPMode = DEFAULT_LP_MODE;
    kernelizerConfig.LPNumPinsToSample = DEFAULT_LP_NUM_PINS_TO_SAMPLE;
    // Deactivate verbose mode to suppress output
    kernelizerConfig.verbose = false;

    // Initialize the kernelizer
    Kernelizer kernelizer(kernelizerConfig);

    for (EdgeIndex k : targetCores)
    {
        mt_kahypar_hypergraph_t coreHypergraphWrapper = coreDecomposition.create_k_core_graph(k, hypergraph);
        // Make sure that the core hypergraph is a static hypergraph
        assert(coreHypergraphWrapper.type == STATIC_HYPERGRAPH);

        // Cast the core hypergraph wrapper to a static core hypergraph
        auto &coreHypergraph = mt_kahypar::utils::cast<StaticHypergraph>(coreHypergraphWrapper);

        EdgeIndex coreNumEdges = coreHypergraph.initialNumEdges() - coreHypergraph.numRemovedHyperedges();
        NodeIndex coreNumNodes = coreHypergraph.initialNumNodes() - coreHypergraph.numRemovedHypernodes();

        std::cout << "core  \t\t\t\t\t" << k << std::endl;
        std::cout << "num_edges  \t\t\t\t" << coreNumEdges << std::endl;
        std::cout << "num_nodes  \t\t\t\t" << coreNumNodes << std::endl;

        // Compute the minimum cut of the hypergraph using kernelization
        // NB: We need to copy the hypergraph, since the kernelizer modifies/contracts it
        StaticHypergraph copyCoreHypergraph = coreHypergraph.copy();
        KernelizerResult result = kernelizer.compute_mincut(copyCoreHypergraph, coreNumNodes, coreNumEdges);

        if (result.minEdgeCut < result.naiveEstimate)
        {
            std::cout << "Mincut value (" << result.minEdgeCut << ") is smaller than naive estimate (" << result.naiveEstimate << ")." << std::endl;

            // Create a string stream buffer
            std::ostringstream buffer;
            std::string outputFileName = config.outputFileName;
            outputFileName += ".core_" + std::to_string(k);
            // Open output file
            std::ofstream f(outputFileName);

            // Iterate over all hyperedges
            for (EdgeID edgeID : coreHypergraph.edges())
            {
                // Ignore disabled hyperedges
                if (!coreHypergraph.edgeIsEnabled(edgeID))
                    continue;

                buffer << coreHypergraph.edgeWeight(edgeID) << " ";

                // Iterate over the pins of the hyperedge
                for (NodeID pinID : coreHypergraph.pins(edgeID))
                {
                    // Ignore disabled pins
                    if (!coreHypergraph.nodeIsEnabled(pinID))
                        continue;

                    buffer << pinID + 1 << " ";
                }

                buffer << std::endl;
            }

            // Iterate over all hypernodes
            for (NodeID nodeID : coreHypergraph.nodes())
            {
                // Ignore disabled nodes
                if (!coreHypergraph.nodeIsEnabled(nodeID))
                    continue;

                buffer << coreHypergraph.nodeWeight(nodeID) << std::endl;
            }

            std::cout << "Writing core hypergraph to: " << outputFileName << std::endl;

            // Write the header to the file
            f << coreNumEdges << " " << coreNumNodes << " 11" << std::endl;
            // Write the hyperedges to the file
            f << buffer.str();
            // Close the file
            f.close();

            // Free the core hypergraph wrapper
            mt_kahypar_free_hypergraph(coreHypergraphWrapper);
            return 0;
        }
        else
        {
            std::cout << "Mincut value (" << result.minEdgeCut << ") is equal to naive estimate (" << result.naiveEstimate << ")." << std::endl;
        }

        std::cout << "===================================================" << std::endl;

        // Free the core hypergraph wrapper
        mt_kahypar_free_hypergraph(coreHypergraphWrapper);
    }

    std::cout << "No core subgraph found where the mincut value is smaller than the naive estimate." << std::endl;

    return 0;
}