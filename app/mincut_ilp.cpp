/******************************************************************************
 * mincut_ilp.cpp
 * *
 * Integer linear programs solving the hypercut minimum cut problem near-exactly
 * by using floating point relaxation.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <iostream>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/solvers/ilp.h"
#include "lib/utils/definitions.h"
#include "lib/io/mt_kahypar_io.h"
// Gurobi headers
#include "gurobi_c++.h"
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
    ILPConfig config;

    int parseReturnCode = ParseParams::parse_parameters_mincut_ilp(argc, argv, config);

    if (parseReturnCode)
        return 1;

    // Initialize the timer
    timer t;

    // Read the hypergraph from the given file
    mt_kahypar_hypergraph_t hypergraphWrapper = MtKaHyParIO::read_hypergraph_from_file(config);

    // Cast the hypergraph wrapper to a static hypergraph
    StaticHypergraph &hypergraph = mt_kahypar::utils::cast<StaticHypergraph>(hypergraphWrapper);

    // Extract the number of nodes and edges
    EdgeIndex inputNumEdges = hypergraph.initialNumEdges() - hypergraph.numRemovedHyperedges();
    NodeIndex inputNumNodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();

    if(config.unweighted) {
        for (EdgeID edgeID : hypergraph.edges())
        if (hypergraph.edgeIsEnabled(edgeID))
        hypergraph.setEdgeWeight(edgeID, 1);
    }

    std::cout << "seed \t\t\t\t\t" << config.seed << std::endl;
    std::cout << "unweighted \t\t\t\t" << config.unweighted << std::endl;
    std::cout << "ilp_timeout \t\t\t\t" << config.ilpTimeout << std::endl;
    std::cout << "ilp_mode \t\t\t\t" << config.ilpMode << std::endl;
    std::cout << "initial_num_edges  \t\t\t" << inputNumEdges << std::endl;
    std::cout << "initial_num_nodes  \t\t\t" << inputNumNodes << std::endl;
    std::cout << "io_time \t\t\t\t" << t.elapsed() << std::endl;

    // Restart timer for measuring the time of the ILP
    t.restart();

    try
    {
        // Create an environment
        GRBEnv env = GRBEnv(true);
        env.start();
        // Set the log file
        env.set("LogFile", "gurobi.log");

        // Initialize the ILP solver
        ILPMincut mincutILPSolver(inputNumNodes, inputNumEdges, env, config.seed, config.ilpTimeout, config.ilpMode, config.numThreads);
        // Add the node variables and constraints to the ILP
        mincutILPSolver.add_node_variables_and_constraints(hypergraph);
        // Add the edge variables and constraints to the ILP
        mincutILPSolver.add_edge_variables_and_constraints(hypergraph);
        // Set the objective function of the ILP
        mincutILPSolver.set_objective();
        // Optimize the ILP
        mincutILPSolver.optimize();
        // Print the results of the ILP
        CutValue minEdgeCut = mincutILPSolver.get_result(hypergraph);
        std::cout << "final_mincut_value \t\t" << minEdgeCut << std::endl;
    }
    catch (GRBException e)
    {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    }

    // Print the total time
    std::cout << "total_computing_time \t\t" << t.elapsed() << std::endl;

    // Free the hypergraph wrapper
    mt_kahypar_free_hypergraph(hypergraphWrapper);

    return 0;
}