/******************************************************************************
 * ilp.h
 * *
 * ILP that solves the mincut problem on a hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_ILP_H
#define SMHM_ILP_H

#include <vector>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/utils/const.h"
// Gurobi headers
#include "gurobi_c++.h"

class ILPMincut
{
private:
    // Environment of the ILP
    GRBEnv env;
    // Model of the ILP
    GRBModel model;

    // Number of nodes in the hypergraph
    const NodeIndex numNodes;
    // Number of edges in the hypergraph
    const EdgeIndex numEdges;

    // Stores the node variables of the ILP
    std::vector<GRBVar> nodes;
    // Stores the edge variables of the ILP
    std::vector<GRBVar> edges;
    // Stores the minimum variables of the ILP (only if ilpMode = MILP)
    std::vector<GRBVar> minimums;
    // Stores the maximum variables of the ILP (only if ilpMode = MILP)
    std::vector<GRBVar> maximums;

    // Linear expression for the total number of nodes in partition 1
    GRBLinExpr nodeTot = 0;
    //  Linear expression for the total weight of the cut
    GRBLinExpr edgeTot = 0;

    // Timeout of the ILP in seconds
    const double ilpTimeout = DEFAULT_TIMEOUT_ILP;
    // Mode of the ILP
    const ILPMode ilpMode = DEFAULT_ILP_MODE;

public:
    ILPMincut(const NodeIndex numNodes, const EdgeIndex numEdges, const GRBEnv env, const int seed, const double ilpTimeout, const ILPMode ilpMode, const size_t numThreads);

    // Add the node variables and constraints to the ILP
    void add_node_variables_and_constraints(const StaticHypergraph &hypergraph);

    // Add the edge variables and constraints to the ILP
    void add_edge_variables_and_constraints(const StaticHypergraph &hypergraph);

    // Set the objective function of the ILP
    void set_objective();

    // Optimize the ILP
    void optimize();

    // Get the result of the ILP
    CutValue get_result(const StaticHypergraph &hypergraph);
};

#endif // end of SMHM_ILP_H