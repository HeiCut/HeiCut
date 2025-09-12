/******************************************************************************
 * ilp.cpp
 * *
 * ILP that solves the mincut problem on a hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "ilp.h"
#include "lib/utils/definitions.h"
// Gurobi headers
#include "gurobi_c++.h"

ILPMincut::ILPMincut(const NodeIndex numNodes, const EdgeIndex numEdges, const GRBEnv env, const int seed, const double ilpTimeout, const ILPMode ilpMode, const size_t numThreads)
    : numNodes(numNodes),
      numEdges(numEdges),
      nodes(numNodes),
      edges(numEdges),
      minimums(numEdges),
      maximums(numEdges),
      env(env),
      model(env),
      ilpTimeout(ilpTimeout),
      ilpMode(ilpMode)
{
    // Set the name of the model
    model.set(GRB_StringAttr_ModelName, "HypergraphMinimumCut");
    // Set the MIP gap of the model to 0 (i.e. accept only exact solution)
    model.set(GRB_DoubleParam_MIPGap, 0);
    // Set smaller tolerances, since the variables are not real bits but rather doubles and
    // rounding errors may lead to non-optimal solutions (e.g. all nodes on one side of the cut)
    model.set(GRB_DoubleParam_IntFeasTol, 1e-7);
    model.set(GRB_DoubleParam_FeasibilityTol, 1e-7);
    // Disable console logging
    model.set(GRB_IntParam_LogToConsole, 0);
    // Set the time limit of the model
    model.set(GRB_DoubleParam_TimeLimit, ilpTimeout);
    // Set the seed
    model.set(GRB_IntParam_Seed, seed);
    // Set the number of threads
    model.set(GRB_IntParam_Threads, numThreads);
};

// Add the node variables and constraints to the ILP
void ILPMincut::add_node_variables_and_constraints(const StaticHypergraph &hypergraph)
{
    for (NodeID nodeID : hypergraph.nodes())
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Add binary variable for the node
        nodes[nodeID] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        // Set the start value of the variable to 0
        nodes[nodeID].set(GRB_DoubleAttr_Start, 0.0);
        // Add the node to the total number of nodes in partition 1 if the node is in the partition 1
        nodeTot += nodes[nodeID];
    }

    // Add balance constraints
    model.addConstr(nodeTot >= 1, "lower_bound_balance");
    model.addConstr(nodeTot <= numNodes - 1, "upper_bound_balance");
};

// Add the edge variables and constraints to the ILP
void ILPMincut::add_edge_variables_and_constraints(const StaticHypergraph &hypergraph)
{
    for (EdgeID edgeID : hypergraph.edges())
    {
        // Ignore disabled edges
        if (!hypergraph.edgeIsEnabled(edgeID))
            continue;

        // Add binary variable for the edge
        edges[edgeID] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
        // Set the start value of the variable to 0
        edges[edgeID].set(GRB_DoubleAttr_Start, 0.0);
        // Add the weight of the edge to the total weight of the cut if the edge is cut
        edgeTot += hypergraph.edgeWeight(edgeID) * edges[edgeID];

        if (ilpMode == ILPMode::MILP)
        {
            // Add binary variable for the minimum pin variable of the edge
            minimums[edgeID] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            // Set the start value of the variable to 0
            minimums[edgeID].set(GRB_DoubleAttr_Start, 0.0);
            // Add binary variable for the maximum pin variable of the edge
            maximums[edgeID] = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS);
            // Set the start value of the variable to 0
            maximums[edgeID].set(GRB_DoubleAttr_Start, 0.0);

            // Add minimum and maximum constraints for the edge
            for (NodeID pinID : hypergraph.pins(edgeID))
                if (hypergraph.nodeIsEnabled(pinID))
                {
                    model.addConstr(minimums[edgeID] <= nodes[pinID]);
                    model.addConstr(maximums[edgeID] >= nodes[pinID]);
                }
            // Add the cut constraint for the edge using the minimum and maximum variables
            model.addConstr(edges[edgeID] >= maximums[edgeID] - minimums[edgeID]);
        }
        else if (ilpMode == ILPMode::BIP)
            // Add cut constraints for each pair of nodes in the edge (i.e. linking node and edge variables)
            for (NodeID firstPinID : hypergraph.pins(edgeID))
                for (NodeID secondPinID : hypergraph.pins(edgeID))
                    if (firstPinID != secondPinID && hypergraph.nodeIsEnabled(firstPinID) && hypergraph.nodeIsEnabled(secondPinID))
                        model.addConstr(edges[edgeID] >= nodes[firstPinID] - nodes[secondPinID]);
    }
};

// Set the objective function of the ILP
void ILPMincut::set_objective()
{
    // Set the objective of the model to minimize the total weight of the cut
    model.setObjective(edgeTot, GRB_MINIMIZE);
};

// Optimize the ILP
void ILPMincut::optimize()
{
    // Optimize the model
    model.optimize();
};

// Print the results of the ILP
CutValue ILPMincut::get_result(const StaticHypergraph &hypergraph)
{
    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL && model.get(GRB_IntAttr_Status) != GRB_TIME_LIMIT)
    {
        std::cout << "ILP: No solution found.\n";
        return std::numeric_limits<CutValue>::max();
    }

    if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
        std::cout << "ILP: Time limit reached.\n";
    else
        std::cout << "ILP: Optimal solution found.\n";

    // for (NodeIndex i = 0; i < numNodes; i++)
    //     std::cout << "Node " << i + 1 << " is in subset S: " << nodes[i].get(GRB_DoubleAttr_X) << std::endl;

    CutValue minEdgeCut = 0;

    for (EdgeIndex j = 0; j < numEdges; j++)
    {
        if (edges[j].get(GRB_DoubleAttr_X) > 0.5)
            minEdgeCut += hypergraph.edgeWeight(j);
        // std::cout << "Edge " << j + 1 << " is cut: " << edges[j].get(GRB_DoubleAttr_X) << std::endl;
    }

    return minEdgeCut;
};
