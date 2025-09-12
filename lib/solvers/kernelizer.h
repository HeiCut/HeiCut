/******************************************************************************
 * kernelizer.h
 * *
 * Computes the (estimated) mincut of a hypergraph. For this, the kernel of the
 * hypergraph is built by coarsening the hypergraph using the label propagation
 * and pruning rules. Afterwards, the mincut of the coarsened hypergraph is
 * computed using a submodular solver or an ILP (with floating point relaxation).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_KERNELIZER_H
#define SMHM_KERNELIZER_H

#include <iostream>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/coarsening/pruner.h"
// KaHIP headers
#include "kahip/timer.h"

struct KernelizerResult
{
    CutValue naiveEstimate;
    CutValue minEdgeCut;
    double time;
};

class Kernelizer
{
private:
    // Initialize the pruner
    Pruner pruner;
    // Store the timer
    timer t;
    // Store the config
    KernelizerConfig config;
    // Store the current estimate of the minimum edge cut
    CutValue minEdgeCut;

    // Update the minimum edge cut value (if it is smaller than the current estimate)
    // NB: We only do this if the coarsening has NOT reduced the entire hypergraph to a single hypernode
    void update_min_cut_value(Pruner &pruner, const StaticHypergraph &hypergraph);
    // Get the current number of edges of the hypergraph
    EdgeIndex get_current_num_edges(const StaticHypergraph &hypergraph);

    // Get the current number of nodes of the hypergraph
    NodeIndex get_current_num_nodes(const StaticHypergraph &hypergraph);

    // Whether we can stop early
    bool can_stop_early(const StaticHypergraph &hypergraph);

    // Print the statistics of the current hypergraph
    double print_stats_and_return_used_time(const StaticHypergraph &hypergraph,
                                            const NodeIndex initialNumNodes,
                                            const EdgeIndex initialNumEdges,
                                            timer &t,
                                            std::string suffix,
                                            std::string timeDescription);

public:
    Kernelizer(KernelizerConfig config);

    // Compute the minimum cut of the hypergraph using kernelization
    KernelizerResult compute_mincut(StaticHypergraph &hypergraph, const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges);
};

// Update the minimum edge cut value (if it is smaller than the current estimate)
// NB: We only do this if the coarsening has NOT reduced the entire hypergraph to a single hypernode
inline void Kernelizer::update_min_cut_value(Pruner &pruner, const StaticHypergraph &hypergraph)
{
    if (get_current_num_nodes(hypergraph) == 1)
        return;

    CutValue newNaiveMinEdgeCut = pruner.compute_naive_mincut_estimate(hypergraph);
    minEdgeCut = std::min(minEdgeCut, newNaiveMinEdgeCut);
}

// Get the current number of edges of the hypergraph
inline EdgeIndex Kernelizer::get_current_num_edges(const StaticHypergraph &hypergraph)
{
    return hypergraph.initialNumEdges() - hypergraph.numRemovedHyperedges();
}

// Get the current number of nodes of the hypergraph
inline NodeIndex Kernelizer::get_current_num_nodes(const StaticHypergraph &hypergraph)
{
    return hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();
}

// Whether we can stop early
inline bool Kernelizer::can_stop_early(const StaticHypergraph &hypergraph)
{
    // Stop early if the current estimate is already 0
    if (minEdgeCut == 0)
    {
        if (config.verbose)
            std::cout << "Mincut is already 0, no base solver necessary." << std::endl;
        return true;
    }

    // Stop early if we coarsened the hypergraph to a single node.
    // In this case, the minimum edge cut is exactly the current estimate and we do not have to use the ILP.
    if (get_current_num_nodes(hypergraph) == 1)
    {
        if (config.verbose)
            std::cout << "Hypergraph contracted to single node, no base solver necessary." << std::endl;
        return true;
    }

    // Stop early if we coarsened the hypergraph to nodes with no edges.
    // In this case, the minimum edge cut is exactly 0 and we do not have to use the ILP.
    if (get_current_num_edges(hypergraph) == 0)
    {
        if (config.verbose)
            std::cout << "Hypergraph contracted to zero edges, no base solver necessary." << std::endl;
        minEdgeCut = 0;
        return true;
    }

    return false;
}

// Print the statistics of the current hypergraph
inline double Kernelizer::print_stats_and_return_used_time(const StaticHypergraph &hypergraph,
                                                           const NodeIndex initialNumNodes,
                                                           const EdgeIndex initialNumEdges,
                                                           timer &t,
                                                           std::string suffix,
                                                           std::string timeDescription)
{
    double elapsedTime = t.elapsed();

    if (config.verbose)
    {
        EdgeIndex currentNumEdges = get_current_num_edges(hypergraph);
        NodeIndex currentNumNodes = get_current_num_nodes(hypergraph);

        std::cout << "abs_num_edges" << suffix << currentNumEdges << std::endl;
        std::cout << "rel_num_edges" << suffix << (100.0 * currentNumEdges) / initialNumEdges << std::endl;
        std::cout << "abs_num_nodes" << suffix << currentNumNodes << std::endl;
        std::cout << "rel_num_nodes" << suffix << (100.0 * currentNumNodes) / initialNumNodes << std::endl;

        std::cout << timeDescription << elapsedTime << std::endl;
        std::cout << "updated_naive_mincut_value \t\t" << minEdgeCut << std::endl;
        std::cout << "===================================================" << std::endl;
    }
    return elapsedTime;
}

#endif // end of SMHM_KERNELIZER_H