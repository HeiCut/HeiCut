/******************************************************************************
 * trimmer.h
 * *
 * Computes a general data structure for a given hypergraph from which we can
 * extract the k-trimmer of the hypergraph for any k.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_TRIMMER_H
#define SMHM_TRIMMER_H

#include <cassert>
#include <vector>
// Own headers
#include "lib/utils/definitions.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"

class Trimmer
{
private:
    // Number of nodes in the hypergraph
    const NodeIndex numNodes;
    // Number of edges in the hypergraph
    const EdgeIndex numEdges;

    // Hypergraph which we want to trim
    const StaticHypergraph &hypergraph;
    // MA-ordering of the nodes
    std::vector<NodeID> &nodeMAOrdering;
    // Head ordering of the edges
    std::vector<EdgeID> &edgeHeadOrdering;
    // Head of the edges
    std::vector<NodeID> &edgeHead;
    // Preset of the generated hypergraph (defined whether the generated hypergraph is static or dynamic)
    const mt_kahypar_preset_type_t preset;

    // Stores the backward edges for each node
    std::vector<std::vector<NodeID>> nodeBackwardEdges;

public:
    Trimmer(const NodeIndex initialNumNodes,
               const NodeIndex numNodes,
               const EdgeIndex numEdges,
               const StaticHypergraph &hypergraph,
               std::vector<NodeID> &nodeMAOrdering,
               std::vector<EdgeID> &edgeHeadOrdering,
               std::vector<NodeID> &edgeHead,
               const mt_kahypar_preset_type_t preset);

    // Build the backward edges of the nodes
    void build_backward_edges();

    // Print the backward edges of the nodes
    void print_backward_edges();

    // Create k-trimmed certificate of the hypergraph
    mt_kahypar_hypergraph_t create_k_trimmed_certificate(const TrimmerValue k);
};

// Create k-trimmed certificate of the hypergraph
inline mt_kahypar_hypergraph_t Trimmer::create_k_trimmed_certificate(const TrimmerValue k)
{
    // Initialize the trimmed edges
    std::vector<std::vector<NodeID>> trimmedEdges(numEdges);
    // Initialize the metrics of the trimmed hypergraph
    NodeIndex numTrimmedNodes = numNodes;
    EdgeIndex numTrimmedEdges = 0;
    NodeIndex numTotalTrimmedEdgeSize = 0;

    // Create the trimmed hypergraph from the backward edges of the nodes
    for (NodeID nodeID : hypergraph.nodes())
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Get the number of remaining backward edges for the node
        TrimmerValue numRemainingBackwardEdges = std::min(k, (TrimmerValue)nodeBackwardEdges[nodeID].size());
        for (EdgeIndex j = 0; j < numRemainingBackwardEdges; j++)
        {
            EdgeID edgeID = nodeBackwardEdges[nodeID][j];
            // Add the header of the trimmed edge if it is empty
            // Performing this operation here ensures that the header is only added
            // to those trimmed edges that have at least one pin besides the header,
            // meaning that we ignore trimmed edges of size 1
            if (trimmedEdges[edgeID].empty())
            {
                trimmedEdges[edgeID].push_back(edgeHead[edgeID]);
                numTrimmedEdges++;
                numTotalTrimmedEdgeSize++;
            }
            trimmedEdges[edgeID].push_back(nodeID);
            numTotalTrimmedEdgeSize++;
        }
    }

    // std::cout << "Trimmed edges:" << std::endl;
    // for (EdgeIndex j = 0; j < numEdges; j++)
    // {
    //     std::cout << "Edge " << j + 1 << ": ";
    //     for (NodeID pinID : trimmedEdges[j])
    //         std::cout << pinID + 1 << " ";
    //     std::cout << std::endl;
    // }

    std::unique_ptr<size_t[]> spasifiedEdgesAdjIndices = std::make_unique<size_t[]>(numTrimmedEdges + 1);
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> spasifiedEdgesAdjList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(numTotalTrimmedEdgeSize);
    // Count the number of edges and pins added to the adjacency list
    EdgeIndex numEdgesAddedToAdjList = 0;
    NodeIndex numPinsAddedToAdjList = 0;

    // Build the adjacency list (will be used to build the hypergraph)
    for (EdgeID edgeID : hypergraph.edges())
    {
        // Ignore disabled edges and trimmed edges of size 1 or less
        if (!hypergraph.edgeIsEnabled(edgeID) || trimmedEdges[edgeID].size() <= 1)
            continue;

        // Store the position of the first pin of the trimmed edge in the adjacency list
        spasifiedEdgesAdjIndices[numEdgesAddedToAdjList++] = numPinsAddedToAdjList;
        // Add the pins of the trimmed edge to the adjacency list
        for (NodeID pinID : trimmedEdges[edgeID])
            spasifiedEdgesAdjList[numPinsAddedToAdjList++] = pinID;
    }
    // Add the sentinel that points to the index behind the last pin
    spasifiedEdgesAdjIndices[numEdgesAddedToAdjList] = numPinsAddedToAdjList;

    // Make sure that the number of edges and pins added to the adjacency list is correct
    assert(numEdgesAddedToAdjList == numTrimmedEdges);
    assert(numPinsAddedToAdjList == numTotalTrimmedEdgeSize);

    // Create and return the k-trimmed hypergraph
    return mt_kahypar_create_hypergraph(preset,
                                        numTrimmedNodes,
                                        numTrimmedEdges,
                                        spasifiedEdgesAdjIndices.get(),
                                        spasifiedEdgesAdjList.get(),
                                        nullptr,
                                        nullptr);
}

#endif // end of SMHM_TRIMMER_H