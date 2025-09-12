/******************************************************************************
 * trimmer.cpp
 * *
 * Computes a general data structure for a given hypergraph from which we can
 * extract the k-trimmer of the hypergraph for any k.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "trimmer.h"
#include "lib/utils/definitions.h"

Trimmer::Trimmer(const NodeIndex initialNumNodes,
                       const NodeIndex numNodes,
                       const EdgeIndex numEdges,
                       const StaticHypergraph &hypergraph,
                       std::vector<NodeID> &nodeMAOrdering,
                       std::vector<EdgeID> &edgeHeadOrdering,
                       std::vector<NodeID> &edgeHead,
                       const mt_kahypar_preset_type_t preset)
    : numNodes(numNodes),
      numEdges(numEdges),
      // NB: It is important to use the initial number of nodes for the allocation here,
      //     because we use the node ID as the index of the vector, which is might be larger
      //     than the number of nodes (e.g. if some nodes were contracted)
      nodeBackwardEdges(initialNumNodes),
      hypergraph(hypergraph),
      preset(preset),
      nodeMAOrdering(nodeMAOrdering),
      edgeHeadOrdering(edgeHeadOrdering),
      edgeHead(edgeHead) {};

// Build the backward edges of the nodes
void Trimmer::build_backward_edges()
{
    // Go over all edges and add the edge to the backward edges of the corresponding pins
    for (EdgeIndex i = 0; i < numEdges; i++)
    {
        EdgeID edgeID = edgeHeadOrdering[i];
        for (NodeID pinID : hypergraph.pins(edgeID))
            if (hypergraph.nodeIsEnabled(pinID) && pinID != edgeHead[edgeID])
                nodeBackwardEdges[pinID].push_back(edgeID);
    }
}

// Print the backward edges of the nodes
void Trimmer::print_backward_edges()
{
    std::cout << "Backward edges of the nodes:" << std::endl;
    for (NodeIndex i = 0; i < numNodes; i++)
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(i))
            continue;

        std::cout << "Node " << i + 1 << ": ";
        for (EdgeID edgeID : nodeBackwardEdges[i])
            std::cout << edgeID + 1 << " ";
        std::cout << std::endl;
    }
}