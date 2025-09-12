/******************************************************************************
 * core_decomposition.cpp
 * *
 * Decomposes a hypergraph into all its (k, 2)-core components, i.e. the maximal
 * subgraphs in which all nodes have at least an unweighted degree of k and all
 * hyperedges have at least a size of 2. This is the direct equivalent of k-core
 * definition on normal graphs, which was implemented in VieCut. This algorithm
 * is a port of the VieCut implementation for normal graphs to hypergraphs. The
 * VieCut implementation is based on the algorithm by Batagelj and Zaversnik.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "core_decomposition.h"
#include "lib/utils/definitions.h"

CoreDecomposition::CoreDecomposition(const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges, const NodeIndex initialNumPins)
    : initialNumNodes(initialNumNodes),
      initialNumEdges(initialNumEdges),
      initialNumPins(initialNumPins),
      degrees(initialNumNodes, 0),
      sortedNodes(initialNumNodes),
      position(initialNumNodes),
      buckets(1),
      edgeSizes(initialNumEdges, 0) {};

// Decompose the hypergraph into its k-core components
void CoreDecomposition::decompose_into_cores(StaticHypergraph &hypergraph)
{
    // Store the sizes of the hyperedges
    for (EdgeID edgeID : hypergraph.edges())
    {
        // Ignore disabled edges
        if (!hypergraph.edgeIsEnabled(edgeID))
            continue;

        edgeSizes[edgeID] = hypergraph.edgeSize(edgeID);
    }

    // Compute the maximum unweighted node degree and fill the buckets accordingly
    EdgeIndex maxUnweightedDegree = 0;
    for (NodeID nodeID : hypergraph.nodes())
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Compute the unweighted node degree of the node
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
            if (hypergraph.edgeIsEnabled(edgeID))
                degrees[nodeID] += 1;
        // Update the minimum cut estimate (if necessary)
        if (degrees[nodeID] > maxUnweightedDegree)
        {
            maxUnweightedDegree = degrees[nodeID];
            buckets.resize(maxUnweightedDegree + 1);
        }
        // Increase the bucket of the node
        buckets[degrees[nodeID]]++;
    }

    // Compute the prefix sum of the buckets
    NodeIndex start = 0;
    for (EdgeIndex i = 0; i < maxUnweightedDegree + 1; ++i)
    {
        NodeIndex num = buckets[i];
        buckets[i] = start;
        start += num;
    }

    // Perform bucket sort on the nodes
    for (NodeID nodeID : hypergraph.nodes())
    {
        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        position[nodeID] = buckets[degrees[nodeID]];
        sortedNodes[position[nodeID]] = nodeID;
        buckets[degrees[nodeID]]++;
    }

    // Reset the buckets
    for (EdgeIndex i = maxUnweightedDegree; i > 0; --i)
        buckets[i] = buckets[i - 1];
    buckets[0] = 0;

    // Iterate over the nodes in the sorted order
    for (NodeIndex i = 0; i < initialNumNodes; ++i)
    {
        NodeID nodeID = sortedNodes[i];

        // Ignore disabled nodes
        if (!hypergraph.nodeIsEnabled(nodeID))
            continue;

        // Iterate over the incident hyperedges of the node
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        {
            // Ignore disabled hyperedges or hyperedges of size less than 2
            if (!hypergraph.edgeIsEnabled(edgeID) || edgeSizes[edgeID] < 2)
                continue;

            // Decrease the size of the hyperedge
            // The idea is that the node cannot be part of the k-core with k > degrees[nodeID] so we remove it from the hyperedge
            edgeSizes[edgeID]--;

            // Check if the hyperedge is now of size 1
            if (edgeSizes[edgeID] == 1)
                // If so, the remaining pin has either a the same or a higher degree than the current node
                // If it has the same degree, there is nothing to do, because the remaining pin is already in the correct bucket
                // If it has a higher degree, we need to move it one bucket down, as the remaining pin has lost an edge (i.e. degree decreased by one)
                for (NodeID pinID : hypergraph.pins(edgeID))
                    if (hypergraph.nodeIsEnabled(pinID) && degrees[pinID] > degrees[nodeID])
                    {
                        // Get the degree of the remaining pin
                        EdgeIndex degree = degrees[pinID];
                        // Get the position of the remaining pin
                        NodeIndex pos = position[pinID];
                        // Get the position of the first node in the bucket of the remaining pin
                        NodeIndex start = buckets[degree];
                        // Get the vertex at the start position
                        NodeID startNodeID = sortedNodes[start];

                        // Check if the remaining pin is not at the start position of its bucket
                        if (pinID != startNodeID)
                        {
                            // Swap the remaining pin with the first node in the bucket
                            sortedNodes[start] = pinID;
                            sortedNodes[pos] = startNodeID;
                            position[pinID] = start;
                            position[startNodeID] = pos;
                        }
                        // Increase the start position of the old bucket (i.e. move first node to previous bucket)
                        buckets[degree]++;
                        // Decrease the degree of the remaining pin
                        degrees[pinID]--;
                    }
        }
    }
};

// Get the target cores of the decomposition
std::vector<EdgeIndex> CoreDecomposition::get_target_cores()
{
    // Initialize the target cores
    std::vector<EdgeIndex> targetCores;
    // Set the minimum core to 2 (i.e. all nodes have at least a degree of 2)
    EdgeIndex minCore = 2;
    // The maximum core is the degree of the last node in the sorted order, i.e. the maximum node degree
    EdgeIndex maxCore = degrees[sortedNodes[initialNumNodes - 1]];

    // Put every core in thg range [minCore, maxCore] into targetCores except of those cores where the corresponding
    // bucket is empty (i.e. its start position is equal to the start position of the previous bucket)
    targetCores.push_back(minCore);
    for (EdgeIndex i = minCore + 1; i <= maxCore; ++i)
        if (buckets[i] != buckets[i - 1])
            targetCores.push_back(i);

    return targetCores;
};