/******************************************************************************
 * core_decomposition.h
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

#ifndef SMHM_CORE_DECOMPOSITION_H
#define SMHM_CORE_DECOMPOSITION_H

#include <vector>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/data_structure/union_find/union_find_sequential.h"
// Mt-KaHyPar headers
#include "mt-kahypar-library/libmtkahypar.h"

class CoreDecomposition
{
private:
    // Stores the degrees of the hypernodes
    std::vector<EdgeIndex> degrees;
    // Stores the nodes sorted by their unweighted degree
    std::vector<NodeID> sortedNodes;
    // Stores the position of the sorted nodes
    std::vector<NodeIndex> position;
    // Bucket that is used for sorting the nodes using the unweighted degree as key
    std::vector<NodeIndex> buckets;
    // Stores the sizes of the hyperedges
    std::vector<NodeIndex> edgeSizes;
    // Initial number of nodes in the hypergraph
    const NodeIndex initialNumNodes;
    // Initial number of edges in the hypergraph
    const NodeIndex initialNumEdges;
    // Initial number of pins in the hypergraph
    const NodeIndex initialNumPins;

public:
    CoreDecomposition(const NodeIndex initialNumNodes, const EdgeIndex initialNumEdges, const NodeIndex initialNumPins);

    // Decompose the hypergraph into its k-core components
    void decompose_into_cores(StaticHypergraph &hypergraph);

    // Get the target cores of the decomposition
    std::vector<EdgeIndex> get_target_cores();

    // Create the k-core hypergraph for a given k
    mt_kahypar_hypergraph_t create_k_core_graph(const EdgeIndex k, StaticHypergraph &hypergraph);
};

// Create the k-core hypergraph for a given k
inline mt_kahypar_hypergraph_t CoreDecomposition::create_k_core_graph(const EdgeIndex k, StaticHypergraph &hypergraph)
{
    // Compute the number of nodes in the k-core hypergraph
    // NB: buckets[k] contains the position of the first node that is part of the k-core
    // Every node with a position >= buckets[k] is part of the k-core
    NodeIndex numNodes = initialNumNodes - buckets[k];

    // Initialize the union-find data structure, which keeps track of the overlapping edges of the hypergraph
    // This is necessary to find the largest SCC in the k-core hypergraph
    UnionFindSequential unionFind(initialNumEdges);

    // For every node, store the first edge that contains the node
    // This is necessary to find the largest SCC in the k-core hypergraph
    std::vector<EdgeID> firstEdge(numNodes, std::numeric_limits<EdgeID>::max());

    for (EdgeID edgeID : hypergraph.edges())
    {
        // Ignore disabled edges
        if (!hypergraph.edgeIsEnabled(edgeID))
            continue;

        // Loop over the pins of the edge
        for (NodeID pinID : hypergraph.pins(edgeID))
            // Check if the pin is enabled and part of the k-core
            if (hypergraph.nodeIsEnabled(pinID) && degrees[pinID] >= k)
            {
                // NB: The new id of the pin is the position of the pin relative to buckets[k], i.e. its distance to buckets[k]
                NodeID newPinID = position[pinID] - buckets[k];
                // Check if the pin was first encountered by the current edge
                if (firstEdge[newPinID] == std::numeric_limits<EdgeID>::max())
                    // If so, store the edge
                    firstEdge[newPinID] = edgeID;
                else
                    // Otherwise, merge the edges
                    unionFind.Union(firstEdge[newPinID], edgeID);
            }
    }

    // Find the largest SCC in the k-core hypergraph
    // We only need to do this if there are more than one SCCs in the k-core hypergraph
    // NB: The ids of the SCCs are edge ids, since we use the union-find data structure to find overlapping edges
    EdgeID largestSCC = std::numeric_limits<EdgeID>::max();
    // Check if there is more than one SCC in the k-core graph
    if (unionFind.n() > 1)
    {
        std::vector<EdgeID> sizeOfSCC(initialNumEdges, 0);
        for (NodeIndex i = 0; i < numNodes; i++)
        {
            // Get the id of the SCC of the node (i.e. the parent edge of the union-find data structure)
            EdgeID nodeSCC = unionFind.Find(firstEdge[i]);
            // Increase the size of the SCC
            sizeOfSCC[nodeSCC]++;
            // Check if the current SCC is the largest one
            if (largestSCC == std::numeric_limits<EdgeID>::max() || sizeOfSCC[nodeSCC] > sizeOfSCC[largestSCC])
                largestSCC = nodeSCC;
        }
    }

    // Initialize the adjacency list of the k-core hypergraph and the list of the weight of the edges and of the nodes
    // NB: As the sizes, we use initialNumEdges and initialNumPins since we do not know the exact number of edges and pins in largest SCC of the k-core hypergraph
    //     The actual number of edges and pins will be computed while building the adjacency list.
    //     The space that was allocated too much will be simply ignored by MT-KaHyPar.
    std::unique_ptr<size_t[]> edgesAdjIndices = std::make_unique<size_t[]>(initialNumEdges + 1);
    std::unique_ptr<mt_kahypar_hyperedge_id_t[]> edgesAdjList = std::make_unique<mt_kahypar_hyperedge_id_t[]>(initialNumPins);
    std::unique_ptr<mt_kahypar_hyperedge_weight_t[]> edgesWeights = std::make_unique<mt_kahypar_hyperedge_weight_t[]>(initialNumEdges);
    std::unique_ptr<mt_kahypar_hypernode_weight_t[]> nodesWeights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(numNodes);

    // Count the number of edges and pins added to the adjacency list
    EdgeIndex numEdgesAddedToAdjList = 0;
    NodeIndex numPinsAddedToAdjList = 0;

    // Build the adjacency list (will be used to build the hypergraph)
    for (EdgeID edgeID : hypergraph.edges())
    {
        // Ignore disabled edges or edges that are not part of the largest SCC (if there is more than one SCC)
        if (!hypergraph.edgeIsEnabled(edgeID) || (unionFind.n() > 1 && unionFind.Find(edgeID) != largestSCC))
            continue;

        // Store the number of pins of the current edge that are part of the k-core
        NodeIndex numPinsOfEdgeAddedToKCore = 0;

        // Loop over the pins of the edge
        for (NodeID pinID : hypergraph.pins(edgeID))
            // Check if the pin is enabled and part of the k-core
            if (hypergraph.nodeIsEnabled(pinID) && degrees[pinID] >= k)
            {
                // If so, add the pin to the adjacency list and store its weight
                // NB: The new id of the pin is the position of the pin relative to buckets[k], i.e. its distance to buckets[k]
                NodeID newPinID = position[pinID] - buckets[k];
                edgesAdjList[numPinsAddedToAdjList++] = newPinID;
                numPinsOfEdgeAddedToKCore++;
                nodesWeights[newPinID] = hypergraph.nodeWeight(pinID);
            }

        // Check if the edge has at least two pins in the k-core
        if (numPinsOfEdgeAddedToKCore >= 2)
        {
            // If so, store the position of the first pin of the edge in the adjacency list
            edgesAdjIndices[numEdgesAddedToAdjList] = numPinsAddedToAdjList - numPinsOfEdgeAddedToKCore;
            // Also store the weight of the edge
            edgesWeights[numEdgesAddedToAdjList] = hypergraph.edgeWeight(edgeID);
            //  Increase the number of edges added to the adjacency list
            numEdgesAddedToAdjList++;
        }
        else
        {
            // If not, remove the pins from the adjacency list
            numPinsAddedToAdjList -= numPinsOfEdgeAddedToKCore;
        }
    }
    // Add the sentinel that points to the index behind the last pin
    edgesAdjIndices[numEdgesAddedToAdjList] = numPinsAddedToAdjList;

    // Create and return the k-core hypergraph
    return mt_kahypar_create_hypergraph(mt_kahypar_preset_type_t::DETERMINISTIC,
                                        numNodes,
                                        numEdgesAddedToAdjList,
                                        edgesAdjIndices.get(),
                                        edgesAdjList.get(),
                                        edgesWeights.get(),
                                        nodesWeights.get());
}

#endif // end of SMHM_CORE_DECOMPOSITION_H