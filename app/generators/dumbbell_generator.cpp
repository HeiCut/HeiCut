/******************************************************************************
 * dumbbell_generator.cpp
 * *
 * Generates a k-uniform dumbbell hypergraph with a given number of nodes per
 * side The dumbbell consists of two k-uniform hypergraphs (either complete or
 * dense depending on the configuration) connected by a single bridging hyperedge.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include <fstream>
#include <iostream>
#include <stack>
// Own headers
#include "lib/parse_parameters/parse_parameters.h"
#include "lib/utils/random.h"

// Generates a dense (but not complete) k-uniform hypergraph from a given range of nodes
// Dense means that each node has a minimum degree of minDegree
EdgeIndex generate_dense_k_uniform_subhypergraph(int start, int end, NodeIndex k, EdgeIndex minDegree, bool weighted, MersenneTwister &altRandEngine, std::ostringstream &buffer)
{
    std::vector<std::pair<NodeID, EdgeIndex>> nodesWithDegrees;
    // Build the nodes of the k-uniform subhypergraph
    for (NodeID i = start; i < end; ++i)
        nodesWithDegrees.push_back({i, 0});

    // Stores the number of edges
    EdgeIndex numEdges = 0;

    // Stores the swaps for undoing them
    std::stack<std::pair<NodeIndex, NodeIndex>> swapStack;

    for (NodeIndex i = 0; i < nodesWithDegrees.size(); ++i)
    {
        // Swap the current node to the first position
        std::swap(nodesWithDegrees[0], nodesWithDegrees[i]);
        swapStack.push(std::make_pair(0, i));
        // Add hyperedges until the minimum degree is reached
        while (nodesWithDegrees[0].second < minDegree)
        {
            // Select k-1 random nodes from the remaining nodes
            for (NodeIndex j = 1; j < k; ++j)
            {
                NodeIndex upperLimit = nodesWithDegrees.size() - 1;
                // Check if we are selecting the first random node for the first outgoing edge of the current node in the subhypergraph.
                // If so, we know the current node is not yet connected to the start node, which we fix by selecting one of the already connected nodes (i.e. upperLimit = i).
                // NB: We need to handle this to ensure that the subhypergraph is connected.
                if (i > 0 && nodesWithDegrees[0].second == 0 && j == 1)
                    upperLimit = i;
                NodeIndex randomNodeIndex = RandomFunctions::get_uniform_random_int_in_bounds<NodeIndex>(j, upperLimit);
                std::swap(nodesWithDegrees[j], nodesWithDegrees[randomNodeIndex]);
                swapStack.push(std::make_pair(j, randomNodeIndex));
            }
            // Write the hyperedge to the buffer
            if (weighted)
                buffer << RandomFunctions::get_uniform_random_int_in_bounds<EdgeWeight>(1, 100, altRandEngine) << " ";
            for (NodeIndex j = 0; j < k; ++j)
            {
                buffer << nodesWithDegrees[j].first + 1 << " ";
                nodesWithDegrees[j].second++;
            }
            buffer << std::endl;
            numEdges++;

            // Undo all but the last swap
            while (swapStack.size() > 1)
            {
                std::pair<NodeIndex, NodeIndex> swap = swapStack.top();
                std::swap(nodesWithDegrees[swap.first], nodesWithDegrees[swap.second]);
                swapStack.pop();
            }
        }
        // Undo the last swap
        std::pair<NodeIndex, NodeIndex> swap = swapStack.top();
        std::swap(nodesWithDegrees[swap.first], nodesWithDegrees[swap.second]);
        swapStack.pop();
    }

    return numEdges;
}

// Generates the complete k-uniform hypergraph from a given range of nodes
EdgeIndex generate_complete_k_uniform_subhypergraph(int start, int end, NodeIndex k, bool weighted, MersenneTwister &altRandEngine, std::ostringstream &buffer)
{
    std::vector<NodeID> nodes;
    // Build the nodes of the k-uniform subhypergraph
    for (NodeID i = start; i < end; ++i)
        nodes.push_back(i);

    // Stores the number of edges
    EdgeIndex numEdges = 0;

    // Stores the current combination/hyperedge of size k
    std::vector<NodeID> combination(k);

    // Dynamic for loop function that generates all possible combinations of size k
    std::function<void(NodeIndex, NodeIndex)> combine = [&](NodeIndex offset, NodeIndex depth)
    {
        // If the combination is complete, write it to the buffer
        if (depth == k)
        {
            if (weighted)
                buffer << RandomFunctions::get_uniform_random_int_in_bounds<EdgeWeight>(1, 100, altRandEngine) << " ";
            for (NodeIndex i = 0; i < k; ++i)
                buffer << combination[i] + 1 << " ";
            buffer << std::endl;
            numEdges++;
            return;
        }
        // Otherwise, add all possible nodes for the current depth
        for (NodeIndex i = offset; i <= nodes.size() - (k - depth); ++i)
        {
            combination[depth] = nodes[i];
            combine(i + 1, depth + 1);
        }
    };

    // Start the combination generation
    combine(0, 0);

    return numEdges;
}

int main(int argc, char *argv[])
{
#ifdef NDEBUG
    std::cout << "build_mode \t\t\t\tRELEASE" << std::endl;
#else
    std::cout << "build_mode \t\t\t\tDEBUG" << std::endl;
#endif

    // Stores the config
    DumbbellGeneratorConfig config;

    int parseReturnCode = ParseParams::parse_parameters_dumbbell_generator(argc, argv, config);

    if (parseReturnCode)
        return 1;

    // Set the seed for the random number generator
    RandomFunctions::set_seed(config.seed);
    // Use an alternative random engine for the weights, so that the hypergraph structure is identical for the unweighted and weighted case when using the same seed
    MersenneTwister altRandEngine;
    altRandEngine.seed(config.seed);

    if (config.numNodesPerSide < config.edgeSize)
    {
        std::cerr << "The number of nodes per side must be at least equal to the edge size of " << config.edgeSize << "." << std::endl;
        return 1;
    }

    if (config.edgeSize < 2)
    {
        std::cerr << "The edge size must be at least 2." << std::endl;
        return 1;
    }

    if (!config.outputFileName)
    {
        std::cerr << "No output file specified." << std::endl;
        return 1;
    }

    // Create a string stream buffer
    std::ostringstream buffer;
    // Open output file
    std::ofstream f(config.outputFileName);

    if (!f)
    {
        std::cerr << "Failed to open output file: " << config.outputFileName << std::endl;
        return 1;
    }

    // Store the number of nodes and edges
    NodeIndex numNodes = 2 * config.numNodesPerSide;
    EdgeIndex numEdges = 0;

    std::cout << "Generating " << config.edgeSize << "-uniform dumbbell hypergraph with " << config.numNodesPerSide << " nodes per side." << std::endl;

    if (config.minDegree > 0)
        std::cout << "Minimum degree per node: " << config.minDegree << std::endl;
    else
        std::cout << "Each side of the dumbbell is a complete hypergraph." << std::endl;

    // Generate the left complete k-uniform hypergraph
    if (config.minDegree > 0)
        numEdges += generate_dense_k_uniform_subhypergraph(0, config.numNodesPerSide, config.edgeSize, config.minDegree, config.weighted, altRandEngine, buffer);
    else
        numEdges += generate_complete_k_uniform_subhypergraph(0, config.numNodesPerSide, config.edgeSize, config.weighted, altRandEngine, buffer);

    // Generate the right complete k-uniform hypergraph
    if (config.minDegree > 0)
        numEdges += generate_dense_k_uniform_subhypergraph(config.numNodesPerSide, 2 * config.numNodesPerSide, config.edgeSize, config.minDegree, config.weighted, altRandEngine, buffer);
    else
        numEdges += generate_complete_k_uniform_subhypergraph(config.numNodesPerSide, 2 * config.numNodesPerSide, config.edgeSize, config.weighted, altRandEngine, buffer);

    // Generate a random bridge hyperedge
    // Initialize the nodes contain 0, 1, 2, ..., numNodes-1
    std::vector<NodeID> nodes(numNodes);
    std::iota(nodes.begin(), nodes.end(), 0);
    // Add a node from the left side the bridge hyperedge
    NodeIndex lefRandomNodeIndex = RandomFunctions::get_uniform_random_int_in_bounds<NodeIndex>(0, config.numNodesPerSide - 1);
    std::swap(nodes[0], nodes[lefRandomNodeIndex]);
    // Add a node from the right side the bridge hyperedge
    NodeIndex rightRandomNodeIndex = RandomFunctions::get_uniform_random_int_in_bounds<NodeIndex>(config.numNodesPerSide, 2 * config.numNodesPerSide - 1);
    std::swap(nodes[1], nodes[rightRandomNodeIndex]);
    // Add the remaining nodes to the bridge hyperedge
    for (NodeIndex i = 2; i < config.edgeSize; ++i)
    {
        NodeIndex randomNodeIndex = RandomFunctions::get_uniform_random_int_in_bounds<NodeIndex>(i, numNodes - 1);
        std::swap(nodes[i], nodes[randomNodeIndex]);
    }
    // Add the bridge hyperedge to the number of edges
    numEdges++;
    // Write the bridge hyperedge to the buffer
    if (config.weighted)
        buffer << 1 << " ";
    for (NodeIndex i = 0; i < config.edgeSize; ++i)
        buffer << nodes[i] + 1 << " ";

    buffer << std::endl;

    // Write the node weights to the buffer if necessary
    if (config.weighted)
        for (NodeIndex i = 0; i < numNodes; ++i)
            buffer << RandomFunctions::get_uniform_random_int_in_bounds<NodeWeight>(1, 100, altRandEngine) << std::endl;

    // Write the header to the file
    f << numEdges << " " << numNodes << (config.weighted ? " 11" : "") << std::endl;
    // Write the hyperedges to the file
    f << buffer.str();
    // Close the file
    f.close();

    std::cout << "Hypergraph has " << numNodes << " nodes and " << numEdges << " edges." << std::endl;
    if (config.weighted)
        std::cout << "Hypergraph is weighted." << std::endl;
    std::cout << "Hypergraph written to " << config.outputFileName << std::endl;

    return 0;
}