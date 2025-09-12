/******************************************************************************
 * submodular.h
 * *
 * Solver that finds the minimum cut of a hypergraph via submodular optimization.
 * More specifically, it implements the algorithms of Klimmek and Wagner, Mak
 * and Wong and Queyranne. Note that one instance of the SubmodularMincut class
 * can only be used for one hypergraph (and its contracted versions).
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_SUBMODULAR_H
#define SMHM_SUBMODULAR_H

#include <cassert>
#include <vector>
// Own headers
#include "lib/utils/definitions.h"
#include "lib/orderer/orderer.h"
#ifdef SMHM_PARALLEL
#include <atomic>
#include <mutex>
#include <condition_variable>
#endif

struct SubmodularMincutResult
{
    CutValue minEdgeCut;
    IterationIndex numIterations;
    double meanContractionsPerIteration;
};

class SubmodularMincut
{
private:
    // Mode of the contraction performed with the ordering
    const OrderingMode orderingMode;
#ifdef SMHM_PARALLEL
    // Orderer for computing the ordering of the nodes for each thread
    // NB: We use a thread-indexed vector instead of enumerable_thread_specific so that we can
    //     reuse the data structures for multiple .solve() calls
    std::vector<std::unique_ptr<OrdererWrapper<DynamicHypergraph>>> ordererPerThread;
    // Store the ordering of the nodes for each thread
    // NB: We use a thread-indexed vector instead of enumerable_thread_specific so that we can
    //     reuse the data structures for multiple .solve() calls
    std::vector<std::vector<NodeID>> nodeOrderingPerThread;
    // Intitialize the hash vector to store the total number of marked pins of each hyperedge for the current last node in the ordering.
    // NB: The first element of the pair is the total number of marked pins of the hyperedge (whose id is equal to index),
    //     the second element is the current last node in the ordering. This allows us to use the same vector for all hyperedges,
    //     since we exactly know when to to reset the values in the vector.
    // NB2: We use a thread-indexed vector instead of enumerable_thread_specific so that we can
    //     reuse the data structures for multiple .solve() calls
    std::vector<std::vector<std::pair<NodeIndex, NodeID>>> markedEdgePinsPerThread;
    // Store the number of threads
    const size_t numThreads;
    // Store the threads (excluding the main thread)
    std::vector<std::thread> threads;
    // Atomically store the number of contractions
    std::atomic<NodeIndex> numContractions{0};
    // Variables for the first barrier
    std::mutex firstMutex;
    std::condition_variable firstConditionVariable;
    int numThreadsReachedFirstBarrier = 0;
    bool readyToLeaveFirstBarrier = false;
    //  Variables for the second barrier
    std::mutex secondMutex;
    std::condition_variable secondConditionalVariable;
    int numThreadsReachedSecondBarrier = 0;
    bool readyToLeaveSecondBarrier = false;
#else
    // Orderer for computing the ordering of the nodes
    Orderer<DynamicHypergraph, EdgeWeight> orderer;
    // Store the ordering of the nodes
    std::vector<NodeID> nodeOrdering;
    // Intitialize the hash vector to store the total number of marked pins of each hyperedge for the current last node in the ordering.
    // NB: The first element of the pair is the total number of marked pins of the hyperedge (whose id is equal to index),
    //     the second element is the current last node in the ordering. This allows us to use the same vector for all hyperedges,
    //     since we exactly know when to to reset the values in the vector.
    std::vector<std::pair<NodeIndex, NodeID>> markedEdgePins;
#endif

    // Compute the cut of the phase (i.e. the cut edges when isolating the last node in the ordering)
    CutValue compute_cut_of_phase(DynamicHypergraph &hypergraph, const NodeID nodeID, std::vector<std::pair<NodeIndex, NodeID>> &markedEdgePins);

    // Try to contract further nodes in the ordering after contracting the last two nodes
    CutValue try_contracting_further_nodes(DynamicHypergraph &hypergraph, const NodeIndex numNodes, const CutValue cutOfThePhase, std::vector<NodeID> &nodeOrdering, std::vector<std::pair<NodeIndex, NodeID>> &markedEdgePins);

#ifdef SMHM_PARALLEL
    // Barrier function to synchronize threads
    void barrier(std::mutex &mutex, std::condition_variable &conditionVariable, int &numThreadsReachedBarrier, bool &readyToLeaveBarrier, std::function<void()> sequentialFunction = nullptr);

    // Solve the minimum cut problem in parallel via submodular optimization using multiple threads
    void solve_in_parallel(const size_t threadIndex, DynamicHypergraph &hypergraph, CutValue &minEdgeCut, NodeIndex &numNodes, IterationIndex &numIterations);
#endif

public:
    SubmodularMincut(const NodeIndex initialNumNodes,
                     const EdgeIndex numEdges,
                     const OrderingType orderingType,
                     const OrderingMode orderingMode,
                     const bool hasWeightedEdges,
                     const size_t numThreads = 1);

    // Solve the minimum cut problem via submodular optimization
    SubmodularMincutResult solve(DynamicHypergraph &hypergraph);
};

// Compute the cut of the phase (i.e. the cut edges when isolating the last node in the ordering)
inline CutValue SubmodularMincut::compute_cut_of_phase(DynamicHypergraph &hypergraph, const NodeID nodeID, std::vector<std::pair<NodeIndex, NodeID>> &markedEdgePins)
{
    // Get the cut of the phase (i.e. the cut edges when isolating the last node in the ordering)
    CutValue cutOfThePhase = 0;
    for (EdgeID edgeID : hypergraph.incidentEdges(nodeID))
        if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
        {
            cutOfThePhase += hypergraph.edgeWeight(edgeID);
            markedEdgePins[edgeID] = std::make_pair(1, nodeID);
        }
    return cutOfThePhase;
}

#ifdef SMHM_PARALLEL
// Barrier function to synchronize threads
inline void SubmodularMincut::barrier(std::mutex &mutex, std::condition_variable &conditionVariable, int &numThreadsReachedBarrier, bool &readyToLeaveBarrier, std::function<void()> sequentialFunction)
{
    // Lock the mutex
    std::unique_lock<std::mutex> lock(mutex);
    // Increase the number of threads that reached the barrier
    numThreadsReachedBarrier++;
    // If all threads reached the barrier, notify all threads and reset the flag
    if (numThreadsReachedBarrier == numThreads)
    {
        // If a sequential function is provided, call it
        if (sequentialFunction)
            sequentialFunction();
        numThreadsReachedBarrier = 0;
        readyToLeaveBarrier = true;
        conditionVariable.notify_all();
    }
    else
    {
        // Wait until all threads reached the barrier
        conditionVariable.wait(lock, [&]
                               { return readyToLeaveBarrier; });
    }
}

// Solve the minimum cut problem in parallel via submodular optimization using multiple threads
inline void SubmodularMincut::solve_in_parallel(const size_t threadIndex, DynamicHypergraph &hypergraph, CutValue &minEdgeCut, NodeIndex &numNodes, IterationIndex &numIterations)
{
    while (numNodes > 1)
    {
        // Get the node ordering for the thread
        std::vector<NodeID> &nodeOrdering = nodeOrderingPerThread[threadIndex];
        // Get the orderer for the thread
        OrdererWrapper<DynamicHypergraph> &orderer = *ordererPerThread[threadIndex];
        // Get the marked edges for the thread
        std::vector<std::pair<NodeIndex, NodeID>> &markedEdgePins = markedEdgePinsPerThread[threadIndex];
        // Store whether the contraction could be registered
        bool couldRegisterContraction = false;
        // Compute the ordering of the nodes in hypergraph (depending on the selected ordering type)
        orderer.compute_ordering(hypergraph, numNodes, &nodeOrdering, nullptr, nullptr);
        // Compute the cut of the phase (i.e. the cut edges when isolating the last node in the ordering)
        CutValue cutOfThePhase = compute_cut_of_phase(hypergraph, nodeOrdering[numNodes - 1], markedEdgePins);
        // Contract the last two nodes in the ordering (if possible)
        if (hypergraph.registerContraction(nodeOrdering[numNodes - 2], nodeOrdering[numNodes - 1]))
        {
            couldRegisterContraction = true;
            if (orderingMode == OrderingMode::MULTI)
                // If the contraction could be registered, we check if the edge list of the last node in the ordering
                // is entirely included in the edge list of the node before it in the ordering, because then the node
                // ordering will not change upon contraction and we can contract also the third last node in the ordering
                // We do this over and over again until the edge lists are not entirely included anymore
                cutOfThePhase = try_contracting_further_nodes(hypergraph, numNodes, cutOfThePhase, nodeOrdering, markedEdgePins);

            // Update the minimum edge cut (if necessary)
            // NB: We need to use a local variable minEdgeCutOfThread here and repeatedly try to update the global minEdgeCut to avoid race conditions
            CutValue minEdgeCutOfThread;
            do
            {
                minEdgeCutOfThread = minEdgeCut;
                // If we cannot even decrease the local minimum edge cut, we can break the loop
                if (cutOfThePhase >= minEdgeCutOfThread)
                    break;
                // Try to update the global minimum edge cut (if it is equal to the local minimum edge cut and thus larger than the cut of the phase)
            } while (!__sync_bool_compare_and_swap(&minEdgeCut, minEdgeCutOfThread, cutOfThePhase));
        }

        // Wait for all threads to finish the computation of the ordering
        barrier(firstMutex, firstConditionVariable, numThreadsReachedFirstBarrier, readyToLeaveFirstBarrier, [&]()
                {
                    // Increase the number of iterations
                    numIterations++;
                    // Reset the ready to leave flag for the second barrier
                    readyToLeaveSecondBarrier = false; });

        // If the contraction could be registered, perform the contraction and increase the number of contractions
        if (couldRegisterContraction)
            numContractions.fetch_add(hypergraph.contract(nodeOrdering[numNodes - 1]), std::memory_order_relaxed);

        // Wait for all threads to finish the contraction
        barrier(secondMutex, secondConditionalVariable, numThreadsReachedSecondBarrier, readyToLeaveSecondBarrier, [&]()
                {
                    // Update the number of nodes
                    numNodes -= numContractions.load(std::memory_order_relaxed);
                    // We should at least have one contraction per iteration
                    assert(numContractions > 0);
                    // Reset the number of contractions for the next iteration
                    numContractions.store(0);
                    // Reset the ready to leave flag for the first barrier
                    readyToLeaveFirstBarrier = false; });
    }
}
#endif

// Try to contract further nodes in the ordering after contracting the last two nodes
inline CutValue SubmodularMincut::try_contracting_further_nodes(DynamicHypergraph &hypergraph, const NodeIndex numNodes, const CutValue cutOfThePhase, std::vector<NodeID> &nodeOrdering, std::vector<std::pair<NodeIndex, NodeID>> &markedEdgePins)
{
    // Store the current cut of the phase
    CutValue currentCutOfThePhase = cutOfThePhase;
    // Store the index of the node in the ordering for which we check if its edge list entirely convers the edge list of the (contracted) nodes coming AFTER it
    NodeIndex nextNodeIndex = numNodes - 2;
    // Get the degree of all of the (contracted) nodes coming AFTER the nextNodeIndex in the ordering
    // NB: For the initalization, this is simply the degree of the last node in the ordering BUT without edges of size 1 or less
    EdgeIndex currentDegree = 0;
    for (EdgeID edgeID : hypergraph.incidentEdges(nodeOrdering[numNodes - 1]))
        if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
            currentDegree++;

    // Try to check for further contractions
    while (nextNodeIndex > 0)
    {
        // NB: We need to use a copy of the cut of the phase here, because we only update it if the new cut of the phase is smaller
        CutValue nextCutOfThePhase = currentCutOfThePhase;
        // NB: We need to use a copy of the degree here, because we only update it if we performed the contraction
        EdgeIndex nextDegree = currentDegree;
        // Store the number of shared edges between the node at nextNodeIndex and the (contracted) nodes coming AFTER it in the ordering
        EdgeIndex numSharedEdges = 0;
        for (EdgeID edgeID : hypergraph.incidentEdges(nodeOrdering[nextNodeIndex]))
            if (hypergraph.edgeIsEnabled(edgeID) && hypergraph.edgeSize(edgeID) > 1)
                // Check if the the edge is marked
                if (markedEdgePins[edgeID].second == nodeOrdering[numNodes - 1])
                {
                    // If so, increase the number of shared edges and the number of marked pins of the edge
                    numSharedEdges++;
                    markedEdgePins[edgeID].first++;
                    // If all pins of the edge are marked, we need to decrease the cut of the phase since the edge will disappear upon contraction
                    if (hypergraph.edgeSize(edgeID) == markedEdgePins[edgeID].first)
                    {
                        nextCutOfThePhase -= hypergraph.edgeWeight(edgeID);
                        nextDegree--;
                        markedEdgePins[edgeID] = std::make_pair(0, std::numeric_limits<NodeID>::max());
                    }
                }
                else
                {
                    // If the edge is not marked, we need to increase the cut of the phase since the edge will be cut upon contraction
                    nextCutOfThePhase += hypergraph.edgeWeight(edgeID);
                    nextDegree++;
                    markedEdgePins[edgeID] = std::make_pair(1, nodeOrdering[numNodes - 1]);
                }
        // Stop if the edge list of the node at nextNodeIndex does not entirely cover the edge list of the (contracted) nodes coming AFTER it
        if (numSharedEdges != currentDegree)
            break;
        // Stop if the contraction could not be registered
        if (!hypergraph.registerContraction(nodeOrdering[nextNodeIndex - 1], nodeOrdering[nextNodeIndex]))
            break;
        // If the contraction could be registered, update the cut of the phase (if necessary) and the degree
        currentCutOfThePhase = std::min(currentCutOfThePhase, nextCutOfThePhase);
        currentDegree = nextDegree;
        --nextNodeIndex;
    }
    // Return the cut of the phase
    return currentCutOfThePhase;
}

// Solve the minimum cut problem via submodular optimization
inline SubmodularMincutResult SubmodularMincut::solve(DynamicHypergraph &hypergraph)
{
    // Extract the number of nodes
    // NB: Since the passed hypergraph could be a contracted version of the original hypergraph,
    //     we cannot use the initialNumNodes passed to the constructor, since this is is the initial
    //     number of nodes of the original hypergraph
    NodeIndex numNodes = hypergraph.initialNumNodes() - hypergraph.numRemovedHypernodes();

    // Store a copy of the number of nodes for computing the mean contractions per iteration
    NodeIndex numNodesBeforeSolver = numNodes;

    // Initialize the minimum edge cut to the maximum possible value
    CutValue minEdgeCut = std::numeric_limits<CutValue>::max();

    // Count the number of iterations
    IterationIndex numIterations = 0;

#ifdef SMHM_PARALLEL
    // Spawn numThreads - 1 threads to solve the minimum cut problem in parallel
    for (size_t threadIndex = 1; threadIndex < numThreads; ++threadIndex)
        threads.emplace_back([&, threadIndex]()
                             { solve_in_parallel(threadIndex, hypergraph, minEdgeCut, numNodes, numIterations); });

    // Also use the main thread to solve the minimum cut problem in parallel
    solve_in_parallel(0, hypergraph, minEdgeCut, numNodes, numIterations);

    // Wait for all threads to finish
    for (auto &t : threads)
        t.join();
#else
    // Repeatedly contract the last two nodes in the node ordering until only one node is left
    while (numNodes > 1)
    {
        numIterations++;
        // Compute the ordering of the nodes in hypergraph (depending on the selected ordering type)
        orderer.compute_ordering(hypergraph, numNodes, &nodeOrdering, nullptr, nullptr);
        // Compute the cut of the phase (i.e. the cut edges when isolating the last node in the ordering)
        CutValue cutOfThePhase = compute_cut_of_phase(hypergraph, nodeOrdering[numNodes - 1], markedEdgePins);
        //  Contract the last two nodes in the ordering (if possible)
        if (hypergraph.registerContraction(nodeOrdering[numNodes - 2], nodeOrdering[numNodes - 1]))
        {
            if (orderingMode == OrderingMode::MULTI)
                // If the contraction could be registered, we check if the edge list of the last node in the ordering
                // is entirely included in the edge list of the node before it in the ordering, because then the node
                // ordering will not change upon contraction and we can contract also the third last node in the ordering
                // We do this over and over again until the edge lists are not entirely included anymore
                cutOfThePhase = try_contracting_further_nodes(hypergraph, numNodes, cutOfThePhase, nodeOrdering, markedEdgePins);
            // If the contraction could be registered, perform the contraction and decrease the number of nodes
            NodeID numContractions = hypergraph.contract(nodeOrdering[numNodes - 1]);
            numNodes -= numContractions;
            // We should at least have one contraction per iteration
            assert(numContractions > 0);

            // Update the minimum edge cut (if necessary)
            if (cutOfThePhase < minEdgeCut)
                minEdgeCut = cutOfThePhase;
        }
    }
#endif
    return {minEdgeCut, numIterations, static_cast<double>(numNodesBeforeSolver - 1) / numIterations};
}

#endif // end of SMHM_SUBMODULAR_H