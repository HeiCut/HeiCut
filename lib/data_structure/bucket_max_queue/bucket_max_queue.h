/******************************************************************************
 * bucket_max_queue.h
 * *
 * Bucket maximum priority queue for integer keys.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_BUCKET_MAX_QUEUE_H
#define SMHM_BUCKET_MAX_QUEUE_H

#include <vector>
// Own headers
#include "lib/utils/definitions.h"

template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
class BucketMaxQueue
{
private:
    // Store the (initial) number of elements that will be inserted into the bucket priority queue
    // NB: It is important to use the initial number of elements for the allocation here,
    //     because we use the element ID as the index of the vector, which is might be larger
    //     than the current number of elements
    ElementIndexType initialNumElements;
    // Store the buckets of the bucket priority queue
    std::vector<std::vector<ElementIDType>> buckets;
    // Store for each node a pair containing first the index of the bucket it is in
    // and second the index of the node inside of the bucket
    // NB: The first index is std::numeric_limits<BucketIndexType>::max() if the node is not in a bucket
    std::vector<std::pair<BucketIndexType, ElementIndexType>> nodeBucketLocation;
    // Store the highest index of a non-empty bucket
    BucketIndexType highestNonEmptyBucketIndex = 0;

public:
    BucketMaxQueue(ElementIndexType initialNumElements);

    // Reset/resize the bucket priority queue and reserve the space of the first bucket to hold all elements
    void reset(BucketIndexType numBuckets, ElementIndexType numElements);

    // Insert an element into the first bucket (after resetting the bucket queue)
    void insertInFirstBucket(ElementIDType elementID, ElementIndexType indexInFirstBucket);

    // Return and remove the top element of the bucket priority queue
    ElementIDType topAndPop(bool shouldDecreaseHighestNonEmptyBucketIndex = true);

    // Check if the bucket priority queue contains an element
    inline bool contains(ElementIDType elementID) const;

    // Increase the key of an element in the bucket priority queue
    inline void increaseByKey(ElementIDType elementID, BucketIndexType numLevelsToElevate);

    // Decrease the highest non-empty bucket index
    inline void decreaseHighestNonEmptyBucketIndex();

    // Mark element as not in a bucket anymore
    inline void markElementAsNotInBucket(ElementIDType elementID);
};

template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::BucketMaxQueue(ElementIndexType initialNumElements) : initialNumElements(initialNumElements), nodeBucketLocation(initialNumElements){};

// Reset/resize the bucket priority queue and reserve the space of the first bucket to hold all elements
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline void BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::reset(BucketIndexType numBuckets, ElementIndexType numElements)
{
    // Resize the bucket queue
    buckets.resize(numBuckets);
    // Insert all nodes into the first bucket
    buckets[0].resize(numElements);
    // Reset the highest non-empty bucket index to 0
    highestNonEmptyBucketIndex = 0;
}

// Insert an element into the first bucket (after resetting the bucket queue)
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline void BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::insertInFirstBucket(ElementIDType elementID, ElementIndexType indexInFirstBucket)
{
    // Insert the element into the first bucket
    buckets[0][indexInFirstBucket] = elementID;
    // Set the bucket location of the element
    nodeBucketLocation[elementID] = {0, indexInFirstBucket};
}

// Return the top element of the bucket priority queue
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline ElementIDType BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::topAndPop(bool shouldDecreaseHighestNonEmptyBucketIndex)
{
    ElementIDType elementID = buckets[highestNonEmptyBucketIndex].back();
    // Remove the node from the bucket queue
    buckets[highestNonEmptyBucketIndex].pop_back();
    // Store that the element is not in a bucket anymore
    markElementAsNotInBucket(elementID);
    // Decrease the highest non-empty bucket index (if necessary)
    if (shouldDecreaseHighestNonEmptyBucketIndex)
        decreaseHighestNonEmptyBucketIndex();
    // Return the element ID
    return elementID;
}

// Check if the bucket priority queue contains an element
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline bool BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::contains(ElementIDType elementID) const
{
    return nodeBucketLocation[elementID].first != std::numeric_limits<BucketIndexType>::max();
}

// Increase the key of an element in the bucket priority queue
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline void BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::increaseByKey(ElementIDType elementID, BucketIndexType numLevelsToElevate)
{
    // Get the index of the bucket and the index of the element inside of the bucket
    BucketIndexType bucketIndex = nodeBucketLocation[elementID].first;
    ElementIndexType elementIndex = nodeBucketLocation[elementID].second;
    // Compute the index of the next bucket to insert the element
    BucketIndexType nextBucketIndex = bucketIndex + numLevelsToElevate;

    // Check if the bucket contains more than one node
    if (buckets[bucketIndex].size() > 1)
    {
        // Update the bucket location of the last element in the bucket (will be swapped)
        nodeBucketLocation[buckets[bucketIndex].back()].second = nodeBucketLocation[elementID].second;
        // Swap the element with the last node in the bucket
        std::swap(buckets[bucketIndex][elementIndex], buckets[bucketIndex].back());
    }

    // Update the bucket location of the element
    nodeBucketLocation[elementID] = {nextBucketIndex, buckets[nextBucketIndex].size()};
    // Remove the element from the old bucket
    buckets[bucketIndex].pop_back();
    // Insert the element into the next bucket
    buckets[nextBucketIndex].push_back(elementID);
    // Increase the highest non-empty bucket index (if necessary)
    if (nextBucketIndex > highestNonEmptyBucketIndex)
        highestNonEmptyBucketIndex = nextBucketIndex;
}

// Decrease the highest non-empty bucket index
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline void BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::decreaseHighestNonEmptyBucketIndex()
{
    while (buckets[highestNonEmptyBucketIndex].empty() && highestNonEmptyBucketIndex > 0)
        highestNonEmptyBucketIndex--;
}

// Mark element as not in a bucket anymore
template <typename ElementIDType, typename ElementIndexType, typename BucketIndexType>
inline void BucketMaxQueue<ElementIDType, ElementIndexType, BucketIndexType>::markElementAsNotInBucket(ElementIDType elementID)
{
    nodeBucketLocation[elementID].first = std::numeric_limits<BucketIndexType>::max();
}

#endif // end of SMHM_BUCKET_MAX_QUEUE_H