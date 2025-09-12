/******************************************************************************
 * random.h
 * *
 * Random number generation functions.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#ifndef SMHM_RANDOM_H
#define SMHM_RANDOM_H

#include <random>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
// Own headers
#include "lib/utils/definitions.h"

typedef std::mt19937 MersenneTwister;

class RandomFunctions
{
private:
    // Initialize the engine for randomness
    // See: https://cplusplus.com/reference/random/mt19937/
    static MersenneTwister randEngine;
    static int randSeed;

public:
    // Set the seed for the random number generator
    static void set_seed(int seed)
    {
        randSeed = seed;
        srand(seed);
        randEngine.seed(seed);
    }

    // Get the seed for the random number generator
    static int get_seed()
    {
        return randSeed;
    }

    // Get the random engine
    static MersenneTwister &get_random_engine()
    {
        return randEngine;
    }

    // Get a random uniformly chosen integer value the range [lb, rb]
    template <typename DataType>
    static DataType get_uniform_random_int_in_bounds(DataType lb, DataType rb)
    {
        return get_uniform_random_int_in_bounds(lb, rb, randEngine);
    }

    // Get a random uniformly chosen integer value the range [lb, rb] with an alternative random engine
    template <typename DataType>
    static DataType get_uniform_random_int_in_bounds(DataType lb, DataType rb, MersenneTwister &altRandEngine)
    {
        std::uniform_int_distribution<DataType> dist(lb, rb);
        return dist(altRandEngine);
    }

    // Get a random number
    static uint32_t get_uniform_random_int()
    {
        return randEngine();
    }

    // Get a random number with an alternative random engine
    static uint32_t get_uniform_random_int(MersenneTwister &altRandEngine)
    {
        return altRandEngine();
    }

    // Permutate a given vector by shuffling local blocks of size 128 (taken and modified from VieCut
    template <typename DataType>
    static void permutate_vector_local(std::vector<DataType> *v, bool init)
    {
        permutate_vector_local(v, init, randEngine);
    }

    // Permutate a given vector by shuffling local blocks of size 128 (taken and modified from VieCut) with an alternative random engine
    template <typename DataType>
    static void permutate_vector_local(std::vector<DataType> *v, bool init, MersenneTwister &altRandEngine)
    {
        std::vector<DataType> &vec = *v;
        // Initialize the vector to contain 0, 1, 2, ..., |vec|-1 (if necessary)
        if (init)
            std::iota(vec.begin(), vec.end(), 0);

        // Shuffle local blocks of size 128
        size_t localsize = 128;
        for (size_t i = 0; i < vec.size(); i += localsize)
        {
            size_t end = std::min(vec.size(), i + localsize);
            std::shuffle(vec.begin() + i, vec.begin() + end, altRandEngine);
        }
    }

    template <typename DataType>
    static void permutate_vector_good(std::vector<DataType> *v, bool init, MersenneTwister &altRandEngine)
    {
        std::vector<DataType> &vec = *v;
        // Initialize the vector to contain 0, 1, 2, ..., |vec|-1 (if necessary)
        if (init)
            std::iota(vec.begin(), vec.end(), 0);

        if (vec.size() < 10)
        {
            permutate_vector_good_small(v, altRandEngine);
            return;
        }
        unsigned int size = vec.size();
        std::uniform_int_distribution<unsigned int> A(0, size - 4);
        std::uniform_int_distribution<unsigned int> B(0, size - 4);

        for (unsigned int i = 0; i < size; i++)
        {
            unsigned int posA = A(altRandEngine);
            unsigned int posB = B(altRandEngine);
            std::swap(vec[posA], vec[posB]);
            std::swap(vec[posA + 1], vec[posB + 1]);
            std::swap(vec[posA + 2], vec[posB + 2]);
            std::swap(vec[posA + 3], vec[posB + 3]);
        }
    }

    template <typename DataType>
    static void permutate_vector_good_small(std::vector<DataType> *v, MersenneTwister &altRandEngine)
    {
        std::vector<DataType> &vec = *v;
        if (vec.size() < 2)
            return;
        unsigned int size = vec.size();
        std::uniform_int_distribution<unsigned int> A(0, size - 1);
        std::uniform_int_distribution<unsigned int> B(0, size - 1);

        for (unsigned int i = 0; i < size; i++)
        {
            unsigned int posA = A(altRandEngine);
            unsigned int posB = B(altRandEngine);
            std::swap(vec[posA], vec[posB]);
        }
    }
};

#endif // end of SMHM_RANDOM_H