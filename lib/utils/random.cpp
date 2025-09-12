/******************************************************************************
 * random.cpp
 * *
 * Random number generation functions.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

#include "random.h"
// Own headers
#include "lib/utils/const.h"

MersenneTwister RandomFunctions::randEngine;
int RandomFunctions::randSeed = DEFAULT_SEED;