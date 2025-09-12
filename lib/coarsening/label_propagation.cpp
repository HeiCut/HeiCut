/******************************************************************************
 * label_propagation.cpp
 * *
 * Label propagation algorithm to coarsen a hypergraph.
 * *
 * Loris Wilwert <loris.wilwert@stud.uni-heidelberg.de>
 *****************************************************************************/

// Own headers
#include "label_propagation.h"
#include "lib/utils/definitions.h"

LabelPropagation::LabelPropagation(const IterationIndex numIterations, const LabelPropagationMode mode, const NodeIndex numPinsToSample)
    : numIterations(numIterations),
      mode(mode),
      numPinsToSample(numPinsToSample) {};