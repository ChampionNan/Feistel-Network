#ifndef BITONIC_H
#define BITONIC_H

#include "common.h"

#include <cstdlib>
#include <cstdint>

void smallBitonicMerge(int64_t *a, int64_t start, int64_t size, int64_t flipped);
void smallBitonicSort(int64_t *a, int64_t start, int64_t size, int64_t flipped);

#endif // !BITONIC_H
