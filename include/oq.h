#ifndef OQ_H
#define OQ_H

#include "common.h"

#include <cstdlib>
#include <cstdint>

// others
__uint128_t prf(__uint128_t a);
int64_t encrypt(int64_t index);
void pseudo_init(int64_t size);
int64_t sumArray(bool *M, int64_t left, int64_t right);
void swapArray(int64_t *D, int64_t i, int64_t j);
void OROffCompact(int64_t *D, bool *M, int64_t left, int64_t right, int64_t z);
void ORCompact(int64_t *D, bool *M, int64_t left, int64_t right);
int64_t assignM(int64_t *arr, bool *M, int64_t left, int64_t right, int64_t pivot);
void obliviousPWayPartition(int64_t *D, bool *M, int64_t low, int64_t high, std::vector<int64_t> pivots, int64_t left, int64_t right, std::vector<int64_t> &partitionIdx);
void internalObliviousSort(int64_t *D, int64_t left, int64_t right, int64_t smallM, double beta, double gamma);
// main sorting part
void floydSampler(int64_t n, int64_t k, std::vector<int64_t> &x);
int64_t Sample(int64_t inStructureId, int64_t sampleSize, std::vector<int64_t> &trustedM2, int64_t is_tight, int64_t is_rec);
void SampleRec(int64_t inStructureId, int64_t sampleId, int64_t sortedSampleId, int64_t is_tight, std::vector<std::vector<int64_t> >& pivots);
void quantileCal(std::vector<int64_t> &samples, int64_t start, int64_t end, int64_t p);
void quantileCal2(std::vector<int64_t> &samples, int64_t start, int64_t end, int64_t p);
int64_t partitionMulti(int64_t *arr, int64_t low, int64_t high, int64_t pivot);
void quickSortMulti(int64_t *arr, int64_t low, int64_t high, std::vector<int64_t> pivots, int64_t left, int64_t right, std::vector<int64_t> &partitionIdx);
std::pair<int64_t, int64_t> OneLevelPartition(int64_t inStructureId, int64_t inSize, std::vector<int64_t> &samples, int64_t sampleSize, int64_t p, int64_t outStructureId1, int64_t is_rec);
std::pair<int64_t, int64_t> TwoLevelPartition(int64_t inStructureId, std::vector<std::vector<int64_t> >& pivots, int64_t p, int64_t outStructureId1, int64_t outStructureId2);
int64_t ObliviousTightSort(int64_t inStructureId, int64_t inSize, int64_t outStructureId1, int64_t outStructureId2);
int64_t ObliviousTightSort2(int64_t inStructureId, int64_t inSize, int64_t sampleId, int64_t sortedSampleId, int64_t outStructureId, int64_t outStructureId2);
std::pair<int64_t, int64_t> ObliviousLooseSort(int64_t inStructureId, int64_t inSize, int64_t outStructureId1, int64_t outStructureId2);
std::pair<int64_t, int64_t> ObliviousLooseSort2(int64_t inStructureId, int64_t inSize, int64_t sampleId, int64_t sortedSampleId, int64_t outStructureId1, int64_t outStructureId2);
void ObliviousLooseSortRec(int64_t sampleId, int64_t sampleSize, int64_t sortedSampleId, std::vector<std::vector<int64_t> >& pivots);

#endif // !OQ_H
