#ifndef BUCKET_H

#include "common.h"

void initMerge(int64_t size);
bool isTargetIterK(int64_t randomKey, int64_t iter, int64_t k, int64_t num);
void mergeSplitHelper(Bucket_x *inputBuffer, int64_t* numRow1, int64_t* numRow2, int64_t* inputId, int64_t* outputId, int64_t iter, int64_t k, int64_t* bucketAddr, int64_t outputStructureId);
void mergeSplit(int64_t inputStructureId, int64_t outputStructureId, int64_t *inputId, int64_t *outputId, int64_t k, int64_t* bucketAddr, int64_t* numRow1, int64_t* numRow2, int64_t iter);
void kWayMergeSort(int64_t inputStructureId, int64_t outputStructureId, int64_t* numRow1, int64_t* bucketAddr, int64_t bucketNum);
void bucketSort(int64_t inputStructureId, int64_t size, int64_t dataStart);
int64_t bucketOSort(int64_t structureId, int64_t size);

#endif // !BUCKET_H
