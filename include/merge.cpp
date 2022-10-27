#include "merge.h"
#include "bucket.h"
#include "common.h"

int merge_sort(int inStructureId, int outStructureId) {
  int64_t bucketNum = ceil(1.0 * N / M);
  std::cout << "Bucket Number: " << bucketNum << std::endl;
  int64_t setSize = floor(1.0 * M / (bucketNum + 1));
  initMerge(setSize);
  std::cout << "HEAP_NODE_SIZE: " << setSize << std::endl;
  int64_t avg = 1.0*N/bucketNum;
  int64_t remainder = 1.0*N - avg*bucketNum;
  int64_t *bucketAddr = (int64_t*)malloc(bucketNum * sizeof(int64_t));
  int64_t *elements = (int64_t*)malloc(bucketNum * sizeof(int64_t));
  int64_t each, total = 0;
  printf("sizeof bucket_x: %lu\n", sizeof(Bucket_x));
  for (int64_t i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    elements[i] = each;
    bucketAddr[i] = total;
    printf("total: %lu\n", bucketAddr[i]);
    total += each;
  }
  // Bucket_x *arr = (Bucket_x*)malloc(M * sizeof(Bucket_x));
  for (int64_t i = 0; i < bucketNum; ++i) {
    printf("bucketSort %lu, elements: %lu, bucketAddr: %lu\n", i, elements[i], bucketAddr[i]);
    bucketSort(inStructureId, elements[i], bucketAddr[i]);
  }
  printf("Begin merge sort\n");
  kWayMergeSort(inStructureId, outStructureId, elements, bucketAddr, bucketNum);
  free(bucketAddr);
  free(elements);
  return outStructureId;
}
