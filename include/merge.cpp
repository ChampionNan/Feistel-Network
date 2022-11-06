#include "merge.h"
#include "bucket.h"
#include "common.h"

int merge_sort(int inStructureId, int outStructureId) {
  int bucketNum = ceil(1.0 * N / M);
  std::cout << "Bucket Number: " << bucketNum << std::endl;
  int setSize = floor(1.0 * M / (bucketNum + 1));
  initMerge(setSize);
  std::cout << "HEAP_NODE_SIZE: " << setSize << std::endl;
  int avg = 1.0*N/bucketNum;
  int remainder = 1.0*N - avg*bucketNum;
  int *bucketAddr = (int*)malloc(bucketNum * sizeof(int));
  int *elements = (int*)malloc(bucketNum * sizeof(int));
  int each, total = 0;
  printf("sizeof bucket_x: %lu\n", sizeof(Bucket_x));
  for (int i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    elements[i] = each;
    bucketAddr[i] = total;
    printf("total: %lu\n", bucketAddr[i]);
    total += each;
  }
  // Bucket_x *arr = (Bucket_x*)malloc(M * sizeof(Bucket_x));
  for (int i = 0; i < bucketNum; ++i) {
    printf("bucketSort %lu, elements: %lu, bucketAddr: %lu\n", i, elements[i], bucketAddr[i]);
    bucketSort(inStructureId, elements[i], bucketAddr[i]);
  }
  printf("Begin merge sort\n");
  kWayMergeSort(inStructureId, outStructureId, elements, bucketAddr, bucketNum);
  free(bucketAddr);
  free(elements);
  return outStructureId;
}
