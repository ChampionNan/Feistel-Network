#include <cmath>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <assert.h>
#include <bitset>
#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include <iomanip>
// TODO: new
#include <unordered_set>
#include <unordered_map>
#include <bitset>

#include "include/common.h"
#include "include/bucket.h"
#include "include/merge.h"
#include "include/oq.h"

// OCALL
void ocall_print_string(const char *str);
void OcallReadBlock(int index, int* buffer, size_t blockSize, int structureId);
void OcallWriteBlock(int index, int* buffer, size_t blockSize, int structureId);
void freeAllocate(int structureIdM, int structureIdF, int size);

// SUPPORT
void callSort(int sortId, int structureId, int paddedSize, int *resId, int *resN);

void init(int **arrayAddr, int structurId, int size);
void print(int* array, int size);
void print(int **arrayAddr, int structureId, int size);
void test(int **arrayAddr, int structureId, int size);
void testWithDummy(int **arrayAddr, int structureId, int size);
int greatestPowerOfTwoLessThan(double n);
int smallestPowerOfKLargerThan(int n, int k);
bool isTargetIterK(int randomKey, int iter, int k, int num);

int *X;
//structureId=3, write back array
int *Y;
//structureId=1, bucket1 in bucket sort; input
Bucket_x *bucketx1;
//structureId=2, bucket 2 in bucket sort
Bucket_x *bucketx2;
int *arrayAddr[NUM_STRUCTURES];
int paddedSize;
int FAN_OUT;
int BUCKET_SIZE;

double IOcost = 0;
int is_tight;

/** procedure test() : verify sort results **/
void init(int **arrayAddr, int structureId, int size) {
  int i;
  if (structureSize[structureId] == 4) {
    int *addr = (int*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      addr[i] = (size - i);
    }
  } else if (structureSize[structureId] == 8) {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (int i = 0; i < size; ++i) {
      addr[i].x = size - i;
    }
  }
}

/* main function */
int main(int argc, const char* argv[]) {
  int ret = 1;
  int *resId = (int*)malloc(sizeof(int));
  int *resN = (int*)malloc(sizeof(int));
  // oe_result_t result;
  // oe_enclave_t* enclave = NULL;
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::seconds duration;
  srand((unsigned)time(NULL));
  
  // 0: OSORT-Tight, 1: OSORT-Loose, 2: bucketOSort, 3: bitonicSort, 4: merge_sort
  int sortId = 2;
  int inputId = 0;

  // step1: init test numbers
  if (sortId == 3) {
    // inputId = 0;
    int addi = 0;
    if (N % BLOCK_DATA_SIZE != 0) {
      addi = ((N / BLOCK_DATA_SIZE) + 1) * BLOCK_DATA_SIZE - N;
    }
    X = (int*)malloc((N + addi) * sizeof(int));
    paddedSize = N + addi;
    arrayAddr[inputId] = X;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 4) {
    inputId = 1;
    arrayAddr[inputId] = X;
    bucketx1 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    bucketx2 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    arrayAddr[inputId] = (int*)bucketx1;
    arrayAddr[inputId+1] = (int*)bucketx2;
    paddedSize = N;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 2) {
    // inputId = 0;
    double z1 = 6 * (KAPPA + log(2*N));
    double z2 = 6 * (KAPPA + log(2.0*N/z1));
    BUCKET_SIZE = BLOCK_DATA_SIZE * ceil(1.0*z2/BLOCK_DATA_SIZE);
    std::cout << "BUCKET_SIZE: " << BUCKET_SIZE << std::endl;
    double thresh = 1.0*M/BUCKET_SIZE;
    std::cout << "Threash: " << thresh << std::endl;
    FAN_OUT = greatestPowerOfTwoLessThan(thresh);
    assert(FAN_OUT >= 2 && "M/Z must greater than 2");
    int bucketNum = smallestPowerOfKLargerThan(ceil(2.0 * N / BUCKET_SIZE), 2);
    int bucketSize = bucketNum * BUCKET_SIZE;
    std::cout << "TOTAL BUCKET SIZE: " << bucketSize << std::endl;
    std::cout << "BUCKET NUMBER: " << bucketNum << std::endl;
    std::cout << "BUCKET SIZE: " << BUCKET_SIZE << std::endl; 
    std::cout << "FAN_OUT: " << FAN_OUT << std::endl;  
    bucketx1 = (Bucket_x*)malloc(bucketSize * sizeof(Bucket_x));
    bucketx2 = (Bucket_x*)malloc(bucketSize * sizeof(Bucket_x));
    memset(bucketx1, 0xff, bucketSize*sizeof(Bucket_x));
    memset(bucketx2, 0xff, bucketSize*sizeof(Bucket_x));
    arrayAddr[1] = (int*)bucketx1;
    arrayAddr[2] = (int*)bucketx2;
    X = (int *) malloc(N * sizeof(int));
    arrayAddr[inputId] = X;
    paddedSize = N;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 0 || sortId == 1) {
    inputId = 3;
    X = (int *)malloc(N * sizeof(int));
    arrayAddr[inputId] = X;
    paddedSize = N;
    init(arrayAddr, inputId, paddedSize);
  }

  // step2: Create the enclave
  // print(arrayAddr, inputId, N);
  
  // step3: call sort algorithms
  start = std::chrono::high_resolution_clock::now();
  if (sortId == 3) {
    std::cout << "Test bitonic sort... " << std::endl;
    callSort(sortId, inputId, paddedSize, resId, resN);
    test(arrayAddr, inputId, paddedSize);
  } else if (sortId == 4) {
    std::cout << "Test merge_sort... " << std::endl;
    callSort(sortId, inputId, paddedSize, resId, resN);
    std::cout << "Result ID: " << *resId << std::endl;
    *resN = N;
    test(arrayAddr, *resId, paddedSize);
  } else if (sortId == 2) {
    std::cout << "Test bucket oblivious sort... " << std::endl;
    callSort(sortId, inputId + 1, paddedSize, resId, resN);
    std::cout << "Result ID: " << *resId << std::endl;
    *resN = N;
    // print(arrayAddr, *resId, N);
    test(arrayAddr, *resId, paddedSize);
  } else if (sortId == 0 || sortId == 1) {
    std::cout << "Test OQSort... " << std::endl;
    callSort(sortId, inputId, paddedSize, resId, resN);
    std::cout << "Result ID: " << *resId << std::endl;
    if (sortId == 0) {
      test(arrayAddr, *resId, paddedSize);
      *resN = N;
    } else {
      // Sample Loose has different test & print
      testWithDummy(arrayAddr, *resId, *resN);
    }
  }
  end = std::chrono::high_resolution_clock::now();
  // step4: std::cout execution time
  duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "Finished. Duration Time: " << duration.count() << " seconds" << std::endl;
  std::cout.precision(4);
  printf("IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  print(arrayAddr, *resId, *resN);
  // step5: exix part
  exit:
    for (int i = 0; i < NUM_STRUCTURES; ++i) {
      if (arrayAddr[i]) {
        free(arrayAddr[i]);
      }
    }
    free(resId);
    free(resN);
    return ret;
}

// trusted function
void callSort(int sortId, int structureId, int paddedSize, int *resId, int *resN) {
  // TODO: Utilize Memory alloction -- structureId
  if (sortId == 0) {
    is_tight = 1;
    if (paddedSize / M <= 128) {
      *resId = ObliviousTightSort(structureId, paddedSize, structureId + 1, structureId);
    } else {
      *resId = ObliviousTightSort2(structureId, paddedSize, structureId+1, structureId+2, structureId+1, structureId);
    }
  } else if (sortId == 1) {
    if (paddedSize / M <= 128) {
      is_tight = 0;
      std::pair<int, int> ans = ObliviousLooseSort(structureId, paddedSize, structureId + 1, structureId);
      *resId = ans.first;
      *resN = ans.second;
    } else {
      std::pair<int, int> ans = ObliviousLooseSort2(structureId, paddedSize, structureId + 1, structureId + 2, structureId + 1, structureId);
      *resId = ans.first;
      *resN = ans.second;
    }
  } else if (sortId == 2) {
    *resId = bucketOSort(structureId, paddedSize);
  } else if (sortId == 4) {
    *resId = merge_sort(structureId, structureId+1, paddedSize);
  }
}

