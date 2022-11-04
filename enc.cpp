#include <iostream>
#include <cmath>
#include <cstdint>  
#include <cstdlib>
#include <stdio.h>
#include <sys/time.h>
#include <assert.h>
#include <random>
#include <chrono>
#include <utility>
#include <algorithm>
#include <iomanip>
// TODO: new
#include <vector>
#include <cstring>
#include <unordered_set>
#include <unordered_map>

#include "include/common.h"
#include "include/bucket.h"
#include "include/merge.h"
#include "include/oq.h"

// Globals
int64_t *X;
//structureId=3, write back array
int64_t *Y;
//structureId=1, bucket1 in bucket sort; input
Bucket_x *bucketx1;
//structureId=2, bucket 2 in bucket sort
Bucket_x *bucketx2;

int64_t *arrayAddr[NUM_STRUCTURES];
int64_t paddedSize;
int64_t FAN_OUT;
int64_t BUCKET_SIZE;
double IOcost = 0;
int is_tight;

/** procedure test() : verify sort results **/
void init(int64_t **arrayAddr, int structureId, int64_t size) {
  int64_t i;
  if (structureSize[structureId] == 8) {
    int64_t *addr = (int64_t*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      addr[i] = size - i;
    }
  } else if (structureSize[structureId] == 16) {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (i = 0; i < size; ++i) {
      // TODO: size overflow, become negative
      addr[i].x = size - i;
    }
  }
}

// trusted function
void callSort(int sortId, int structureId, int64_t paddedSize, int *resId, int64_t *resN) {
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
      std::pair<int, int64_t> ans = ObliviousLooseSort(structureId, paddedSize, structureId + 1, structureId);
      *resId = ans.first;
      *resN = ans.second;
    } else {
      std::pair<int, int64_t> ans = ObliviousLooseSort2(structureId, paddedSize, structureId + 1, structureId + 2, structureId + 1, structureId);
      *resId = ans.first;
      *resN = ans.second;
    }
  } else if (sortId == 2) {
    *resId = bucketOSort(structureId, paddedSize);
  } else if (sortId == 4) {
    *resId = merge_sort(structureId, structureId+1);
  } else {
    Bucket_x *test = (Bucket_x*)malloc(sizeof(Bucket_x)*N);
    Bucket_x *test2 = (Bucket_x*)malloc(sizeof(Bucket_x)*2*N);
    for (int i = 0; i < N; ++i) {
      test[i].x = N - i;
      test[i].key = N - i;
    }
    aes_init();
    freeAllocate(1, 1, 4*N);
    nonEnc = 1;
    printf("Before encrypt: \n");
    print((int64_t*)test, 2*N);
    for (int i = 0; i < 1; ++i) {
      opOneLinearScanBlock(0, (int64_t*)test, N, 1, 1, 0);
      opOneLinearScanBlock(0, (int64_t*)test2, N, 1, 0, 0);
    }
    printf("After decrypt: \n");
    print((int64_t*)test2, 2*N);
    free(test);
    free(test2);
  }
}

/* main function */
int main(int argc, const char* argv[]) {
  int ret = 1;
  int *resId = (int*)malloc(sizeof(int));
  int64_t *resN = (int64_t*)malloc(sizeof(int64_t));
  // oe_result_t result;
  // oe_enclave_t* enclave = NULL;
  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::seconds duration;
  srand((unsigned)time(NULL));
  // freopen("/home/data/bchenba/errors.txt", "w+", stdout); 

  // 0: OQSORT-Tight, 1: OQSORT-Loose, 2: bucketOSort, 3: bitonicSort, 4: merge_sort
  int sortId = 5;
  int inputId = 0;

  // double beta = calParams();
  // double alpha = (KAPPA+1+log(1.0*N))*4*(1+beta)*(1+2*beta)/beta/beta/M;
  // int p = ceil((1+2*beta)*N/M);
  // printf("Parameters: alpha: %f, beta: %f, p: %d\n", alpha, beta, p);
  // step1: init test numbers
  if (sortId == 3) {
    // inputId = 0;
    int64_t addi = 0;
    if (N % BLOCK_DATA_SIZE != 0) {
      addi = ((N / BLOCK_DATA_SIZE) + 1) * BLOCK_DATA_SIZE - N;
    }
    X = (int64_t*)malloc((N + addi) * sizeof(int64_t));
    paddedSize = N + addi;
    arrayAddr[inputId] = X;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 4) {
    inputId = 1;
    // arrayAddr[inputId] = X;
    bucketx1 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    bucketx2 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    arrayAddr[inputId] = (int64_t*)bucketx1;
    arrayAddr[inputId+1] = (int64_t*)bucketx2;
    paddedSize = N;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 2) {
    // inputId = 0;
    double z1 = 6 * (KAPPA + log(2.0*N));
    double z2 = 6 * (KAPPA + log(2.0*N/z1));
    BUCKET_SIZE = BLOCK_DATA_SIZE * ceil(1.0*z2/BLOCK_DATA_SIZE);
    std::cout << "BUCKET_SIZE: " << BUCKET_SIZE << std::endl;
    double thresh = 1.0*M/BUCKET_SIZE;
    std::cout << "Threash: " << thresh << std::endl;
    FAN_OUT = greatestPowerOfTwoLessThan(thresh)/2;
    assert(FAN_OUT >= 2 && "M/Z must greater than 2");
    int64_t bucketNum = smallestPowerOfKLargerThan(ceil(2.0 * N / BUCKET_SIZE), 2);
    int64_t bucketSize = bucketNum * BUCKET_SIZE;
    std::cout << "TOTAL BUCKET SIZE: " << bucketSize << std::endl;
    std::cout << "BUCKET NUMBER: " << bucketNum << std::endl;
    std::cout << "BUCKET SIZE: " << BUCKET_SIZE << std::endl; 
    std::cout << "FAN_OUT: " << FAN_OUT << std::endl;  
    bucketx1 = (Bucket_x*)malloc(bucketSize * sizeof(Bucket_x));
    bucketx2 = (Bucket_x*)malloc(bucketSize * sizeof(Bucket_x));
    memset(bucketx1, 0xff, bucketSize*sizeof(Bucket_x));
    memset(bucketx2, 0xff, bucketSize*sizeof(Bucket_x));
    std::cout << "After bucket malloc\n";
    arrayAddr[1] = (int64_t*)bucketx1;
    arrayAddr[2] = (int64_t*)bucketx2;
    X = (int64_t*) malloc(N * sizeof(int64_t));
    arrayAddr[inputId] = X;
    paddedSize = N;
    init(arrayAddr, inputId, paddedSize);
  } else if (sortId == 0 || sortId == 1) {
    inputId = 3;
    X = (int64_t*)malloc(N * sizeof(int64_t));
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
    if (*resId == -1) {
      std::cout << "TEST Failed\n";
    } else if (sortId == 0) {
      test(arrayAddr, *resId, paddedSize);
      *resN = N;
    } else if (sortId == 1) {
      // Sample Loose has different test & print
      testWithDummy(arrayAddr, *resId, *resN);
    }
  } else {
    callSort(sortId, inputId, paddedSize, resId, resN);
  }
  end = std::chrono::high_resolution_clock::now();
  // step4: std::cout execution time
  duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "Finished. Duration Time: " << duration.count() << " seconds" << std::endl;
  // std::cout.precision(4);
  printf("IOcost: %f\n", 1.0*IOcost/N*(BLOCK_DATA_SIZE/2));
  print(arrayAddr, *resId, *resN);
  // step5: exix part
  exit:
    for (int i = 0; i < NUM_STRUCTURES; ++i) {
      if (arrayAddr[i]) {
        // free(arrayAddr[i]);
      }
    }
    free(resId);
    free(resN);
    return ret;
}
