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
int *X;
//structureId=3, write back array
int *Y;
//structureId=1, bucket1 in bucket sort; input
Bucket_x *bucketx1;
//structureId=2, bucket 2 in bucket sort
Bucket_x *bucketx2;

int *arrayAddr[NUM_STRUCTURES];
int paddedSize;
double IOcost;
int is_tight;

/** procedure test() : verify sort results **/
void init(int **arrayAddr, int structureId, int size) {
  int i;
  if (structureSize[structureId] == 4) {
    int *addr = (int*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      addr[i] = size - i;
    }
  } else if (structureSize[structureId] == 8) {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (i = 0; i < size; ++i) {
      // TODO: size overflow, become negative
      addr[i].x = size - i;
    }
  }
}

// size: real numbers, no encrypt but still EncB structure
void initEnc(int **arrayAddr, int structureId, int size) {
  int i = 0, j, blockNumber;
  EncBlock *addr = (EncBlock*)arrayAddr[structureId];
  if (structureSize[structureId] == 4) {
    blockNumber = (int)ceil(1.0*size/BLOCK_DATA_SIZE);
    for (i = 0; i < blockNumber; i++) {
      int *addx = (int*)(&addr[i]);
      for (j = 0; j < BLOCK_DATA_SIZE; ++j) {
        addx[j] = size - (i * BLOCK_DATA_SIZE + j);
      }
    }
  } else if (structureSize[structureId] == 8) {
    blockNumber = (int)ceil(1.0*size/BLOCK_DATA_SIZE/2); 
    for (i = 0; i < blockNumber; ++i) {
      Bucket_x *addx = (Bucket_x*)(&addr[i]);
      for (j = 0; j < BLOCK_DATA_SIZE/2; ++j) {
        addx[j].x = size - (i * BLOCK_DATA_SIZE/2 + j);
      }
    }
  }
}

void testIO(int IOnum, int inId) {
  int totalB = ceil(IOnum / 2);
  printf("Testing IO: %d\n", IOnum);
  aes_init();
  int *buffer = (int*)malloc(PAGEBYTES);
  int pageNum = PAGEBYTES / 4;
  freeAllocate(inId, inId, totalB*PAGEBYTES/4);
  // std::random_device dev;
  // std::mt19937 rng(dev());
  // std::uniform_int_distribution<int> dist{0, totalB-1};
  int index;
  nonEnc = 0;
  for (int i = 0; i < totalB; ++i) {
    // index = dist(rng);
    opOneLinearScanBlock(i*pageNum, buffer, PAGEBYTES/8, inId, 0, 0);
  }
  for (int i = 0; i < totalB; ++i) {
    // index = dist(rng);
    opOneLinearScanBlock(i*pageNum, buffer, PAGEBYTES/8, inId, 1, 0);
  }
  free(buffer);
}

// trusted function
void callSort(int sortId, int structureId, int paddedSize, int *resId, int *resN) {
  // TODO: Utilize Memory alloction -- structureId
  if (sortId == 0) {
    is_tight = 1;
    if (paddedSize / M <= 128) {
      *resId = ObliviousTightSort(structureId, paddedSize, structureId + 1, structureId);
    }
  } else if (sortId == 1) {
    if (paddedSize / M <= 128) {
      is_tight = 0;
      std::pair<int, int> ans = ObliviousLooseSort(structureId, paddedSize, structureId + 1, structureId);
      *resId = ans.first;
      *resN = ans.second;
    }
  } else if (sortId == 2) {
    *resId = bucketOSort(structureId, paddedSize);
  } else if (sortId == 4) {
    *resId = merge_sort(structureId, structureId+1);
  } else if (sortId == 5){
    testIO(2, structureId);
  } else {
    std::cout << "In testing: \n";
    aes_init();
    // 0, 3, 4: int , 1, 2: bucket
    freeAllocate(0, 0, ceil(1.0*N/BLOCK_DATA_SIZE)); 
    freeAllocate(3, 3, ceil(1.0*N/BLOCK_DATA_SIZE));
    // freeAllocate(1, 1, ceil(1.0*N/2));
    // freeAllocate(2, 2, ceil(1.0*N/2));
    initEnc(arrayAddr, 0, N);
    // initEnc(arrayAddr, 1, N); 
    printf("Initial\n");
    printEnc(arrayAddr, 0, N);
    int *read = (int*)malloc(N*sizeof(int));
    // int *read = (int*)malloc(N*sizeof(Bucket_x));
    nonEnc = 1;
    opOneLinearScanBlock(0, read, N, 0, 0, 0);
    std::cout << "After Read in nonEnc\n";
    print(read,N); // pass
    nonEnc = 0;
    printf("Before write Enc\n");
    // write & Enc
    opOneLinearScanBlock(0, read, N, 3, 1, 0);
    std::cout << "After write in Enc\n";
    printEnc(arrayAddr, 3, N);
    // read & decrypt
    opOneLinearScanBlock(0, read, N, 3, 0, 0);
    std::cout << "After Read in Enc\n";
    print(read, N);
    nonEnc = 1;
    opOneLinearScanBlock(0, read, N, 0, 1, 0);
    std::cout << "After write in Enc2\n";
    printEnc(arrayAddr, 0, N);
    opOneLinearScanBlock(0, read, N, 0, 0, 0);
    std::cout << "After Read in Enc2\n";
    print(read, N);
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
  IOcost = 0;
  // freopen("/home/data/bchenba/errors.txt", "w+", stdout); 

  // 0: OQSORT-Tight, 1: OQSORT-Loose, 2: bucketOSort, 3: bitonicSort, 4: merge_sort(x)
  int sortId = 2;
  int inputId = 0;

  // double beta = calParams();
  // double alpha = (KAPPA+1+log(1.0*N))*4*(1+beta)*(1+2*beta)/beta/beta/M;
  // int p = ceil((1+2*beta)*N/M);
  // printf("Parameters: alpha: %f, beta: %f, p: %d\n", alpha, beta, p);
  // step1: init test numbers
  int FAN_OUT;
  int BUCKET_SIZE;
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
    // arrayAddr[inputId] = X;
    bucketx1 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    bucketx2 = (Bucket_x*)malloc(N * sizeof(Bucket_x));
    arrayAddr[inputId] = (int*)bucketx1;
    arrayAddr[inputId+1] = (int*)bucketx2;
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
    int bucketNum = smallestPowerOfKLargerThan(ceil(2.0 * N / BUCKET_SIZE), 2);
    int bucketSize = bucketNum * BUCKET_SIZE;
    std::cout << "TOTAL BUCKET SIZE: " << bucketSize << std::endl;
    std::cout << "BUCKET NUMBER: " << bucketNum << std::endl;
    std::cout << "BUCKET SIZE: " << BUCKET_SIZE << std::endl; 
    std::cout << "FAN_OUT: " << FAN_OUT << std::endl;  
    freeAllocate(1, 1, ceil(1.0*bucketSize/2));
    freeAllocate(2, 2, ceil(1.0*bucketSize/2));
    std::cout << "After bucket malloc\n";
    freeAllocate(inputId, inputId, ceil(1.0*N/4));
    paddedSize = N;
    // TODO: 
    initEnc(arrayAddr, inputId, paddedSize);
  } else if (sortId == 0 || sortId == 1) {
    inputId = 3;
    freeAllocate(inputId, inputId, ceil(1.0*N/4));
    paddedSize = N;
    initEnc(arrayAddr, inputId, paddedSize);
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
    testEnc(arrayAddr, *resId, paddedSize);
  } else if (sortId == 0 || sortId == 1) {
    std::cout << "Test OQSort... " << std::endl;
    callSort(sortId, inputId, paddedSize, resId, resN);
    std::cout << "Result ID: " << *resId << std::endl;
    if (*resId == -1) {
      std::cout << "TEST Failed\n";
    } else if (sortId == 0) {
      // test(arrayAddr, *resId, paddedSize);
      *resN = N;
    } else if (sortId == 1) {
      // Sample Loose has different test & print
      // testWithDummy(arrayAddr, *resId, *resN);
    }
  } else {
    callSort(sortId, inputId, paddedSize, resId, resN);
  }
  end = std::chrono::high_resolution_clock::now();
  // step4: std::cout execution time
  duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
  std::cout << "Finished. Duration Time: " << duration.count() << " seconds" << std::endl;
  // std::cout.precision(4);
  int multi = (sortId == 2 || sortId == 4) ? 2 : 1;
  printf("IOcost: %f, %f\n", 1.0*IOcost/N*(BLOCK_DATA_SIZE/multi), IOcost);
  printEnc(arrayAddr, *resId, *resN);
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
