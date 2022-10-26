#include "bucket.h"
#include "common.h"

// int finalFlag = 0;
int oldK;
int HEAP_NODE_SIZE;
int MERGE_BATCH_SIZE;
extern int FAN_OUT;
extern int BUCKET_SIZE;
int merSHflag1 = 0;
int merSHflag2 = 0;
std::random_device dev3;
std::mt19937 rng3(dev3());
// BUCKET_SORT
void initMerge(int size) {
  HEAP_NODE_SIZE = size;
  MERGE_BATCH_SIZE = size;
}

bool isTargetIterK(int randomKey, int iter, int k, int num) {
  /*
  while (iter) {
    randomKey = randomKey / oldK;
    iter--;
  }*/
  // return (randomKey & (0x01 << (iter - 1))) == 0 ? false : true;
  randomKey = (randomKey >> iter * oldK);
  return (randomKey % k) == num;
}

void mergeSplitHelper(Bucket_x *inputBuffer, int* numRow1, int* numRow2, int* inputId, int* outputId, int iter, int k, int* bucketAddr, int outputStructureId) {
  // int batchSize = BUCKET_SIZE; // 8192
  // TODO: FREE these malloc
  printf("In mergesplitHep!\n");
  Bucket_x **buf = (Bucket_x**)malloc(k * sizeof(Bucket_x*));
  for (int i = 0; i < k; ++i) {
    buf[i] = (Bucket_x*)malloc(BUCKET_SIZE * sizeof(Bucket_x));
  }
  
  // int counter0 = 0, counter1 = 0;
  int randomKey;
  int *counter = (int*)malloc(k * sizeof(int));
  memset(counter, 0, k * sizeof(int));
  printf("Before assign bins!\n");
  merSHflag1 = 1;
  for (int i = 0; i < k * BUCKET_SIZE; ++i) {
    if (inputBuffer[i].x != DUMMY) {
      randomKey = inputBuffer[i].key;
      for (int j = 0; j < k; ++j) {
        if (isTargetIterK(randomKey, iter, k, j)) {
          buf[j][counter[j] % BUCKET_SIZE] = inputBuffer[i];
          counter[j]++;
          // std::cout << "couter j: " << counter[j] << std::endl;
          if (counter[j] % BUCKET_SIZE == 0) {
            opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] +  numRow2[outputId[j]]), (int*)buf[j], (size_t)BUCKET_SIZE, outputStructureId, 1, 0);
            numRow2[outputId[j]] += BUCKET_SIZE;
          }
        }
      }
    }
  }
  merSHflag1 = 0;
  printf("Before final assign!\n");
  merSHflag2 = 1;
  for (int j = 0; j < k; ++j) {
    if (numRow2[outputId[j]] + (counter[j] % BUCKET_SIZE) > BUCKET_SIZE) {
      printf("overflow error during merge split!\n");
    }
    opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] + numRow2[outputId[j]]), (int*)buf[j], (size_t)(counter[j] % BUCKET_SIZE), outputStructureId, 1, BUCKET_SIZE-numRow2[outputId[j]]-(counter[j] % BUCKET_SIZE));
    numRow2[outputId[j]] += counter[j] % BUCKET_SIZE;
    // padWithDummy(outputStructureId, bucketAddr[outputId[j]], numRow2[outputId[j]], BUCKET_SIZE);
    free(buf[j]);
  }
  free(counter);
  free(buf);
  merSHflag1 = 0;
  printf("Finish mergesplitHep!\n");
}

void mergeSplit(int inputStructureId, int outputStructureId, int *inputId, int *outputId, int k, int* bucketAddr, int* numRow1, int* numRow2, int iter) {
  // step1. Read k buckets together
  printf("In mergesplit!\n");
  Bucket_x *inputBuffer = (Bucket_x*)malloc(k * sizeof(Bucket_x) * BUCKET_SIZE);
  for (int i = 0; i < k; ++i) {
    opOneLinearScanBlock(2 * bucketAddr[inputId[i]], (int*)(&inputBuffer[i * BUCKET_SIZE]), BUCKET_SIZE, inputStructureId, 0, 0);
  }
  // step2. process k buckets
  mergeSplitHelper(inputBuffer, numRow1, numRow2, inputId, outputId, iter, k, bucketAddr, outputStructureId);
  free(inputBuffer);
  printf("Finish mergesplit!\n");
}

void kWayMergeSort(int inputStructureId, int outputStructureId, int* numRow1, int* bucketAddr, int bucketNum) {
  int mergeSortBatchSize = HEAP_NODE_SIZE; // 256
  int writeBufferSize = MERGE_BATCH_SIZE; // 8192
  int numWays = bucketNum;
  printf("mergeSortBatchSize %d, %d\n", mergeSortBatchSize, writeBufferSize);
  // HeapNode inputHeapNodeArr[numWays];
  HeapNode *inputHeapNodeArr = (HeapNode*)malloc(numWays * sizeof(HeapNode));
  int totalCounter = 0;
  
  int *readBucketAddr = (int*)malloc(sizeof(int) * numWays);
  memcpy(readBucketAddr, bucketAddr, sizeof(int) * numWays);
  int writeBucketAddr = 0;
  int j = 0;
  printf("Initialize\n");
  for (int i = 0; i < numWays; ++i) {
    // TODO: 数据0跳过
    if (numRow1[i] == 0) {
      continue;
    }
    HeapNode node;
    node.data = (Bucket_x*)malloc(mergeSortBatchSize * sizeof(Bucket_x));
    node.bucketIdx = i;
    node.elemIdx = 0;
    opOneLinearScanBlock(2 * readBucketAddr[i], (int*)node.data, (size_t)std::min(mergeSortBatchSize, numRow1[i]), inputStructureId, 0, 0);
    inputHeapNodeArr[j++] = node;
    readBucketAddr[i] += std::min(mergeSortBatchSize, numRow1[i]);
  }
  
  Heap heap(inputHeapNodeArr, j, mergeSortBatchSize);
  Bucket_x *writeBuffer = (Bucket_x*)malloc(writeBufferSize * sizeof(Bucket_x));
  int writeBufferCounter = 0;
  printf("Start merging\n");
  while (1) {
    HeapNode *temp = heap.getRoot();
    memcpy(writeBuffer + writeBufferCounter, temp->data + temp->elemIdx % mergeSortBatchSize, sizeof(Bucket_x));
    writeBufferCounter ++;
    totalCounter ++;
    temp->elemIdx ++;
    
    if (writeBufferCounter == writeBufferSize) {
      opOneLinearScanBlock(2 * writeBucketAddr, (int*)writeBuffer, (size_t)writeBufferSize, outputStructureId, 1, 0);
      writeBucketAddr += writeBufferSize;
      // numRow2[temp->bucketIdx] += writeBufferSize;
      writeBufferCounter = 0;
      // print(arrayAddr, outputStructureId, numWays * BUCKET_SIZE);
    }
    
    if (temp->elemIdx < numRow1[temp->bucketIdx] && (temp->elemIdx % mergeSortBatchSize) == 0) {
      opOneLinearScanBlock(2 * readBucketAddr[temp->bucketIdx], (int*)(temp->data), (size_t)std::min(mergeSortBatchSize, numRow1[temp->bucketIdx]-temp->elemIdx), inputStructureId, 0, 0);
      
      readBucketAddr[temp->bucketIdx] += std::min(mergeSortBatchSize, numRow1[temp->bucketIdx]-temp->elemIdx);
      heap.Heapify(0);
      
    } else if (temp->elemIdx >= numRow1[temp->bucketIdx]) {
      bool res = heap.reduceSizeByOne();
      if (!res) {
        break;
      }
    } else {
      heap.Heapify(0);
    }
  }
  opOneLinearScanBlock(2 * writeBucketAddr, (int*)writeBuffer, (size_t)writeBufferCounter, outputStructureId, 1, 0);
  // numRow2[0] += writeBufferCounter;
  // TODO: ERROR writeBuffer
  free(writeBuffer);
  free(readBucketAddr);
  free(inputHeapNodeArr);
}

void bucketSort(int inputStructureId, int size, int dataStart) {
  Bucket_x *arr = (Bucket_x*)malloc(size * sizeof(Bucket_x));
  opOneLinearScanBlock(2 * dataStart, (int*)arr, (size_t)size, inputStructureId, 0, 0);
  quickSort(arr, 0, size - 1);
  opOneLinearScanBlock(2 * dataStart, (int*)arr, (size_t)size, inputStructureId, 1, 0);
  free(arr);
}

// int inputTrustMemory[BLOCK_DATA_SIZE];
int bucketOSort(int structureId, int size) {
  int k = FAN_OUT;
  int bucketNum = smallestPowerOfKLargerThan(ceil(2.0 * size / BUCKET_SIZE), 2);
  int mem1 = k * BUCKET_SIZE;
  HEAP_NODE_SIZE = fmax(floor(1.0 * M / (bucketNum+1)), 1);
  MERGE_BATCH_SIZE = HEAP_NODE_SIZE;
  std::cout << "Mem1's memory: " << mem1 << std::endl;
  std::cout << "Mem2's memory: " << HEAP_NODE_SIZE*(bucketNum+1) << std::endl;
  if (mem1 > M || MERGE_BATCH_SIZE > M) {
    std::cout << "Memory exceed\n";
  }
  // TODO: Change to ceil, not interger divide
  int ranBinAssignIters = ceil(1.0*log2(bucketNum)/log2(k));
  std::cout << "Iteration times: " << ceil(1.0*log2(bucketNum)/log2(k)) << std::endl;
  srand((unsigned)time(NULL));
  // std::cout << "Iters:" << ranBinAssignIters << std::endl;
  int *bucketAddr = (int*)malloc(bucketNum * sizeof(int));
  for (int i = 0; i < bucketNum; ++i) {
    bucketAddr[i] = i * BUCKET_SIZE;
  }
  int *numRow1 = (int*)malloc(bucketNum * sizeof(int));
  memset(numRow1, 0, bucketNum * sizeof(int));
  int *numRow2 = (int*)malloc(bucketNum * sizeof(int));
  memset(numRow2, 0, bucketNum * sizeof(int));
  
  Bucket_x *trustedMemory = (Bucket_x*)malloc(BUCKET_SIZE * sizeof(Bucket_x));
  int *inputTrustMemory = (int*)malloc(BUCKET_SIZE * sizeof(int));
  int total = 0, readStart = 0;
  int each;
  int avg = size / bucketNum;
  int remainder = size % bucketNum;
  std::uniform_int_distribution<int> dist{0, bucketNum-1};

  for (int i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    opOneLinearScanBlock(readStart, inputTrustMemory, each, structureId - 1, 0, 0);
    readStart += each;
    int randomKey;
    for (int j = 0; j < each; ++j) {
      // randomKey = (int)oe_rdrand();
      // randomKey = rand();
      randomKey = dist(rng3);
      trustedMemory[j].x = inputTrustMemory[j];
      trustedMemory[j].key = randomKey;
    }
    opOneLinearScanBlock(2 * bucketAddr[i], (int*)trustedMemory, each, structureId, 1, BUCKET_SIZE-each);
    numRow1[i] += each;
  }

  free(trustedMemory);
  free(inputTrustMemory);
  /*
  for (int i = 0; i < bucketNum; ++i) {
    //printf("currently bucket %d has %d records/%d\n", i, numRow1[i], BUCKET_SIZE);
    padWithDummy(structureId, bucketAddr[i], numRow1[i], BUCKET_SIZE);
  }*/
  // print(arrayAddr, structureId, bucketNum * BUCKET_SIZE);
  // std::cout << "Iters:" << ranBinAssignIters << std::endl;
  // std::cout << "k:" << k << std::endl;
  int *inputId = (int*)malloc(k * sizeof(int));
  int *outputId = (int*)malloc(k *sizeof(int));
  int outIdx = 0, expo, tempk, jboundary, jjboundary;
  // TODO: change k to tempk, i != last level, k = FAN_OUT
  // last level k = k % k1, 2N/Z < 2^k, 2^k1 < M/Z
  int k1 = (int)log2(k);
  // oldK = (int)pow(2, k1);
  oldK = k1;
  int tempk_i, tempk_i1;
  printf("Before RBA!\n");
  printf("Total bits: %d, each level solved %d bits\n", (int)log2(bucketNum), k1);
  for (int i = 0; i < ranBinAssignIters; ++i) {
    expo = std::min(k1, (int)log2(bucketNum)-i*k1);
    tempk = (int)pow(2, expo);
    tempk_i = (int)pow(tempk, i);
    tempk_i1 = tempk_i * tempk;
    printf("level%d read&write %d buckets\n", i, tempk);
    jboundary = (i != (ranBinAssignIters-1)) ? (bucketNum / tempk_i1) : 1;
    // jboundary = bucketNum / tempk_i1; 
    jjboundary = (i != (ranBinAssignIters-1)) ? tempk_i : (bucketNum/tempk);
    printf("jboundary: %d, jjboundary: %d\n", jboundary, jjboundary);
    // finalFlag = (i != (ranBinAssignIters-1)) ? 0 : 1;
    if (i % 2 == 0) {
      for (int j = 0; j < jboundary; ++j) {
        // pass (i-1) * k^i
        //printf("j: %d\n", j);
        for (int jj = 0; jj < jjboundary; ++jj) {
          //printf("jj: %d\n", jj);
          for (int m = 0; m < tempk; ++m) {
            //printf("j, jj, m: %d, %d, %d\n", j, jj, m);
            inputId[m] = j * tempk_i1+ jj + m * jjboundary;
            outputId[m] = (outIdx * tempk + m) % bucketNum;
            // printf("input, output: %d, %d\n", inputId[m], outputId[m]);
          }/*
          for (int m = 0; m < tempk; ++m) {
            printf("input%d: %d\n", m, inputId[m]);
          }
          for (int m = 0; m < tempk; ++m) {
            printf("output%d: %d\n", m, outputId[m]);
          }*/
          printf("Before mergeSplit\n");
          mergeSplit(structureId, structureId + 1, inputId, outputId, tempk, bucketAddr, numRow1, numRow2, i);
          outIdx ++;
        }
      }
      int count = 0;
      for (int n = 0; n < bucketNum; ++n) {
        numRow1[n] = 0;
        count += numRow2[n];
      }
      printf("after %dth merge split, we have %d tuples\n", i, count);
      // testRealNum(structureId + 1);
      // cnt.clear();
      outIdx = 0;
      //print(arrayAddr, structureId + 1, bucketNum * BUCKET_SIZE);
    } else {
      for (int j = 0; j < jboundary; ++j) {
        //printf("j: %d\n", j);
        for (int jj = 0; jj < jjboundary; ++jj) {
          //printf("jj: %d\n", jj);
          for (int m = 0; m < tempk; ++m) {
            //printf("j, jj, m: %d, %d, %d\n", j, jj, m);
            inputId[m] = j * tempk_i1+ jj + m * jjboundary;
            outputId[m] = (outIdx * tempk + m) % bucketNum;
            //printf("input, output: %d, %d\n", inputId[m], outputId[m]);
          }/*
          for (int m = 0; m < tempk; ++m) {
            printf("input%d: %d\n", m, inputId[m]);
          }
          for (int m = 0; m < tempk; ++m) {
            printf("output%d: %d\n", m, outputId[m]);
          }*/
          printf("Before mergeSplit\n");
          mergeSplit(structureId + 1, structureId, inputId, outputId, tempk, bucketAddr, numRow2, numRow1, i);
          outIdx ++;
        }
      }
      int count = 0;
      for (int n = 0; n < bucketNum; ++n) {
        numRow2[n] = 0;
        count += numRow1[n];
      }
      printf("after %dth merge split, we have %d tuples\n", i, count);
      // testRealNum(structureId);
      // cnt.clear();
      outIdx = 0;
      //print(arrayAddr, structureId, bucketNum * BUCKET_SIZE);
    }
    std::cout << "----------------------------------------\n";
    printf("\n\n Finish random bin assignment iter%dth out of %d\n\n", i, ranBinAssignIters);
    std::cout << "----------------------------------------\n";
  }
  free(inputId);
  free(outputId);
  int resultId = 0;
  printf("Finish RBA!\n");
  if (ranBinAssignIters % 2 == 0) {
    print(arrayAddr, structureId, bucketNum * BUCKET_SIZE);
    for (int i = 0; i < bucketNum; ++i) {
      bucketSort(structureId, numRow1[i], bucketAddr[i]);
    }
    kWayMergeSort(structureId, structureId + 1, numRow1, bucketAddr, bucketNum);
    resultId = structureId + 1;
  } else {
    for (int i = 0; i < bucketNum; ++i) {
      bucketSort(structureId + 1, numRow2[i], bucketAddr[i]);
    }
    kWayMergeSort(structureId + 1, structureId, numRow2, bucketAddr, bucketNum);
    resultId = structureId;
  }
  // test(arrayAddr, resultId, N);
  // print(arrayAddr, resultId, N);
  free(bucketAddr);
  free(numRow1);
  free(numRow2);
  return resultId;
}