#include "bucket.h"
#include "common.h"

// int finalFlag = 0;
int64_t oldK;
int64_t HEAP_NODE_SIZE;
int64_t MERGE_BATCH_SIZE;
extern int64_t FAN_OUT;
extern int64_t BUCKET_SIZE;
std::random_device dev3;
std::mt19937 rng3(dev3());
// BUCKET_SORT
void initMerge(int64_t size) {
  HEAP_NODE_SIZE = size;
  MERGE_BATCH_SIZE = size;
}

bool isTargetIterK(int64_t randomKey, int64_t iter, int64_t k, int64_t num) {
  /*
  while (iter) {
    randomKey = randomKey / oldK;
    iter--;
  }*/
  // return (randomKey & (0x01 << (iter - 1))) == 0 ? false : true;
  randomKey = (randomKey >> iter * oldK);
  return (randomKey % k) == num;
}

void mergeSplitHelper(Bucket_x *inputBuffer, int64_t* numRow1, int64_t* numRow2, int64_t* inputId, int64_t* outputId, int64_t iter, int64_t k, int64_t* bucketAddr, int outputStructureId) {
  // int batchSize = BUCKET_SIZE; // 8192
  // TODO: FREE these malloc
  // printf("In mergesplitHep!\n");
  Bucket_x **buf = (Bucket_x**)malloc(k * sizeof(Bucket_x*));
  for (int64_t i = 0; i < k; ++i) {
    buf[i] = (Bucket_x*)malloc(BUCKET_SIZE * sizeof(Bucket_x));
  }
  // int counter0 = 0, counter1 = 0;
  int64_t randomKey;
  int64_t *counter = (int64_t*)malloc(k * sizeof(int64_t));
  memset(counter, 0, k * sizeof(int64_t));
  // printf("Before assign bins!\n");
  for (int64_t i = 0; i < k * BUCKET_SIZE; ++i) {
    if (inputBuffer[i].x != DUMMY) {
      randomKey = inputBuffer[i].key;
      for (int64_t j = 0; j < k; ++j) {
        if (isTargetIterK(randomKey, iter, k, j)) {
          buf[j][counter[j] % BUCKET_SIZE] = inputBuffer[i];
          counter[j]++;
          // std::cout << "couter j: " << counter[j] << std::endl;
          if (counter[j] % BUCKET_SIZE == 0) {
            opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] +  numRow2[outputId[j]]), (int64_t*)buf[j], BUCKET_SIZE, outputStructureId, 1, 0);
            numRow2[outputId[j]] += BUCKET_SIZE;
          }
        }
      }
    }
  }
  // printf("Before final assign!\n");
  for (int64_t j = 0; j < k; ++j) {
    if (numRow2[outputId[j]] + (counter[j] % BUCKET_SIZE) > BUCKET_SIZE) {
      printf("overflow error during merge split!\n");
    }
    opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] + numRow2[outputId[j]]), (int64_t*)buf[j], counter[j] % BUCKET_SIZE, outputStructureId, 1, BUCKET_SIZE-numRow2[outputId[j]]-(counter[j] % BUCKET_SIZE));
    numRow2[outputId[j]] += counter[j] % BUCKET_SIZE;
    // padWithDummy(outputStructureId, bucketAddr[outputId[j]], numRow2[outputId[j]], BUCKET_SIZE);
    free(buf[j]);
  }
  free(counter);
  free(buf);
  // printf("Finish mergesplitHep!\n");
}

void mergeSplit(int inputStructureId, int outputStructureId, int64_t *inputId, int64_t *outputId, int64_t k, int64_t* bucketAddr, int64_t* numRow1, int64_t* numRow2, int64_t iter) {
  // step1. Read k buckets together
  // printf("In mergesplit!\n");
  Bucket_x *inputBuffer = (Bucket_x*)malloc(k * sizeof(Bucket_x) * BUCKET_SIZE);
  for (int64_t i = 0; i < k; ++i) {
    opOneLinearScanBlock(2 * bucketAddr[inputId[i]], (int64_t*)(&inputBuffer[i * BUCKET_SIZE]), BUCKET_SIZE, inputStructureId, 0, 0);
  }
  // step2. process k buckets
  mergeSplitHelper(inputBuffer, numRow1, numRow2, inputId, outputId, iter, k, bucketAddr, outputStructureId);
  free(inputBuffer);
  // printf("Finish mergesplit!\n");
}

void kWayMergeSort(int inputStructureId, int outputStructureId, int64_t* numRow1, int64_t* bucketAddr, int64_t bucketNum) {
  int64_t mergeSortBatchSize = HEAP_NODE_SIZE; // 256
  int64_t writeBufferSize = MERGE_BATCH_SIZE; // 8192
  int64_t numWays = bucketNum;
  printf("mergeSortBatchSize %d, %d\n", mergeSortBatchSize, writeBufferSize);
  // HeapNode inputHeapNodeArr[numWays];
  HeapNode *inputHeapNodeArr = (HeapNode*)malloc(numWays * sizeof(HeapNode));
  int64_t totalCounter = 0;
  
  int64_t *readBucketAddr = (int64_t*)malloc(sizeof(int64_t) * numWays);
  memcpy(readBucketAddr, bucketAddr, sizeof(int64_t) * numWays);
  int64_t writeBucketAddr = 0;
  int64_t j = 0;
  printf("Initialize\n");
  for (int64_t i = 0; i < numWays; ++i) {
    // TODO: 数据0跳过
    if (numRow1[i] == 0) {
      continue;
    }
    HeapNode node;
    node.data = (Bucket_x*)malloc(mergeSortBatchSize * sizeof(Bucket_x));
    node.bucketIdx = i;
    node.elemIdx = 0;
    opOneLinearScanBlock(2 * readBucketAddr[i], (int64_t*)node.data, std::min(mergeSortBatchSize, numRow1[i]), inputStructureId, 0, 0);
    inputHeapNodeArr[j++] = node;
    readBucketAddr[i] += std::min(mergeSortBatchSize, numRow1[i]);
  }
  
  Heap heap(inputHeapNodeArr, j, mergeSortBatchSize);
  Bucket_x *writeBuffer = (Bucket_x*)malloc(writeBufferSize * sizeof(Bucket_x));
  int64_t writeBufferCounter = 0;
  printf("Start merging\n");
  while (1) {
    HeapNode *temp = heap.getRoot();
    memcpy(writeBuffer + writeBufferCounter, temp->data + temp->elemIdx % mergeSortBatchSize, sizeof(Bucket_x));
    writeBufferCounter ++;
    totalCounter ++;
    temp->elemIdx ++;
    
    if (writeBufferCounter == writeBufferSize) {
      opOneLinearScanBlock(2 * writeBucketAddr, (int64_t*)writeBuffer, writeBufferSize, outputStructureId, 1, 0);
      writeBucketAddr += writeBufferSize;
      // numRow2[temp->bucketIdx] += writeBufferSize;
      writeBufferCounter = 0;
      // print(arrayAddr, outputStructureId, numWays * BUCKET_SIZE);
    }
    
    if (temp->elemIdx < numRow1[temp->bucketIdx] && (temp->elemIdx % mergeSortBatchSize) == 0) {
      opOneLinearScanBlock(2 * readBucketAddr[temp->bucketIdx], (int64_t*)(temp->data), std::min(mergeSortBatchSize, numRow1[temp->bucketIdx]-temp->elemIdx), inputStructureId, 0, 0);
      
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
  opOneLinearScanBlock(2 * writeBucketAddr, (int64_t*)writeBuffer, writeBufferCounter, outputStructureId, 1, 0);
  // numRow2[0] += writeBufferCounter;
  // TODO: ERROR writeBuffer
  free(writeBuffer);
  free(readBucketAddr);
  free(inputHeapNodeArr);
}

void bucketSort(int inputStructureId, int64_t size, int64_t dataStart) {
  Bucket_x *arr = (Bucket_x*)malloc(size * sizeof(Bucket_x));
  opOneLinearScanBlock(2 * dataStart, (int64_t*)arr, size, inputStructureId, 0, 0);
  quickSort(arr, 0, size - 1);
  opOneLinearScanBlock(2 * dataStart, (int64_t*)arr, size, inputStructureId, 1, 0);
  free(arr);
}

// int inputTrustMemory[BLOCK_DATA_SIZE];
int bucketOSort(int structureId, int64_t size) {
  int64_t k = FAN_OUT;
  int64_t bucketNum = smallestPowerOfKLargerThan(ceil(2.0 * size / BUCKET_SIZE), 2);
  int64_t mem1 = k * BUCKET_SIZE;
  HEAP_NODE_SIZE = fmax(floor(1.0 * M / (bucketNum+1)), 1);
  MERGE_BATCH_SIZE = HEAP_NODE_SIZE;
  std::cout << "Mem1's memory: " << mem1 << std::endl;
  std::cout << "Mem2's memory: " << HEAP_NODE_SIZE*(bucketNum+1) << std::endl;
  if (mem1 > M || MERGE_BATCH_SIZE > M) {
    std::cout << "Memory exceed\n";
  }
  // TODO: Change to ceil, not interger divide
  int64_t ranBinAssignIters = ceil(1.0*log2(bucketNum)/log2(k));
  std::cout << "Iteration times: " << ceil(1.0*log2(bucketNum)/log2(k)) << std::endl;
  srand((unsigned)time(NULL));
  // std::cout << "Iters:" << ranBinAssignIters << std::endl;
  int64_t *bucketAddr = (int64_t*)malloc(bucketNum * sizeof(int64_t));
  for (int64_t i = 0; i < bucketNum; ++i) {
    bucketAddr[i] = i * BUCKET_SIZE;
  }
  int64_t *numRow1 = (int64_t*)malloc(bucketNum * sizeof(int64_t));
  memset(numRow1, 0, bucketNum * sizeof(int64_t));
  int64_t *numRow2 = (int64_t*)malloc(bucketNum * sizeof(int64_t));
  memset(numRow2, 0, bucketNum * sizeof(int64_t));
  
  Bucket_x *trustedMemory = (Bucket_x*)malloc(BUCKET_SIZE * sizeof(Bucket_x));
  int64_t *inputTrustMemory = (int64_t*)malloc(BUCKET_SIZE * sizeof(int64_t));
  int64_t total = 0, readStart = 0;
  int64_t each;
  int64_t avg = size / bucketNum;
  int64_t remainder = size % bucketNum;
  std::uniform_int_distribution<int64_t> dist{0, bucketNum-1};

  for (int64_t i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    opOneLinearScanBlock(readStart, inputTrustMemory, each, structureId - 1, 0, 0);
    readStart += each;
    int64_t randomKey;
    for (int64_t j = 0; j < each; ++j) {
      // randomKey = (int)oe_rdrand();
      // randomKey = rand();
      randomKey = dist(rng3);
      trustedMemory[j].x = inputTrustMemory[j];
      trustedMemory[j].key = randomKey;
    }
    opOneLinearScanBlock(2 * bucketAddr[i], (int64_t*)trustedMemory, each, structureId, 1, BUCKET_SIZE-each);
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
  int64_t *inputId = (int64_t*)malloc(k * sizeof(int64_t));
  int64_t *outputId = (int64_t*)malloc(k *sizeof(int64_t));
  int64_t outIdx = 0, expo, tempk, jboundary, jjboundary;
  // TODO: change k to tempk, i != last level, k = FAN_OUT
  // last level k = k % k1, 2N/Z < 2^k, 2^k1 < M/Z
  int64_t k1 = (int64_t)log2(k);
  // oldK = (int)pow(2, k1);
  oldK = k1;
  int64_t tempk_i, tempk_i1;
  printf("Before RBA!\n");
  printf("Total bits: %lu, each level solved %lu bits\n", log2(bucketNum), k1);
  for (int64_t i = 0; i < ranBinAssignIters; ++i) {
    expo = std::min((int64_t)k1, (int64_t)(log2(bucketNum)-i*k1));
    tempk = pow(2, expo);
    tempk_i = pow(tempk, i);
    tempk_i1 = tempk_i * tempk;
    printf("level %lu read&write %lu buckets\n", i, tempk);
    jboundary = (i != (ranBinAssignIters-1)) ? (bucketNum / tempk_i1) : 1;
    // jboundary = bucketNum / tempk_i1; 
    jjboundary = (i != (ranBinAssignIters-1)) ? tempk_i : (bucketNum/tempk);
    printf("jboundary: %lu, jjboundary: %lu\n", jboundary, jjboundary);
    // finalFlag = (i != (ranBinAssignIters-1)) ? 0 : 1;
    if (i % 2 == 0) {
      for (int64_t j = 0; j < jboundary; ++j) {
        // pass (i-1) * k^i
        //printf("j: %d\n", j);
        for (int64_t jj = 0; jj < jjboundary; ++jj) {
          //printf("jj: %d\n", jj);
          for (int64_t m = 0; m < tempk; ++m) {
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
      int64_t count = 0;
      for (int64_t n = 0; n < bucketNum; ++n) {
        numRow1[n] = 0;
        count += numRow2[n];
      }
      printf("after %luth merge split, we have %lu tuples\n", i, count);
      // testRealNum(structureId + 1);
      // cnt.clear();
      outIdx = 0;
      //print(arrayAddr, structureId + 1, bucketNum * BUCKET_SIZE);
    } else {
      for (int64_t j = 0; j < jboundary; ++j) {
        //printf("j: %d\n", j);
        for (int64_t jj = 0; jj < jjboundary; ++jj) {
          //printf("jj: %d\n", jj);
          for (int64_t m = 0; m < tempk; ++m) {
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
      int64_t count = 0;
      for (int64_t n = 0; n < bucketNum; ++n) {
        numRow2[n] = 0;
        count += numRow1[n];
      }
      printf("after %luth merge split, we have %lu tuples\n", i, count);
      // testRealNum(structureId);
      // cnt.clear();
      outIdx = 0;
      //print(arrayAddr, structureId, bucketNum * BUCKET_SIZE);
    }
    std::cout << "----------------------------------------\n";
    printf("\n\n Finish random bin assignment iter%luth out of %lu\n\n", i, ranBinAssignIters);
    std::cout << "----------------------------------------\n";
  }
  free(inputId);
  free(outputId);
  int64_t resultId = 0;
  printf("Finish RBA!\n");
  if (ranBinAssignIters % 2 == 0) {
    print(arrayAddr, structureId, bucketNum * BUCKET_SIZE);
    for (int64_t i = 0; i < bucketNum; ++i) {
      bucketSort(structureId, numRow1[i], bucketAddr[i]);
    }
    kWayMergeSort(structureId, structureId + 1, numRow1, bucketAddr, bucketNum);
    resultId = structureId + 1;
  } else {
    for (int64_t i = 0; i < bucketNum; ++i) {
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