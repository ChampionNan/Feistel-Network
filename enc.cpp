#include <mbedtls/aes.h>
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"

#include <set>
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
#include <fstream>
#include <algorithm>
#include <iomanip>
// TODO: new
#include <unordered_set>
#include <unordered_map>
#include <bitset>

#define N 268435456// 1073741824L
#define M 16777216// 16777216
#define BLOCK_DATA_SIZE 4
#define NUM_STRUCTURES 10
// #define MEM_IN_ENCLAVE 5
#define DUMMY 0xffffffff
#define NULLCHAR '\0'
#define MY_RAND_MAX 2147483647

#define _ALPHA -1
#define _BETA -1
#define _P -1
#define ALPHA 0.0344927628034036
#define BETA 0.206333806761951
#define P 12

#define KAPPA 28
// #define MERGE_BATCH_SIZE 20 // merge split hepler

typedef struct {
  int x;
  int key;
} Bucket_x;
// OCALL
void ocall_print_string(const char *str);
void OcallReadBlock(int index, int* buffer, size_t blockSize, int structureId);
void OcallWriteBlock(int index, int* buffer, size_t blockSize, int structureId);
void freeAllocate(int structureIdM, int structureIdF, int size);
void opOneLinearScanBlock(int index, int* block, size_t blockSize, int structureId, int write, int dummyNum);
// OQSORT
void floydSampler(int n, int k, std::vector<int> &x);
int Sample(int inStructureId, int sampleSize, std::vector<int> &trustedM2, int is_tight, int is_rec);
void SampleRec(int inStructureId, int sampleId, int sortedSampleId, int is_tight, std::vector<std::vector<int> >& pivots);
void quantileCal(std::vector<int> &samples, int start, int end, int p);
int partitionMulti(int *arr, int low, int high, int pivot);
void quickSortMulti(int *arr, int low, int high, std::vector<int> pivots, int left, int right, std::vector<int> &partitionIdx);
std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec);
std::pair<int, int> TwoLevelPartition(int inStructureId, std::vector<std::vector<int> >& pivots, int p, int outStructureId1, int outStructureId2);
int ObliviousTightSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2);
int ObliviousTightSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId, int outStructureId2);
std::pair<int, int> ObliviouLooseSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2);
std::pair<int, int> ObliviousLooseSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2);
void ObliviousLooseSortRec(int sampleId, int sampleSize, int sortedSampleId, std::vector<std::vector<int> >& pivots);
// SUPPORT
void prf(unsigned char *right, char key, int tweak, unsigned char *ret);
void round(unsigned char* data, char key, int tweak, unsigned char* newData);
int encrypt(int index, char key[8], int rounds);
void callSort(int sortId, int structureId, int paddedSize, int *resId, int *resN);
int myrand();
int Hypergeometric(int NN, int Msize, int n_prime);
void shuffle(int *array, int n);
void padWithDummy(int structureId, int start, int realNum, int secSize);
int moveDummy(int *a, int size);
void swapRow(int *a, int *b);
bool cmpHelper(int *a, int *b);
bool cmpHelper(Bucket_x *a, Bucket_x *b);
int partition(int *arr, int low, int high);
void quickSort(int *arr, int low, int high);
int partition(Bucket_x *arr, int low, int high);
void quickSort(Bucket_x *arr, int low, int high);
void init(int **arrayAddr, int structurId, int size);
void print(int* array, int size);
void print(int **arrayAddr, int structureId, int size);
void test(int **arrayAddr, int structureId, int size);
void testWithDummy(int **arrayAddr, int structureId, int size);
int greatestPowerOfTwoLessThan(double n);
int smallestPowerOfKLargerThan(int n, int k);
bool isTargetIterK(int randomKey, int iter, int k, int num);
void mergeSplitHelper(Bucket_x *inputBuffer, int* numRow1, int* numRow2, int* inputId, int* outputId, int iter, int k, int* bucketAddr, int outputStructureId);
void mergeSplit(int inputStructureId, int outputStructureId, int *inputId, int *outputId, int k, int* bucketAddr, int* numRow1, int* numRow2, int iter);
void kWayMergeSort(int inputStructureId, int outputStructureId, int* numRow1, int* bucketAddr, int bucketNum);
void bucketSort(int inputStructureId, int size, int dataStart);
int bucketOSort(int structureId, int size);

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
int HEAP_NODE_SIZE;
int MERGE_BATCH_SIZE;

double IOcost = 0;
int is_tight = 1;
int finalFlag = 0;
int oldK;

struct HeapNode {
  Bucket_x *data;
  int bucketIdx;
  int elemIdx;
};


class Heap {
  HeapNode *harr;
  int heapSize;
  int batchSize;
public:
  Heap(HeapNode *a, int size, int bsize);
  void Heapify(int i);
  int left(int i);
  int right (int i);
  void swapHeapNode(HeapNode *a, HeapNode *b);
  HeapNode *getRoot();
  int getHeapSize();
  bool reduceSizeByOne();
  void replaceRoot(HeapNode x);
};



Heap::Heap(HeapNode *a, int size, int bsize) {
  heapSize = size;
  harr = a;
  int i = (heapSize - 1) / 2;
  batchSize = bsize;
  while (i >= 0) {
    Heapify(i);
    i --;
  }
}

void Heap::Heapify(int i) {
  int l = left(i);
  int r = right(i);
  int target = i;

  if (l < heapSize && cmpHelper(harr[i].data + harr[i].elemIdx % batchSize, harr[l].data + harr[l].elemIdx % batchSize)) {
    target = l;
  }
  if (r < heapSize && cmpHelper(harr[target].data + harr[target].elemIdx % batchSize, harr[r].data + harr[r].elemIdx % batchSize)) {
    target = r;
  }
  if (target != i) {
    swapHeapNode(&harr[i], &harr[target]);
    Heapify(target);
  }
}

int Heap::left(int i) {
  return (2 * i + 1);
}

int Heap::right(int i) {
  return (2 * i + 2);
}

void Heap::swapHeapNode(HeapNode *a, HeapNode *b) {
  HeapNode temp = *a;
  *a = *b;
  *b = temp;
}

HeapNode* Heap::getRoot() {
  return &harr[0];
}

int Heap::getHeapSize() {
  return heapSize;
}

bool Heap::reduceSizeByOne() {
  free(harr[0].data);
  heapSize --;
  if (heapSize > 0) {
    harr[0] = harr[heapSize];
    Heapify(0);
    return true;
  } else {
    return false;
  }
}

void Heap::replaceRoot(HeapNode x) {
  harr[0] = x;
  Heapify(0);
}

// TODO: set up structure size
const int structureSize[NUM_STRUCTURES] = {sizeof(int),
  2 * sizeof(int), 2 * sizeof(int),
  sizeof(int), sizeof(int), sizeof(int), sizeof(int)};

mbedtls_aes_context aes;
unsigned char key[16];
int base;
int max_num;
int ROUND = 3;

__uint128_t prf(__uint128_t a) {
  unsigned char input[16] = {0};
  unsigned char encrypt_output[16] = {0};
  for (int i = 0; i < 16; ++i) {
    input[i] = (a >> (120 - i * 8)) & 0xFF;
  }
  mbedtls_aes_crypt_ecb(&aes, MBEDTLS_AES_ENCRYPT, input, encrypt_output);
  __uint128_t res = 0;
  for (int i = 0; i < 16; ++i) {
    res |= encrypt_output[i] << (120 - i * 8);
  }
  return res;
}

int encrypt(int index) {
  int l = index / (1 << base);
  int r = index % (1 << base);
  __uint128_t e;
  int temp, i = 1;
  while (i <= ROUND) {
    e = prf((r << 16 * 8 - base) + i);
    temp = r;
    r = l ^ (e >> 16 * 8 - base);
    l = temp;
    i += 1;
  }
  return (l << base) + r;
}

void pseudo_init(int size) {
  mbedtls_aes_init(&aes);
  mbedtls_ctr_drbg_context ctr_drbg;
  mbedtls_entropy_context entropy;
  char *pers = "aes generate key";
  int ret;
  mbedtls_entropy_init(&entropy);
  mbedtls_ctr_drbg_init(&ctr_drbg);
  if((ret = mbedtls_ctr_drbg_seed(&ctr_drbg, mbedtls_entropy_func, &entropy, (unsigned char *) pers, strlen(pers))) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_init returned -0x%04x\n", -ret);
    return;
  }
  if((ret = mbedtls_ctr_drbg_random(&ctr_drbg, key, 16)) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_random returned -0x%04x\n", -ret);
    return ;
  }
  mbedtls_aes_setkey_enc(&aes, key, 128);
  base = ceil(1.0 * log2(size) / 2);
  max_num = 1 << 2 * base;
}

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

void print(int* array, int size) {
  int i;
  for (i = 0; i < size; i++) {
    printf("%d ", array[i]);
    if ((i != 0) && (i % 5 == 0)) {
      printf("\n");
    }
  }
  printf("\n");
}

void print(int **arrayAddr, int structureId, int size) {
  int i;
  std::ofstream fout("/home/data/bchenba/output.txt");
  if(structureSize[structureId] == 4) {
    int *addr = (int*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      // printf("%d ", addr[i]);
      fout << addr[i] << " ";
      if ((i != 0) && (i % 8 == 0)) {
        // printf("\n");
        fout << std::endl;
      }
    }
  } else if (structureSize[structureId] == 8) {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      // printf("(%d, %d) ", addr[i].x, addr[i].key);
      fout << "(" << addr[i].x << ", " << addr[i].key << ") ";
      if ((i != 0) && (i % 5 == 0)) {
        // printf("\n");
        fout << std::endl;
      }
    }
  }
  // printf("\n");
  fout << std::endl;
  fout.close();
}

// TODO: change nt types
void test(int **arrayAddr, int structureId, int size) {
  int pass = 1;
  int i;
  // print(structureId);
  if(structureSize[structureId] == 4) {
    for (i = 1; i < size; i++) {
      pass &= ((arrayAddr[structureId])[i-1] <= (arrayAddr[structureId])[i]);
      if (!pass) {
        std::cout << (arrayAddr[structureId])[i-1] << ' ' << (arrayAddr[structureId])[i];
        break;
      }
    }
  }
  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}

void testWithDummy(int **arrayAddr, int structureId, int size) {
  int i = 0;
  int j = 0;
  // print(structureId);
  if(structureSize[structureId] == 4) {
    for (i = 0; i < size; ++i) {
      if ((arrayAddr[structureId])[i] != DUMMY) {
        break;
      }
    }
    if (i == size - 1) {
      printf(" TEST PASSed\n");
      return;
    }
    while (i < size && j < size) {
      for (j = i + 1; j < size; ++j) {
        if ((arrayAddr[structureId])[j] != DUMMY) {
          break;
        }
      }
      if (j == size) { // Only 1 element not dummy
        printf(" TEST PASSed\n");
        return;
      }
      if ((arrayAddr[structureId])[i] < (arrayAddr[structureId])[j]) {
        i = j;
      } else {
        printf(" TEST FAILed\n");
        return;
      }
    }
    printf(" TEST PASSed\n");
    return;
  }
}

/* OCall functions */
void ocall_print_string(const char *str) {
  /* Proxy/Bridge will check the length and null-terminate
   * the input string to prevent buffer overflow.
   */
  printf("%s", str);
  fflush(stdout);
}

void OcallReadBlock(int index, int* buffer, size_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(buffer, arrayAddr[structureId] + index, blockSize * structureSize[structureId]);
  memcpy(buffer, arrayAddr[structureId] + index, blockSize);
  IOcost += ceil(1.0*blockSize/structureSize[structureId]/BLOCK_DATA_SIZE);
}

void OcallWriteBlock(int index, int* buffer, size_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(arrayAddr[structureId] + index, buffer, blockSize * structureSize[structureId]);
  memcpy(arrayAddr[structureId] + index, buffer, blockSize);
  IOcost += ceil(1.0*blockSize/structureSize[structureId]/BLOCK_DATA_SIZE);
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

int greatestPowerOfTwoLessThan(double n) {
    int k = 1;
    while (k > 0 && k < n) {
        k = k << 1;
    }
    return k >> 1;
}

int smallestPowerOfKLargerThan(int n, int k) {
  int num = 1;
  while (num > 0 && num < n) {
    num = num * k;
  }
  return num;
}

// TODO: Set this function as OCALL
void freeAllocate(int structureIdM, int structureIdF, int size) {
  // 1. Free arrayAddr[structureId]
  if (arrayAddr[structureIdF]) {
    free(arrayAddr[structureIdF]);
  }
  // 2. malloc new asked size (allocated in outside)
  if (size <= 0) {
    return;
  }
  int *addr = (int*)malloc(size * sizeof(int));
  memset(addr, DUMMY, size * sizeof(int));
  // 3. assign malloc address to arrayAddr
  arrayAddr[structureIdM] = addr;
  return ;
}

// Functions x crossing the enclave boundary, unit: BLOCK_DATA_SIZE
// TODO: FIX This function ? need checking
void opOneLinearScanBlock(int index, int* block, size_t blockSize, int structureId, int write, int dummyNum=0) {
  if (blockSize + dummyNum == 0) {
    return ;
  }
  int boundary = (int)((blockSize + BLOCK_DATA_SIZE - 1 )/ BLOCK_DATA_SIZE);
  int Msize, i;
  int multi = structureSize[structureId] / sizeof(int);
  if (!write) {
    // OcallReadBlock(index, block, blockSize * structureSize[structureId], structureId);
    OcallReadBlock(index, block, blockSize * structureSize[structureId], structureId);
  } else {
    // OcallWriteBlock(index, block, blockSize * structureSize[structureId], structureId);
    OcallWriteBlock(index, block, blockSize * structureSize[structureId], structureId);
    if (dummyNum) {
      int *junk = (int*)malloc(dummyNum * multi * sizeof(int));
      for (int j = 0; j < dummyNum * multi; ++j) {
        junk[j] = DUMMY;
      }
      int startIdx = index + multi * blockSize;
      OcallWriteBlock(startIdx , junk, dummyNum * structureSize[structureId], structureId);
    }
  }
  return;
}


void floydSampler(int n, int k, std::vector<int> &x) {
  std::unordered_set<int> H;
  for (int i = n - k; i < n; ++i) {
    x.push_back(i);
  }
  unsigned int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine e(seed1);
  int r, j, temp;
  for (int i = 0; i < k; ++i) {
    std::uniform_int_distribution<int> dist{0, n-k+1+i};
    r = dist(e); // get random numbers with PRNG
    if (H.count(r)) {
      std::uniform_int_distribution<int> dist2{0, i};
      j = dist2(e);
      temp = x[i];
      x[i] = x[j];
      x[j] = temp;
      H.insert(n-k+i);
    } else {
      x[i] = r;
      H.insert(r);
    }
  }
  sort(x.begin(), x.end());
}

// TODO: Unity Different Sample
int Sample(int inStructureId, int sampleSize, std::vector<int> &trustedM2, int is_tight, int is_rec=0) {
  std::cout << "In sample\n";
  int N_prime = sampleSize;
  double alpha = (!is_rec) ? ALPHA : _ALPHA;
  int n_prime = ceil(1.0 * alpha * N_prime);
  int boundary = ceil(1.0 * N_prime / BLOCK_DATA_SIZE);
  int j = 0, Msize;
  int *trustedM1 = (int*)malloc(BLOCK_DATA_SIZE * sizeof(int));
  std::vector<int> sampleIdx;
  floydSampler(N_prime, n_prime, sampleIdx);
  for (int i = 0; i < boundary; ++i) {
    if (is_tight) {
      Msize = std::min(BLOCK_DATA_SIZE, N_prime - i * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(i * BLOCK_DATA_SIZE, trustedM1, Msize, inStructureId, 0, 0);
      while ((j < n_prime) && (sampleIdx[j] >= i * BLOCK_DATA_SIZE) && (sampleIdx[j] < (i+1) * BLOCK_DATA_SIZE)) {
        trustedM2.push_back(trustedM1[sampleIdx[j] % BLOCK_DATA_SIZE]);
        j += 1;
      }
    } else if ((!is_tight) && (sampleIdx[j] >= i * BLOCK_DATA_SIZE) && (sampleIdx[j] < (i+1) * BLOCK_DATA_SIZE)) {
      Msize = std::min(BLOCK_DATA_SIZE, N_prime - i * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(i * BLOCK_DATA_SIZE, trustedM1, Msize, inStructureId, 0, 0);
      while ((sampleIdx[j] >= i * BLOCK_DATA_SIZE) && (sampleIdx[j] < (i+1) * BLOCK_DATA_SIZE)) {
        trustedM2.push_back(trustedM1[sampleIdx[j] % BLOCK_DATA_SIZE]);
        j += 1;
        if (j >= n_prime) break;
      }
      if (j >= n_prime) break;
    }
  }
  sort(trustedM2.begin(), trustedM2.end());
  return n_prime;
}

// TODO: What's return value?
void SampleRec(int inStructureId, int sampleId, int sortedSampleId, int is_tight, std::vector<std::vector<int> >& pivots) {
  int N_prime = N;
  int n_prime = ceil(1.0 * ALPHA * N_prime);
  int boundary = ceil(1.0 * N / M);
  int realNum = 0;
  int readStart = 0;
  int *trustedM1 = (int*)malloc(M * sizeof(int));
  int m = 0, Msize;
  freeAllocate(sampleId, sampleId, n_prime);
  for (int i = 0; i < boundary; ++i) {
    Msize = std::min(M, N - i * M);
    // TODO: USing boost library
    m = Hypergeometric(N_prime, Msize, n_prime);
    if (is_tight || (!is_tight && m > 0)) {
      opOneLinearScanBlock(readStart, trustedM1, Msize, inStructureId, 0, 0);
      readStart += Msize;
      shuffle(trustedM1, Msize);
      opOneLinearScanBlock(realNum, trustedM1, m, sampleId, 1, 0);
      realNum += m;
      n_prime -= m;
    }
    N_prime -= Msize;
  }
  std::cout << "Till Sample IOcost: " << 1.0*IOcost/N*BLOCK_DATA_SIZE << std::endl;
  if (realNum > M) {
    ObliviousLooseSortRec(sampleId, realNum, sortedSampleId, pivots);
  }
  return ;
}

// TODO: vector quantile
void quantileCal(std::vector<int> &samples, int start, int end, int p) {
  int sampleSize = end - start;
  for (int i = 1; i < p; ++i) {
    samples[i] = samples[i * sampleSize / p];
  }
  samples[0] = INT_MIN;
  samples[p] = INT_MAX;
  samples.resize(p+1);
  samples.shrink_to_fit();
  return ;
}

// TODO: Add new multi-pivots quicksort
int partitionMulti(int *arr, int low, int high, int pivot) {
  int i = low - 1;
  for (int j = low; j < high + 1; ++j) {
    if (pivot > arr[j]) {
      i += 1;
      swapRow(arr + i, arr + j);
    }
  }
  return i;
}

void quickSortMulti(int *arr, int low, int high, std::vector<int> pivots, int left, int right, std::vector<int> &partitionIdx) {
  int pivotIdx, pivot, mid;
  if (right >= left) {
    pivotIdx = (left + right) >> 1;
    pivot = pivots[pivotIdx];
    mid = partitionMulti(arr, low, high, pivot);
    partitionIdx.push_back(mid);
    quickSortMulti(arr, low, mid, pivots, left, pivotIdx-1, partitionIdx);
    quickSortMulti(arr, mid+1, high, pivots, pivotIdx+1, right, partitionIdx);
  }
}

int BSFirstGreaterEqual(std::vector<int> &nums,int target)
{
    int left = 0, right = nums.size() - 1, res = right;
    while(left <= right)
    {
        int mid = left + (right - left) / 2 ;
        if(nums[mid] >= target)
        {
            res = mid; right = mid - 1;
        }else
        {
            left = mid + 1;
        }
 
    }
    if(target<=nums[res]) return res;
    return -1;
}

// Partition num to =target & <target
// Return the first element index greater than target
int partitionEqual(int *num, int size, int target) {
  int i = -1;
  for (int j = 0; j < size; ++j) {
    if (num[j] == target) {
      i++;
      if (i != j) {
        swapRow(num+i, num+j);
      }
    }
  }
  return i;
}

/*
std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec=0, int is_duplicate=0) {
  if (inSize <= M) {
    return {inSize, 1};
  }
  double beta = (!is_rec) ? BETA : _BETA;
  int hatN = ceil(1.0 * (1 + 2 * beta) * inSize);
  int M_prime = ceil(1.0 * M / (1 + 2 * beta));
  int r = ceil(1.0 * log(hatN / M) / log(p));
  int p0 = ceil(1.0 * hatN / (M * pow(p, r - 1)));
  quantileCal(samples, 0, sampleSize, p0);
  for (int i = 0; i < samples.size(); ++i) {
    std::cout << samples[i] << ' ';
  }
  mbedtls_aes_setkey_enc(&aes, key, 256);
  std::set<int> hash;
  int i = 0, res;
  while (hash.size() < 2000000) {
    res = encrypt(i, key, 3);
    if (hash.count(res)) {
      std::cout << "Duplicate!\n";
      return 0;
    }
    if (res >= 0 && res < 2000000 && !hash.count(res)) {
      hash.insert(res);
    }
  }
  // TODO: count each memory block # duplicate keys, find out why misssing some pivot numbers
  // TODO: Find FFSEM implementation in c++
  for (int i = 0; i < boundary1; ++i) {
    for (int j = 0; j < boundary2; ++j) {
      if (total_blocks - 1 - blocks_done == 0) {
        k = 0;
      } else {
        k = rand() % (total_blocks - blocks_done);
      }
      Msize1 = std::min(BLOCK_DATA_SIZE, inSize - k * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(k * BLOCK_DATA_SIZE, &trustedM3[j*BLOCK_DATA_SIZE], Msize1, inStructureId, 0, 0);
      memset(shuffleB, DUMMY, sizeof(int) * BLOCK_DATA_SIZE);
      Msize2 = std::min(BLOCK_DATA_SIZE, inSize - (total_blocks-1-blocks_done) * BLOCK_DATA_SIZE);
      opOneLinearScanBlock((total_blocks-1-blocks_done) * BLOCK_DATA_SIZE, shuffleB, Msize2, inStructureId, 0, 0);
      opOneLinearScanBlock(k * BLOCK_DATA_SIZE, shuffleB, BLOCK_DATA_SIZE, inStructureId, 1, 0);
      blocks_done += 1;
      if (blocks_done == total_blocks) {
        break;
      }
    }
    int blockNum = moveDummy(trustedM3, dataBoundary);
    quickSortMulti(trustedM3, 0, blockNum-1, samples, 1, p0, partitionIdx);
    sort(partitionIdx.begin(), partitionIdx.end());
    partitionIdx.insert(partitionIdx.begin(), -1);
    for (int j = 0; j < p0; ++j) {
      index1 = partitionIdx[j]+1;
      index2 = partitionIdx[j+1];
      // std::cout << trustedM3[index1] << ' ' << trustedM3[index2] << std::endl;
      writeBackNum = index2 - index1 + 1;
      if (writeBackNum == 0) {
        continue;
      }
      // Find out #elements=pivots=samples[j]
      equalNum = 0;
      eachNum = 0;
      // TODO: partition trustedM3[index1]
      index3 = partitionEqual(&trustedM3[index1], writeBackNum, samples[j]);
      if (j != 0 && index3 != -1) {
        // std::cout << trustedM3[index1+index3] << std::endl;
        
        index3 += index1 + 1;
        equalNum = index3 - index1;
        int parts = dupPivots[samples[j]];
        int average = equalNum / parts;
        int remainder = equalNum % parts;
        int startJ = BSFirstGreaterEqual(samples, samples[j]);
        int readM3Idx = 0;
        for (int m = 0; m < parts; ++m) {
          eachNum = average + ((m < remainder) ? 1 : 0);
          if (eachNum > smallSectionSize) {
            std::cout << "Overflow in small section1 M/p0: " << eachNum << std::endl;
          }
          opOneLinearScanBlock((startJ+m)*bucketSize0+i*smallSectionSize, &trustedM3[index1+readM3Idx], eachNum, outStructureId1, 1, smallSectionSize-eachNum);
          readM3Idx += eachNum;
        }
      }
      // TODO: Change reading start
      if (eachNum+writeBackNum-equalNum > smallSectionSize) {
        std::cout << "Overflow in small section2 M/p0: " << eachNum+writeBackNum-equalNum << std::endl;
      }
      opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize + eachNum, &trustedM3[index1+equalNum], writeBackNum-equalNum, outStructureId1, 1, smallSectionSize-eachNum-(writeBackNum-equalNum));
    }
    memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
    partitionIdx.clear();
  }
  free(trustedM3);
  free(shuffleB);
  if (bucketSize0 > M) {
    std::cout << "Each section size is greater than M, adjst parameters: " << bucketSize0 << ", " << M << std::endl;
  }
  return {bucketSize0, p0};
}*/

std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec) {
  if (inSize <= M) {
    return {inSize, 1};
  }
  std::cout << "In Onelevel Partition\n";
  double beta = (!is_rec) ? BETA : _BETA;
  int hatN = ceil(1.0 * (1 + 2 * beta) * inSize);
  int M_prime = ceil(1.0 * M / (1 + 2 * beta));
  int r = ceil(1.0 * log(hatN / M) / log(p));
  int p0 = ceil(1.0 * hatN / (M * pow(p, r - 1)));
  quantileCal(samples, 0, sampleSize, p0);
  int boundary1 = ceil(1.0 * inSize / M_prime);
  int boundary2 = ceil(1.0 * M_prime / BLOCK_DATA_SIZE);
  int dataBoundary = boundary2 * BLOCK_DATA_SIZE;
  int smallSectionSize = M / p0;
  int bucketSize0 = boundary1 * smallSectionSize;
  freeAllocate(outStructureId1, outStructureId1, boundary1 * smallSectionSize * p0);
  
  int Msize1, Msize2, index1, index2, writeBackNum;
  int total_blocks = ceil(1.0 * inSize / BLOCK_DATA_SIZE);
  int *trustedM3 = (int*)malloc(sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  int *shuffleB = (int*)malloc(sizeof(int) * BLOCK_DATA_SIZE);
  std::vector<int> partitionIdx;
  // Finish FFSEM implementation in c++
  pseudo_init(total_blocks);
  int index_range = max_num;
  int k = 0, read_index;
  for (int i = 0; i < boundary1; ++i) {
    for (int j = 0; j < boundary2; ++j) {
      read_index = encrypt(k);
      while (read_index >= total_blocks) {
        k += 1;
        if (k == index_range) {
          k = -1;
          break;
        }
        read_index = encrypt(k);
      }
      if (k == -1) {
        break;
      }
      Msize1 = std::min(BLOCK_DATA_SIZE, inSize - read_index * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(read_index * BLOCK_DATA_SIZE, &trustedM3[j*BLOCK_DATA_SIZE], Msize1, inStructureId, 0, 0);
      k += 1;
      if (k == index_range) {
        break;
      }
    }
    int blockNum = moveDummy(trustedM3, dataBoundary);
    quickSortMulti(trustedM3, 0, blockNum-1, samples, 1, p0, partitionIdx);
    sort(partitionIdx.begin(), partitionIdx.end());
    partitionIdx.insert(partitionIdx.begin(), -1);
    for (int j = 0; j < p0; ++j) {
      index1 = partitionIdx[j]+1;
      index2 = partitionIdx[j+1];
      writeBackNum = index2 - index1 + 1;
      if (writeBackNum > smallSectionSize) {
        printf("Overflow in small section M/p0: %d", writeBackNum);
      }
      opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
    }
    memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
    partitionIdx.clear();
  }
  free(trustedM3);
  free(shuffleB);
  if (bucketSize0 > M) {
    printf("Each section size is greater than M, adjst parameters: %d, %d", bucketSize0, M);
  }
  return {bucketSize0, p0};
}

// TODO: Add TwoLevelPartition
std::pair<int, int> TwoLevelPartition(int inStructureId, std::vector<std::vector<int> >& pivots, int p, int outStructureId1, int outStructureId2) {
  int M_prime = ceil(1.0 * M / (1 + 2 * BETA));
  int p0 = p;
  int boundary1 = ceil(1.0 * N / M_prime);
  int boundary2 = ceil(1.0 * M_prime / BLOCK_DATA_SIZE);
  int dataBoundary = boundary2 * BLOCK_DATA_SIZE;
  int smallSectionSize = M / p0;
  int bucketSize0 = boundary1 * smallSectionSize;
  freeAllocate(outStructureId1, outStructureId1, boundary1 * smallSectionSize * p0);
  int k, Msize1, Msize2, index1, index2, writeBackNum;
  int blocks_done = 0;
  int total_blocks = ceil(1.0 * N / BLOCK_DATA_SIZE);
  int *trustedM3 = (int*)malloc(sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  int *shuffleB = (int*)malloc(sizeof(int) * BLOCK_DATA_SIZE);
  std::vector<int> partitionIdx;
  for (int i = 0; i < boundary1; ++i) {
    for (int j = 0; j < boundary2; ++j) {
      if (total_blocks - 1 - blocks_done == 0) {
        k = 0;
      } else {
        k = rand() % (total_blocks - blocks_done);
      }
      Msize1 = std::min(BLOCK_DATA_SIZE, N - k * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(k * BLOCK_DATA_SIZE, &trustedM3[j*BLOCK_DATA_SIZE], Msize1, inStructureId, 0, 0);
      memset(shuffleB, DUMMY, sizeof(int) * BLOCK_DATA_SIZE);
      Msize2 = std::min(BLOCK_DATA_SIZE, N - (total_blocks-1-blocks_done) * BLOCK_DATA_SIZE);
      opOneLinearScanBlock((total_blocks-1-blocks_done) * BLOCK_DATA_SIZE, shuffleB, Msize2, inStructureId, 0, 0);
      opOneLinearScanBlock(k * BLOCK_DATA_SIZE, shuffleB, BLOCK_DATA_SIZE, inStructureId, 1, 0);
      blocks_done += 1;
      if (blocks_done == total_blocks) {
        break;
      }
    }
    int blockNum = moveDummy(trustedM3, dataBoundary);
    quickSortMulti(trustedM3, 0, blockNum-1, pivots[0], 1, p0, partitionIdx);
    sort(partitionIdx.begin(), partitionIdx.end());
    partitionIdx.insert(partitionIdx.begin(), -1);
    for (int j = 0; j < p0; ++j) {
      index1 = partitionIdx[j]+1;
      index2 = partitionIdx[j+1];
      writeBackNum = index2 - index1 + 1;
      if (writeBackNum > smallSectionSize) {
        std::cout << "Overflow in small section M/p0: " << writeBackNum << std::endl;
      }
      opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
    }
    memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
    partitionIdx.clear();
  }
  free(trustedM3);
  free(shuffleB);
  // Level2
  int p1 = p0 * p, readSize, readSize2, k1, k2;
  int boundary3 = ceil(1.0 * bucketSize0 / M);
  int bucketSize1 = boundary3 * smallSectionSize;
  freeAllocate(outStructureId2, outStructureId2, boundary3 * smallSectionSize * p0 * p);
  // std::vector<int> trustedM1;
  int *trustedM2 = (int*)malloc(sizeof(int) * M);
  // TODO: Change memory excceeds? use &trustedM2[M-1-p]
  int *trustedM2_part = (int*)malloc(sizeof(int) * (p+1));
  for (int j = 0; j < p0; ++j) {
    // trustedM1 = pivots[1+j];
    for (int k = 0; k < boundary3; ++k) {
      Msize1 = std::min(M, bucketSize0 - k * M);
      readSize = (Msize1 < (p+1)) ? Msize1 : (Msize1-(p+1));
      opOneLinearScanBlock(j*bucketSize0+k*M, trustedM2, readSize, outStructureId1, 0, 0);
      k1 = moveDummy(trustedM2, readSize);
      readSize2 = (Msize1 < (p+1)) ? 0 : (p+1);
      opOneLinearScanBlock(j*bucketSize0+k*M+readSize, trustedM2_part, readSize2, outStructureId1, 0, 0);
      k2 = moveDummy(trustedM2_part, readSize2);
      memcpy(&trustedM2[k1], trustedM2_part, sizeof(int) * k2);
      quickSortMulti(trustedM2, 0, k1+k2-1, pivots[1+j], 1, p, partitionIdx);
      sort(partitionIdx.begin(), partitionIdx.end());
      partitionIdx.insert(partitionIdx.begin(), -1);
      for (int ll = 0; ll < p; ++ll) {
        index1 = partitionIdx[ll]+1;
        index2 = partitionIdx[ll+1];
        writeBackNum = index2 - index1 + 1;
        if (writeBackNum > smallSectionSize) {
          std::cout << "Overflow in small section M/p0: " << writeBackNum << std::endl;
        }
        opOneLinearScanBlock((j*p0+ll)*bucketSize1+k*smallSectionSize, &trustedM2[index1], writeBackNum, outStructureId2, 1, smallSectionSize-writeBackNum);
      }
      memset(trustedM2, DUMMY, sizeof(int) * M);
      partitionIdx.clear();
    }
  }
  if (bucketSize1 > M) {
    std::cout << "Each section size is greater than M, adjust parameters: " << bucketSize1 << std::endl;
  }
  return {bucketSize1, p1};
}


int ObliviousTightSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2) {
  int *trustedM;
  std::cout << "In ObliviousTightSort\n";
  if (inSize <= M) {
    trustedM = (int*)malloc(sizeof(int) * M);
    opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
    quickSort(trustedM, 0, inSize - 1);
    freeAllocate(outStructureId1, outStructureId1, inSize);
    opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
    free(trustedM);
    return outStructureId1;
  }
  std::vector<int> trustedM2;
  int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
  std::pair<int, int> section= OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  int sectionSize = section.first;
  int sectionNum = section.second;
  // TODO: IN order to reduce memory, can replace outStructureId2 with inStructureId
  freeAllocate(outStructureId2, outStructureId2, inSize);
  trustedM = (int*)malloc(sizeof(int) * M);
  int j = 0, k;
  std::cout << "In final\n";
  for (int i = 0; i < sectionNum; ++i) {
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k-1);
    opOneLinearScanBlock(j, trustedM, k, outStructureId2, 1, 0);
    j += k;
    if (j > inSize) {
      std::cout << "Final error" << std::endl;
    }
  }
  free(trustedM);
  return outStructureId2;
}

// TODO: TightSort2
int ObliviousTightSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2) {
  std::cout << "In ObliviousTightSort2 && In SampleRec\n";
  std::vector<std::vector<int> > pivots;
  SampleRec(inStructureId, sampleId, sortedSampleId, 1, pivots);
  std::cout << "In TwoLevelPartition\n";
  std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
  std::cout << "Till Partition IOcost: " << IOcost/N*BLOCK_DATA_SIZE << std::endl;
  int sectionSize = section.first;
  int sectionNum = section.second;
  int *trustedM = (int*)malloc(sizeof(int) * M);
  int j = 0, k;
  for (int i = 0; i < sectionNum; ++i) {
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k-1);
    opOneLinearScanBlock(j, trustedM, k, outStructureId1, 1, 0);
    j += k;
    if (j > inSize) {
      std::cout << "Final error2\n";
    }
  }
  std::cout << "Till Final IOcost: " << IOcost/N*BLOCK_DATA_SIZE << std::endl;
  return outStructureId1;
}

std::pair<int, int> ObliviousLooseSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2) {
  std::cout << "In ObliviousLooseSort\n";
  int *trustedM;
  if (inSize <= M) {
    trustedM = (int*)malloc(sizeof(int) * M);
    opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
    quickSort(trustedM, 0, inSize - 1);
    opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
    free(trustedM);
    return {outStructureId1, inSize};
  }
  std::vector<int> trustedM2;
  int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
  std::pair<int, int> section = OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  int sectionSize = section.first;
  int sectionNum = section.second;
  int totalLevelSize = sectionNum * sectionSize;
  int k;
  freeAllocate(outStructureId2, outStructureId2, totalLevelSize);
  trustedM = (int*)malloc(sizeof(int) * M);
  for (int i = 0; i < sectionNum; ++i) {
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k - 1);
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 1, 0);
  }
  return {outStructureId2, totalLevelSize};
}

// TODO: Finish
std::pair<int, int> ObliviousLooseSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2) {
  std::cout << "In ObliviousLooseSort2 && In SasmpleRec\n";
  std::vector<std::vector<int> > pivots;
  SampleRec(inStructureId, sampleId, sortedSampleId, 0, pivots);
  std::cout << "Till Sample IOcost: " << IOcost/N*BLOCK_DATA_SIZE << std::endl;
  std::cout << "In TwoLevelPartition\n";
  std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
  int sectionSize = section.first;
  int sectionNum = section.second;
  int totalLevelSize = sectionNum * sectionSize;
  int k;
  freeAllocate(outStructureId1, outStructureId1, totalLevelSize);
  int *trustedM = (int*)malloc(sizeof(int) * M);
  for (int i = 0; i < sectionNum; ++i) {
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k - 1);
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 1, 0);
  }
  std::cout << "Till Final IOcost: " << IOcost/N*BLOCK_DATA_SIZE << std::endl;
  return {outStructureId1, totalLevelSize};
}

// TODO: Finish, what's return value
void ObliviousLooseSortRec(int sampleId, int sampleSize, int sortedSampleId, std::vector<std::vector<int> >& pivots) {
  std::cout << "In ObliviousLooseSortRec\n";
  std::vector<int> trustedM2;
  int realNum = Sample(sampleId, sampleSize, trustedM2, 0, 1);
  std::cout << "In OneLevelPartition\n";
  std::pair<int, int> section = OneLevelPartition(sampleId, sampleSize, trustedM2, realNum, _P, sortedSampleId, 1);
  int sectionSize = section.first;
  int sectionNum = section.second;
  int j = 0, k = 0, total = 0;
  int outj = 0, inj = 0;
  int *trustedM = (int*)malloc(sizeof(int) * M);
  std::vector<int> quantileIdx;
  for (int i = 1; i < P; ++i) {
    quantileIdx.push_back(i * sampleSize / P);
  }
  int size = ceil(1.0 * sampleSize / P);
  std::vector<std::vector<int> > quantileIdx2;
  std::vector<int> index;
  for (int i = 0; i < P; ++i) {
    for (int j = 1; j < P; ++j) {
      index.push_back(i * size + j * size / P);
    }
    quantileIdx2.push_back(index);
    index.clear();
  }
  std::vector<int> pivots1, pivots2_part;
  for (int i = 0; i < sectionNum; ++i) {
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, sortedSampleId, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k-1);
    total += k;
    // Cal Level1 pivots
    while ((j < P-1) && (quantileIdx[j] < total)) {
      pivots1.push_back(trustedM[quantileIdx[j]-(total-k)]);
      j += 1;
    }
    // Cal Level2 pivots
    while (outj < P) {
      while ((inj < P-1) && (quantileIdx2[outj][inj] < total)) {
        pivots2_part.push_back(trustedM[quantileIdx2[outj][inj]-(total-k)]);
        inj += 1;
        if (inj == P-1) {
          inj = 0;
          outj += 1;
          pivots.push_back(pivots2_part);
          pivots2_part.clear();
          break;
        }
      }
      if (outj == P || quantileIdx2[outj][inj] >= total) {
        break;
      }
    }
    if ((j >= P-1) && (outj >= P)) {
      break;
    }
  }
  pivots1.insert(pivots1.begin(), INT_MIN);
  pivots1.push_back(INT_MAX);
  for (int i = 0; i < P; ++i) {
    pivots[i].insert(pivots[i].begin(), INT_MIN);
    pivots[i].push_back(INT_MAX);
  }
  pivots.insert(pivots.begin(), pivots1);
}

// BUCKET_SORT
void mergeSplitHelper(Bucket_x *inputBuffer, int* numRow1, int* numRow2, int* inputId, int* outputId, int iter, int k, int* bucketAddr, int outputStructureId) {
  // int batchSize = BUCKET_SIZE; // 8192
  // TODO: FREE these malloc
  Bucket_x **buf = (Bucket_x**)malloc(k * sizeof(Bucket_x*));
  for (int i = 0; i < k; ++i) {
    buf[i] = (Bucket_x*)malloc(BUCKET_SIZE * sizeof(Bucket_x));
  }
  
  // int counter0 = 0, counter1 = 0;
  int randomKey;
  int *counter = (int*)malloc(k * sizeof(int));
  memset(counter, 0, k * sizeof(int));
  
  for (int i = 0; i < k * BUCKET_SIZE; ++i) {
    if ((inputBuffer[i].key != DUMMY) && (inputBuffer[i].x != DUMMY)) {
      randomKey = inputBuffer[i].key;
      for (int j = 0; j < k; ++j) {
        if (isTargetIterK(randomKey, iter, k, j)) {
          buf[j][counter[j] % BUCKET_SIZE] = inputBuffer[i];
          counter[j]++;
          // std::cout << "couter j: " << counter[j] << std::endl;
          if (counter[j] % BUCKET_SIZE == 0) {
            opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] +  numRow2[outputId[j]]), (int*)buf[j], (size_t)BUCKET_SIZE, outputStructureId, 1);
            numRow2[outputId[j]] += BUCKET_SIZE;
          }
        }
      }
    }
  }
  
  for (int j = 0; j < k; ++j) {
    opOneLinearScanBlock(2 * (bucketAddr[outputId[j]] + numRow2[outputId[j]]), (int*)buf[j], (size_t)(counter[j] % BUCKET_SIZE), outputStructureId, 1);
    numRow2[outputId[j]] += counter[j] % BUCKET_SIZE;
    padWithDummy(outputStructureId, bucketAddr[outputId[j]], numRow2[outputId[j]], BUCKET_SIZE);
    if (numRow2[outputId[j]] > BUCKET_SIZE) {
      printf("overflow error during merge split!\n");
    }
    free(buf[j]);
  }
  free(counter);
}

void mergeSplit(int inputStructureId, int outputStructureId, int *inputId, int *outputId, int k, int* bucketAddr, int* numRow1, int* numRow2, int iter) {
  // step1. Read k buckets together
  Bucket_x *inputBuffer = (Bucket_x*)malloc(k * sizeof(Bucket_x) * BUCKET_SIZE);
  for (int i = 0; i < k; ++i) {
    opOneLinearScanBlock(2 * bucketAddr[inputId[i]], (int*)(&inputBuffer[i * BUCKET_SIZE]), BUCKET_SIZE, inputStructureId, 0);
  }
  // step2. process k buckets
  mergeSplitHelper(inputBuffer, numRow1, numRow2, inputId, outputId, iter, k, bucketAddr, outputStructureId);
  free(inputBuffer);
  
  
}

void kWayMergeSort(int inputStructureId, int outputStructureId, int* numRow1, int* bucketAddr, int bucketNum) {
  int mergeSortBatchSize = HEAP_NODE_SIZE; // 256
  int writeBufferSize = MERGE_BATCH_SIZE; // 8192
  int numWays = bucketNum;
  // HeapNode inputHeapNodeArr[numWays];
  HeapNode *inputHeapNodeArr = (HeapNode*)malloc(numWays * sizeof(HeapNode));
  int totalCounter = 0;
  
  int *readBucketAddr = (int*)malloc(sizeof(int) * numWays);
  memcpy(readBucketAddr, bucketAddr, sizeof(int) * numWays);
  int writeBucketAddr = 0;
  int j = 0;
  
  for (int i = 0; i < numWays; ++i) {
    // TODO: 数据0跳过
    if (numRow1[i] == 0) {
      continue;
    }
    HeapNode node;
    node.data = (Bucket_x*)malloc(mergeSortBatchSize * sizeof(Bucket_x));
    node.bucketIdx = i;
    node.elemIdx = 0;
    opOneLinearScanBlock(2 * readBucketAddr[i], (int*)node.data, (size_t)std::min(mergeSortBatchSize, numRow1[i]), inputStructureId, 0);
    inputHeapNodeArr[j++] = node;
    readBucketAddr[i] += std::min(mergeSortBatchSize, numRow1[i]);
  }
  
  Heap heap(inputHeapNodeArr, j, mergeSortBatchSize);
  Bucket_x *writeBuffer = (Bucket_x*)malloc(writeBufferSize * sizeof(Bucket_x));
  int writeBufferCounter = 0;

  while (1) {
    HeapNode *temp = heap.getRoot();
    memcpy(writeBuffer + writeBufferCounter, temp->data + temp->elemIdx % mergeSortBatchSize, sizeof(Bucket_x));
    writeBufferCounter ++;
    totalCounter ++;
    temp->elemIdx ++;
    
    if (writeBufferCounter == writeBufferSize) {
      opOneLinearScanBlock(2 * writeBucketAddr, (int*)writeBuffer, (size_t)writeBufferSize, outputStructureId, 1);
      writeBucketAddr += writeBufferSize;
      // numRow2[temp->bucketIdx] += writeBufferSize;
      writeBufferCounter = 0;
      // print(arrayAddr, outputStructureId, numWays * BUCKET_SIZE);
    }
    
    if (temp->elemIdx < numRow1[temp->bucketIdx] && (temp->elemIdx % mergeSortBatchSize) == 0) {
      opOneLinearScanBlock(2 * readBucketAddr[temp->bucketIdx], (int*)(temp->data), (size_t)std::min(mergeSortBatchSize, numRow1[temp->bucketIdx]-temp->elemIdx), inputStructureId, 0);
      
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
  opOneLinearScanBlock(2 * writeBucketAddr, (int*)writeBuffer, (size_t)writeBufferCounter, outputStructureId, 1);
  // numRow2[0] += writeBufferCounter;
  // TODO: ERROR writeBuffer
  free(writeBuffer);
  free(readBucketAddr);
  free(inputHeapNodeArr);
}

void bucketSort(int inputStructureId, int size, int dataStart) {
  Bucket_x *arr = (Bucket_x*)malloc(size * sizeof(Bucket_x));
  opOneLinearScanBlock(2 * dataStart, (int*)arr, (size_t)size, inputStructureId, 0);
  quickSort(arr, 0, size - 1);
  opOneLinearScanBlock(2 * dataStart, (int*)arr, (size_t)size, inputStructureId, 1);
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
  int ranBinAssignIters = ceil(1.0*log(bucketNum)/log(k));
  std::cout << "Iteration times: " << ceil(1.0*log(bucketNum)/log(k)) << std::endl;
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

  for (int i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    opOneLinearScanBlock(readStart, inputTrustMemory, each, structureId - 1, 0);
    readStart += each;
    int randomKey;
    for (int j = 0; j < each; ++j) {
      // randomKey = (int)oe_rdrand();
      randomKey = rand();
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
  oldK = (int)pow(2, k1);
  for (int i = 0; i < ranBinAssignIters; ++i) {
    expo = std::min(k1, (int)log2(bucketNum)-i*k1);
    tempk = (int)pow(2, expo);
    jboundary = (i != (ranBinAssignIters-1)) ? (bucketNum / (int)pow(tempk, i+1)) : 1;
    jjboundary = (i != (ranBinAssignIters-1)) ? ((int)pow(tempk, i)) : (bucketNum/tempk);
    finalFlag = (i != (ranBinAssignIters-1)) ? 0 : 1;
    if (i % 2 == 0) {
      for (int j = 0; j < jboundary; ++j) {
        // pass (i-1) * k^i
        //printf("j: %d\n", j);
        for (int jj = 0; jj < jjboundary; ++jj) {
          //printf("jj: %d\n", jj);
          for (int m = 0; m < tempk; ++m) {
            //printf("j, jj, m: %d, %d, %d\n", j, jj, m);
            inputId[m] = j * (int)pow(tempk, i+1)+ jj + m * jjboundary;
            outputId[m] = (outIdx * tempk + m) % bucketNum;
            // printf("input, output: %d, %d\n", inputId[m], outputId[m]);
          }/*
          for (int m = 0; m < tempk; ++m) {
            printf("input%d: %d\n", m, inputId[m]);
          }
          for (int m = 0; m < tempk; ++m) {
            printf("output%d: %d\n", m, outputId[m]);
          }*/
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
            inputId[m] = j * (int)pow(tempk, i+1)+ jj + m * jjboundary;
            outputId[m] = (outIdx * tempk + m) % bucketNum;
            //printf("input, output: %d, %d\n", inputId[m], outputId[m]);
          }/*
          for (int m = 0; m < tempk; ++m) {
            printf("input%d: %d\n", m, inputId[m]);
          }
          for (int m = 0; m < tempk; ++m) {
            printf("output%d: %d\n", m, outputId[m]);
          }*/
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

int merge_sort(int inStructureId, int outStructureId, int size) {
  int bucketNum = ceil(1.0 * size / M);
  std::cout << "Bucket Number: " << bucketNum << std::endl;
  HEAP_NODE_SIZE = floor(1.0 * M / (bucketNum + 1));
  MERGE_BATCH_SIZE = HEAP_NODE_SIZE;
  std::cout << "HEAP_NODE_SIZE: " << HEAP_NODE_SIZE << std::endl;
  int avg = size / bucketNum;
  int remainder = size % bucketNum;
  int *bucketAddr = (int*)malloc(bucketNum * sizeof(int));
  int *elements = (int*)malloc(bucketNum * sizeof(int));
  int each, total = 0;
  for (int i = 0; i < bucketNum; ++i) {
    each = avg + ((i < remainder) ? 1 : 0);
    elements[i] = each;
    bucketAddr[i] = total;
    total += each;
  }
  Bucket_x *arr = (Bucket_x*)malloc(M * sizeof(Bucket_x));
  for (int i = 0; i < bucketNum; ++i) {
    bucketSort(inStructureId, elements[i], bucketAddr[i]);
  }
  kWayMergeSort(inStructureId, outStructureId, elements, bucketAddr, bucketNum);
  return outStructureId;
}

// SUPPORTINGs
int myrand() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> dist(0, INT_MAX);
  return dist(mt);
}

int Hypergeometric(int NN, int Msize, int n_prime) {
  int m = 0;
  srand((unsigned)time(0));
  double rate = double(n_prime) / NN;
  for (int j = 0; j < Msize; ++j) {
    if (rand() / double(INT_MAX) < rate) {
      m += 1;
      n_prime -= 1;
    }
    NN -= 1;
    rate = double(n_prime) / double(NN);
  }
  return m;
}

void shuffle(int *array, int n) {
  if (n > 1) {
    for (int i = 0; i < n - 1; ++i) {
      int j = i + rand() / (INT_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

void padWithDummy(int structureId, int start, int realNum, int secSize) {
  int len = secSize - realNum;
  if (len <= 0) {
    return ;
  }
  
  if (structureSize[structureId] == 4) {
    int *junk = (int*)malloc(len * sizeof(int));
    for (int i = 0; i < len; ++i) {
      junk[i] = DUMMY;
    }
    opOneLinearScanBlock(start + realNum, (int*)junk, len, structureId, 1);
    free(junk);
  
  } else if (structureSize[structureId] == 8) {
    Bucket_x *junk = (Bucket_x*)malloc(len * sizeof(Bucket_x));
    for (int i = 0; i < len; ++i) {
      junk[i].x = DUMMY;
      junk[i].key = DUMMY;
    }
    opOneLinearScanBlock(2 * (start + realNum), (int*)junk, len, structureId, 1);
    free(junk);
  }
}

int moveDummy(int *a, int size) {
  // k: #elem != DUMMY
  int k = 0;
  for (int i = 0; i < size; ++i) {
    if (a[i] != DUMMY) {
      if (i != k) {
        swapRow(&a[i], &a[k++]);
      } else {
        k++;
      }
    }
  }
  return k;
}

bool isTargetIterK(int randomKey, int iter, int k, int num) {
  while (iter) {
    randomKey = randomKey / oldK;
    iter--;
  }
  // return (randomKey & (0x01 << (iter - 1))) == 0 ? false : true;
  return (randomKey % k) == num;
}

void swapRow(int *a, int *b) {
  int *temp = (int*)malloc(sizeof(int));
  memmove(temp, a, sizeof(int));
  memmove(a, b, sizeof(int));
  memmove(b, temp, sizeof(int));
  free(temp);
}

void swapRow(Bucket_x *a, Bucket_x *b) {
  Bucket_x *temp = (Bucket_x*)malloc(sizeof(Bucket_x));
  memmove(temp, a, sizeof(Bucket_x));
  memmove(a, b, sizeof(Bucket_x));
  memmove(b, temp, sizeof(Bucket_x));
  free(temp);
}

bool cmpHelper(int *a, int *b) {
  return (*a > *b) ? true : false;
}

bool cmpHelper(Bucket_x *a, Bucket_x *b) {
  return (a->x > b->x) ? true : false;
}

// TODO: Later change to vector, using internal sort
int partition(int *arr, int low, int high) {
  // TODO: random version
  // srand(unsigned(time(NULL)));
  int randNum = rand() % (high - low + 1) + low;
  swapRow(arr + high, arr + randNum);
  int *pivot = arr + high;
  int i = low - 1;
  for (int j = low; j <= high - 1; ++j) {
    if (cmpHelper(pivot, arr + j)) {
      i++;
      if (i != j) {
        swapRow(arr + i, arr + j);
      }
    }
  }
  if (i + 1 != high) {
    swapRow(arr + i + 1, arr + high);
  }
  return (i + 1);
}

int partition(Bucket_x *arr, int low, int high) {
  int randNum = rand() % (high - low + 1) + low;
  swapRow(arr + high, arr + randNum);
  Bucket_x *pivot = arr + high;
  int i = low - 1;
  for (int j = low; j <= high - 1; ++j) {
    if (cmpHelper(pivot, arr + j)) {
      i++;
      if (i != j) {
        swapRow(arr + i, arr + j);
      }
    }
  }
  if (i + 1 != high) {
    swapRow(arr + i + 1, arr + high);
  }
  return (i + 1);
}


void quickSort(int *arr, int low, int high) {
  if (high > low) {
    int mid = partition(arr, low, high);
    quickSort(arr, low, mid - 1);
    quickSort(arr, mid + 1, high);
  }
}

void quickSort(Bucket_x *arr, int low, int high) {
  if (high > low) {
    int mid = partition(arr, low, high);
    quickSort(arr, low, mid - 1);
    quickSort(arr, mid + 1, high);
  }
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
    if (paddedSize / M < 100) {
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

