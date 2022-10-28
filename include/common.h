#ifndef COMMON_H
#define COMMON_H
#include <mbedtls/aes.h>
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"

#include <random>
#include <iostream>
#include <cstring>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <string.h>
#include <vector>
#include <chrono>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <bitset>

#define N 134217728// 2147483648L// 536870912// 1073741824
#define M 4194304 // 16777216 // 2097152
#define BLOCK_DATA_SIZE 4
#define NUM_STRUCTURES 10
// #define MEM_IN_ENCLAVE 5
#define DUMMY -1
#define NULLCHAR '\0'
#define MY_RAND_MAX 2147483647

#define _ALPHA -1
#define _BETA -1
#define _P -1
#define ALPHA 0.0161416823830321
#define BETA 0.0575360376415013
#define P 36
#define KAPPA 27.8

extern int is_tight;
extern double IOcost;
extern int64_t *arrayAddr[NUM_STRUCTURES];

// TODO: set up structure size
const int structureSize[NUM_STRUCTURES] = {sizeof(int64_t),
  2 * sizeof(int64_t), 2 * sizeof(int64_t),
  sizeof(int64_t), sizeof(int64_t), sizeof(int64_t), sizeof(int64_t)};

typedef struct {
  int64_t x;
  int64_t key;
} Bucket_x;

struct HeapNode {
  Bucket_x *data;
  int64_t bucketIdx;
  int64_t elemIdx;
};

class Heap {
  HeapNode *harr;
  int64_t heapSize;
  int64_t batchSize;
public:
  Heap(HeapNode *a, int64_t size, int64_t bsize);
  void Heapify(int64_t i);
  int64_t left(int64_t i);
  int64_t right (int64_t i);
  void swapHeapNode(HeapNode *a, HeapNode *b);
  HeapNode *getRoot();
  int64_t getHeapSize();
  bool reduceSizeByOne();
  void replaceRoot(HeapNode x);
};

void opOneLinearScanBlock(int64_t index, int64_t* block, int64_t blockSize, int structureId, int write, int64_t dummyNum);
void ocall_print_string(const char *str);
void OcallReadBlock(int64_t index, int64_t* buffer, int64_t blockSize, int structureId);
void OcallWriteBlock(int64_t index, int64_t* buffer, int64_t blockSize, int structureId);
void freeAllocate(int structureIdM, int structureIdF, int64_t size);
int myrand();
int64_t Hypergeometric(int64_t NN, int64_t Msize, int64_t n_prime);
void shuffle(int64_t *array, int64_t n);
// void padWithDummy(int structureId, int start, int realNum, int secSize);
int64_t moveDummy(int64_t *a, int64_t size);
void swapRow(int64_t *a, int64_t *b);
bool cmpHelper(int64_t *a, int64_t *b);
bool cmpHelper(Bucket_x *a, Bucket_x *b);
int64_t partition(int64_t *arr, int64_t low, int64_t high);
void quickSort(int64_t *arr, int64_t low, int64_t high);
int64_t partition(Bucket_x *arr, int64_t low, int64_t high);
void quickSort(Bucket_x *arr, int64_t low, int64_t high);
int64_t greatestPowerOfTwoLessThan(double n);
int64_t smallestPowerOfKLargerThan(int64_t n, int64_t k);
void print(int64_t* array, int64_t size);
void print(int64_t **arrayAddr, int structureId, int64_t size);
void test(int64_t **arrayAddr, int structureId, int64_t size);
void testWithDummy(int64_t **arrayAddr, int structureId, int64_t size);

#endif // !COMMON_H
