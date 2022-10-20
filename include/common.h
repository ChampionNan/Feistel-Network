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

extern int is_tight;
extern double IOcost;
extern int *arrayAddr[NUM_STRUCTURES];

// TODO: set up structure size
const int structureSize[NUM_STRUCTURES] = {sizeof(int),
  2 * sizeof(int), 2 * sizeof(int),
  sizeof(int), sizeof(int), sizeof(int), sizeof(int)};

typedef struct {
  int x;
  int key;
} Bucket_x;


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

void opOneLinearScanBlock(int index, int* block, size_t blockSize, int structureId, int write, int dummyNum);
void ocall_print_string(const char *str);
void OcallReadBlock(int index, int* buffer, size_t blockSize, int structureId);
void OcallWriteBlock(int index, int* buffer, size_t blockSize, int structureId);
void freeAllocate(int structureIdM, int structureIdF, int size);
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
int greatestPowerOfTwoLessThan(double n);
int smallestPowerOfKLargerThan(int n, int k);
void print(int* array, int size);
void print(int **arrayAddr, int structureId, int size);
void test(int **arrayAddr, int structureId, int size);
void testWithDummy(int **arrayAddr, int structureId, int size);

#endif // !COMMON_H
