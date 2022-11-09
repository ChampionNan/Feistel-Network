#ifndef COMMON_H
#define COMMON_H
#include <mbedtls/aes.h>
#include <mbedtls/des.h>
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"

#include <random>
#include <iostream>
#include <cstring>
#include <climits>
#include <cstdlib>
#include <cstdint>
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

#define N 167772160// 83886080// 2147483648L// 536870912// 1073741824
#define M 8388608// 4194304 // 16777216 // 2097152
#define BLOCK_DATA_SIZE 4
#define NUM_STRUCTURES 10
// #define MEM_IN_ENCLAVE 5
#define DUMMY -1
#define NULLCHAR '\0'
#define MY_RAND_MAX 2147483647

#define _ALPHA -1
#define _BETA -1
#define _P -1
#define ALPHA 0.0302070506959582
#define BETA 0.0286274342516871
#define P 22
#define KAPPA 27.8

extern int is_tight;
extern double IOcost;
extern int *arrayAddr[NUM_STRUCTURES];
extern int nonEnc;

typedef struct {
  int x;
  int key;
} Bucket_x;

typedef struct {
  __uint128_t x;
  __uint128_t iv;
} EncBlock;

// TODO: set up structure size
// BOS: 0, 1, 2; ODS: 3, 4
const int structureSize[NUM_STRUCTURES] = {sizeof(int), sizeof(Bucket_x), sizeof(Bucket_x), sizeof(int), sizeof(int)};

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

void fyShuffle(int structureId, int size, int B);
void opOneLinearScanBlock(int index, int* block, int blockSize, int structureId, int write, int dummyNum);
void ocall_print_string(const char *str);
void aes_init();
void aes2_init();
void cbc_encrypt(int* buffer, int blockSize);
void cbc_decrypt(int* buffer, int blockSize);
void OcallRB(int index, int* buffer, int blockSize, int structureId);
void OcallReadBlock(int startIdx, int offset, int* buffer, int blockSize, int structureId);
void OcallWB(int index, int offset, int* buffer, int blockSize, int structureId);
void OcallWriteBlock(int startIdx, int offset, int* buffer, int blockSize, int structureId);
void freeAllocate(int structureIdM, int structureIdF, int size);
int myrand();
int Hypergeometric(int NN, int Msize, int n_prime);
void shuffle(int *array, int n);
// void padWithDummy(int structureId, int start, int realNum, int secSize);
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
void printEnc(int **arrayAddr, int structureId, int size);
void print(int **arrayAddr, int structureId, int size);
void test(int **arrayAddr, int structureId, int size);
void testEnc(int **arrayAddr, int structureId, int size);
void testWithDummy(int **arrayAddr, int structureId, int size);

#endif // !COMMON_H
