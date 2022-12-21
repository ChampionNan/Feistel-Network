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

#define SMALLM 19954
#define PAGEBYTES 10
#define BLOCK_DATA_SIZE 255
#define NUM_STRUCTURES 10
// #define MEM_IN_ENCLAVE 5
#define DUMMY -1
#define NULLCHAR '\0'
#define MY_RAND_MAX 2147483647

#define _ALPHA -1
#define _BETA -1
#define _P -1
#define ALPHA 0.01
#define BETA 0.04
#define GAMMA 0.08
#define P 565
#define KAPPA 27.8
#define N 4194304000L // 83886080// 2147483648L// 536870912
const int M = 8388608; // 4194304 // 16777216 // 2097152

extern int64_t is_tight;
extern double IOcost;
extern int64_t *arrayAddr[NUM_STRUCTURES];
extern int64_t nonEnc;

typedef struct {
  int64_t x;
  int64_t key;
} Bucket_x;

typedef struct {
  uint32_t x[BLOCK_DATA_SIZE];
  __uint128_t iv;
} EncBlock;

// TODO: set up structure size
// BOS: 0, 1, 2; ODS: 3, 4
const int64_t structureSize[NUM_STRUCTURES] = {sizeof(int64_t), sizeof(Bucket_x), sizeof(Bucket_x), sizeof(int64_t), sizeof(int64_t)};

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

void fyShuffle(int64_t structureId, int64_t size, int64_t B);
void opOneLinearScanBlock(int64_t index, int64_t* block, int64_t blockSize, int64_t structureId, int64_t write, int64_t dummyNum);
void ocall_print_string(const char *str);
void aes_init();
void aes2_init();
void cbc_encrypt(int64_t* buffer, int64_t blockSize);
void cbc_decrypt(int64_t* buffer, int64_t blockSize);
void OcallRB(int64_t index, int64_t* buffer, int64_t blockSize, int64_t structureId);
void OcallReadBlock(int64_t startIdx, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId);
void OcallWB(int64_t index, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId);
void OcallWriteBlock(int64_t startIdx, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId);
void freeAllocate(int64_t structureIdM, int64_t structureIdF, int64_t size);
int64_t myrand();
int64_t Hypergeometric(int64_t NN, int64_t Msize, int64_t n_prime);
void shuffle(int64_t *array, int64_t n);
// void padWithDummy(int64_t structureId, int64_t start, int64_t realNum, int64_t secSize);
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
void printEnc(int64_t **arrayAddr, int64_t structureId, int64_t size);
void print(int64_t **arrayAddr, int64_t structureId, int64_t size);
void test(int64_t **arrayAddr, int64_t structureId, int64_t size);
void testEnc(int64_t **arrayAddr, int64_t structureId, int64_t size);
void testWithDummy(int64_t **arrayAddr, int64_t structureId, int64_t size);

#endif // !COMMON_H
