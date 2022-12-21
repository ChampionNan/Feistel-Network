#include "common.h"

std::random_device dev2;
std::mt19937 rng2(dev2());
unsigned char key2[16];
mbedtls_aes_context aes2;
mbedtls_ctr_drbg_context ctr_drbg;
mbedtls_entropy_context entropy;
size_t iv_offset, iv_offset1;

Heap::Heap(HeapNode *a, int64_t size, int64_t bsize) {
  heapSize = size;
  harr = a;
  int64_t i = (heapSize - 1) / 2;
  batchSize = bsize;
  while (i >= 0) {
    Heapify(i);
    i --;
  }
}

void Heap::Heapify(int64_t i) {
  int64_t l = left(i);
  int64_t r = right(i);
  int64_t target = i;

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

int64_t Heap::left(int64_t i) {
  return (2 * i + 1);
}

int64_t Heap::right(int64_t i) {
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

int64_t Heap::getHeapSize() {
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

// support int64_t version only
void fyShuffle(int64_t structureId, int64_t size, int64_t B) {
  if (size % B != 0) {
    printf("Error! Not BLOCK times.\n");
  }
  int64_t total_blocks = ceil(size / B);
  EncBlock *trustedM3 = (EncBlock*)malloc(sizeof(EncBlock));
  int64_t k;
  std::random_device dev;
  std::mt19937 rng(dev()); 
  // switch block i & block k
  for (int64_t i = total_blocks-1; i >= 0; i--) {
    std::uniform_int_distribution<int64_t> dist(0, i);
    k = dist(rng);
    memcpy(trustedM3, (EncBlock*)(arrayAddr[structureId]) + k, sizeof(EncBlock));
    memcpy((EncBlock*)(arrayAddr[structureId]) + k, (EncBlock*)(arrayAddr[structureId]) + i, sizeof(EncBlock));
    memcpy((EncBlock*)(arrayAddr[structureId]) + i, trustedM3, sizeof(EncBlock));
  }
  // std::cout << "Finished floyd shuffle\n";
}

/* OCall functions */
void ocall_print_string(const char *str) {
  /* Proxy/Bridge will check the length and null-terminate
   * the input string to prevent buffer overflow.
   */
  printf("%s", str);
  fflush(stdout);
}

void aes_init() {
  mbedtls_aes_init(&aes2);
  char *pers = "aes2 generate key";
  int64_t ret;
  mbedtls_entropy_init(&entropy);
  mbedtls_ctr_drbg_init(&ctr_drbg);
  if((ret = mbedtls_ctr_drbg_seed(&ctr_drbg, mbedtls_entropy_func, &entropy, (unsigned char *) pers, strlen(pers))) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_init returned -0x%04x\n", -ret);
    return;
  }
  if((ret = mbedtls_ctr_drbg_random(&ctr_drbg, key2, 16)) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_random returned -0x%04x\n", -ret);
    return ;
  }
  mbedtls_aes_setkey_enc(&aes2, key2, 128);
  mbedtls_aes_setkey_dec(&aes2, key2, 128);
}

// Assume blockSize = 16 * k
void cbc_encrypt(EncBlock* buffer, int64_t blockSize) {
  // std::cout<< "In cbc_encrypt\n";
  mbedtls_ctr_drbg_random(&ctr_drbg, (uint8_t*)(&(buffer->iv)), 16);
  // // // std::cout<< "memcpy iv\n";
  unsigned char iv[16];
  iv_offset = 0;
  memcpy(iv, (uint8_t*)(&(buffer->iv)), 16);
  mbedtls_aes_crypt_ofb(&aes2, blockSize, &iv_offset, (uint8_t*)(&(buffer->iv)), (uint8_t*)buffer, (uint8_t*)buffer);
  memcpy((uint8_t*)(&(buffer->iv)), iv, 16);
  return;
}

// Assume blockSize = 16 * k
void cbc_decrypt(EncBlock* buffer, int64_t blockSize) {
  // std::cout<< "In cbc_decrypt\n";
  // mbedtls_aes_crypt_cfb8(&aes2, MBEDTLS_AES_DECRYPT, blockSize, (uint8_t*)(&(buffer->iv)), (uint8_t*)buffer, (uint8_t*)buffer);
  iv_offset1 = 0;
  mbedtls_aes_crypt_ofb(&aes2, blockSize, &iv_offset1, (uint8_t*)(&(buffer->iv)), (uint8_t*)buffer, (uint8_t*)buffer);
  return;
}

// index: EncBlock, offset: int64_t, blockSize: bytes
void OcallRB(int64_t index, int64_t* buffer, int64_t blockSize, int64_t structureId) {
  // std::cout<< "In OcallRB\n";
  memcpy(buffer, &(((EncBlock*)(arrayAddr[structureId]))[index]), blockSize); 
  IOcost += 1;
}

// startIdx: index of blocks, offset: data offset in the block, blockSize: bytes of real data
void OcallReadBlock(int64_t startIdx, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId) {
  if (blockSize == 0) {
    return;
  }
  // std::cout<< "In OcallReadBlock\n";
  EncBlock *readBuffer = (EncBlock*)malloc(sizeof(EncBlock));
  if (nonEnc) {
    OcallRB(startIdx, (int64_t*)readBuffer, sizeof(EncBlock), structureId);
    int64_t *a = (int64_t*)readBuffer;
    memcpy(buffer, (int64_t*)readBuffer+offset*(structureSize[structureId]/sizeof(int64_t)), blockSize);
  } else {
    OcallRB(startIdx, (int64_t*)readBuffer, sizeof(EncBlock), structureId);
    cbc_decrypt(readBuffer, sizeof(int64_t)*BLOCK_DATA_SIZE);
    memcpy(buffer, (int64_t*)readBuffer+offset*(structureSize[structureId]/sizeof(int64_t)), blockSize);
  }
  free(readBuffer);
}

// index: EncBlock, offset: int64_t, blockSize: bytes(only non-enc mode could use)
void OcallWB(int64_t index, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId) {
  // std::cout<< "In OcallWB\n";
  memcpy((int64_t*)(&(((EncBlock*)(arrayAddr[structureId]))[index]))+offset, buffer, blockSize);
  IOcost += 1;
}

// startIdx: index of blocks, offset(int64_t): data offset in the block, blockSize: bytes of real data
void OcallWriteBlock(int64_t startIdx, int64_t offset, int64_t* buffer, int64_t blockSize, int64_t structureId) {
  // std::cout<< "In OcallWriteBlock\n";
  if (blockSize == 0) {
    return;
  }
  EncBlock* writeBuf = (EncBlock*)malloc(sizeof(EncBlock));
  if (nonEnc) {
    OcallWB(startIdx, offset*(structureSize[structureId]/sizeof(int64_t)), buffer, blockSize, structureId);
  } else {
    if (offset == 0) { // could write the whole block
      // printf("In OcallWriteBlock1\n");
      memcpy((int64_t*)writeBuf, buffer, blockSize);
      cbc_encrypt(writeBuf, sizeof(int64_t)*BLOCK_DATA_SIZE);
      OcallWB(startIdx, 0, (int64_t*)writeBuf, sizeof(EncBlock), structureId);
    } else { // read&decrypt first, later write&encrypt
      // printf("In OcallWriteBlock2\n");
      OcallRB(startIdx, (int64_t*)writeBuf, sizeof(EncBlock), structureId);
      cbc_decrypt(writeBuf, sizeof(int64_t)*BLOCK_DATA_SIZE);
      memcpy((int64_t*)writeBuf+offset*(structureSize[structureId]/sizeof(int64_t)), buffer, blockSize);
      cbc_encrypt(writeBuf, sizeof(int64_t)*BLOCK_DATA_SIZE);
      OcallWB(startIdx, 0, (int64_t*)writeBuf, sizeof(EncBlock), structureId);
    }
  }
  free(writeBuf);
}

// TODO: Set this function as OCALL
// allocate encB for outside memory
void freeAllocate(int64_t structureIdM, int64_t structureIdF, int64_t size) {
  // 1. Free arrayAddr[structureId]
  if (arrayAddr[structureIdF]) {
    free(arrayAddr[structureIdF]);
  }
  // 2. malloc new asked size (allocated in outside)
  if (size <= 0) {
    return;
  }
  int64_t *addr = (int64_t*)malloc(size * sizeof(EncBlock));
  memset(addr, DUMMY, size * sizeof(EncBlock));
  // 3. assign malloc address to arrayAddr
  arrayAddr[structureIdM] = addr;
  return ;
}

// index: start index counted by elements (count from 0), blockSize: #elements
void opOneLinearScanBlock(int64_t index, int64_t* block, int64_t blockSize, int64_t structureId, int64_t write, int64_t dummyNum=0) {
  // std::cout<< "In opOneLinearScanBlock\n";
  if (blockSize + dummyNum == 0) {
    return ;
  }
  if (dummyNum < 0) {
    printf("Dummy padding error!");
    return ;
  }
  int64_t multi = structureSize[structureId] / sizeof(int64_t);
  int64_t B = BLOCK_DATA_SIZE / multi;
  int64_t startBIdx = index / B;
  int64_t startOffset = index % B;
  int64_t firstSize = B - startOffset;
  int64_t boundary = ceil(1.0 * (blockSize - firstSize) / B) + 1;
  int64_t opStart = 0;
  int64_t Msize, offset = 0;
  if (!write) {
    // std::cout<< "In reading B: \n";
    for (int64_t i = 0; i < boundary; ++i) {
      if (i != 0 && (i != boundary - 1)) {
        Msize = B;
      } else if (i == 0) {
        Msize = firstSize;
      } else {
        Msize = blockSize - opStart;
      }
      offset = (i == 0) ? startOffset : 0;
      // printf("Start Idx: %d, offset: %d, opstart: %d, Msize: %d\n", startBIdx + i, offset, opStart, Msize);
      OcallReadBlock(startBIdx + i, offset, &block[opStart * multi], Msize * structureSize[structureId], structureId);
      opStart += Msize;
    }
  } else {
    // printf("In opwriteReal\n");
    for (int64_t i = 0; i < boundary; ++i) {
      if (i != 0 && (i != boundary - 1)) {
        Msize = B;
      } else if (i == 0) {
        Msize = firstSize;
      } else {
        Msize = blockSize - opStart;
      }
      offset = (i == 0) ? startOffset : 0;
      OcallWriteBlock(startBIdx + i, offset, &block[opStart * multi], Msize * structureSize[structureId], structureId);
      opStart += Msize;
    }
    if (dummyNum > 0) {
      // printf("In opwriteDummy\n");
      int64_t *junk = (int64_t*)malloc(dummyNum * multi * sizeof(int64_t));
      for (int64_t j = 0; j < dummyNum * multi; ++j) {
        junk[j] = DUMMY;
      }
      int64_t dummyStart = (index + blockSize) / B;
      int64_t dummyOffset = (index + blockSize) % B;
      int64_t dummyFirstSize = B - dummyOffset;
      int64_t dummyBoundary = ceil(1.0 * (dummyNum - dummyFirstSize) / B) + 1;
      int64_t dummyOpStart = 0;
      for (int64_t j = 0; j < dummyBoundary; ++j) {
        if (j != 0 && (j != dummyBoundary - 1)) {
          Msize = B;
        } else if (j == 0) {
          Msize = dummyFirstSize;
        } else {
          Msize = dummyNum - dummyOpStart;
        }
        offset = (j == 0) ? dummyOffset : 0;
        OcallWriteBlock(dummyStart + j, offset, &junk[dummyOpStart * multi], Msize * structureSize[structureId], structureId);
        dummyOpStart += Msize; 
      }
    }
  }
}

// SUPPORTINGs
int64_t myrand() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int64_t> dist(0, INT_MAX);
  return dist(mt);
}

int64_t Hypergeometric(int64_t NN, int64_t Msize, int64_t n_prime) {
  int64_t m = 0;
  srand((unsigned)time(0));
  double rate = double(n_prime) / NN;
  for (int64_t j = 0; j < Msize; ++j) {
    if (rand() / double(INT_MAX) < rate) {
      m += 1;
      n_prime -= 1;
    }
    NN -= 1;
    rate = double(n_prime) / double(NN);
  }
  return m;
}

void shuffle(int64_t *array, int64_t n) {
  int64_t j, t;
  if (n > 1) {
    for (int64_t i = 0; i < n - 1; ++i) {
      // TODO: shuffle error
      j = i + rand() / (INT_MAX / (n - i) + 1);
      t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
  }
}

int64_t moveDummy(int64_t *a, int64_t size) {
  // k: #elem != DUMMY
  int64_t k = 0;
  for (int64_t i = 0; i < size; ++i) {
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

void swapRow(int64_t *a, int64_t *b) {
  int64_t *temp = (int64_t*)malloc(sizeof(int64_t));
  memmove(temp, a, sizeof(int64_t));
  memmove(a, b, sizeof(int64_t));
  memmove(b, temp, sizeof(int64_t));
  free(temp);
}

void swapRow(Bucket_x *a, Bucket_x *b) {
  Bucket_x *temp = (Bucket_x*)malloc(sizeof(Bucket_x));
  memmove(temp, a, sizeof(Bucket_x));
  memmove(a, b, sizeof(Bucket_x));
  memmove(b, temp, sizeof(Bucket_x));
  free(temp);
}

bool cmpHelper(int64_t *a, int64_t *b) {
  return (*a > *b) ? true : false;
}

bool cmpHelper(Bucket_x *a, Bucket_x *b) {
  return (a->x > b->x) ? true : false;
}

// TODO: Later change to vector, using internal sort
int64_t partition(int64_t *arr, int64_t low, int64_t high) {
  // TODO: random version
  // srand(unsigned(time(NULL)));
  std::uniform_int_distribution<int64_t> dist{low, high};
  int64_t randNum = dist(rng2);
  swapRow(arr + high, arr + randNum);
  int64_t *pivot = arr + high;
  int64_t i = low - 1;
  for (int64_t j = low; j <= high - 1; ++j) {
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

int64_t partition(Bucket_x *arr, int64_t low, int64_t high) {
  std::uniform_int_distribution<int64_t> dist{low, high};
  int64_t randNum = dist(rng2);
  // int64_t randNum = rand() % (high - low + 1) + low;
  swapRow(arr + high, arr + randNum);
  Bucket_x *pivot = arr + high;
  int64_t i = low - 1;
  for (int64_t j = low; j <= high - 1; ++j) {
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


void quickSort(int64_t *arr, int64_t low, int64_t high) {
  if (high > low) {
    int64_t mid = partition(arr, low, high);
    quickSort(arr, low, mid - 1);
    quickSort(arr, mid + 1, high);
  }
}

void quickSort(Bucket_x *arr, int64_t low, int64_t high) {
  if (high > low) {
    // printf("In quickSort: %ld, %ld\n", low, high);
    int64_t mid = partition(arr, low, high);
    quickSort(arr, low, mid - 1);
    quickSort(arr, mid + 1, high);
  }
}

int64_t greatestPowerOfTwoLessThan(double n) {
    int64_t k = 1;
    while (k > 0 && k < n) {
        k = k << 1;
    }
    return k >> 1;
}

int64_t smallestPowerOfKLargerThan(int64_t n, int64_t k) {
  int64_t num = 1;
  while (num > 0 && num < n) {
    num = num * k;
  }
  return num;
}

// SUPPORT
void print(int64_t* array, int64_t size) {
  int64_t i;
  for (i = 0; i < size; i++) {
    printf("%d ", array[i]);
    if ((i != 0) && (i % 6 == 0)) {
      printf("\n");
    }
  }
  printf("\n");
}

void print(int64_t **arrayAddr, int64_t structureId, int64_t size) {
  int64_t i;
  std::ofstream fout("/home/data/bchenba/output.txt");
  if(structureSize[structureId] == 8) {
    int64_t *addr = (int64_t*)arrayAddr[structureId];
    for (i = 0; i < size; i++) {
      // printf("%d ", addr[i]);
      fout << addr[i] << " ";
      if ((i != 0) && (i % 8 == 0)) {
        // printf("\n");
        fout << std::endl;
      }
    }
  } else if (structureSize[structureId] == 16) {
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

// size: real numbers
void printEnc(int64_t **arrayAddr, int64_t structureId, int64_t size) {
  int64_t i, j, blockNumber;
  std::ofstream fout("/home/data/bchenba/res.txt");
  printf("Print size: %d\n", size);
  EncBlock *addr = (EncBlock*)arrayAddr[structureId];
  if(structureSize[structureId] == 4) {
    // std::cout << "In printEnc\n";
    blockNumber = (int64_t)ceil(1.0*size/BLOCK_DATA_SIZE);
    for (i = 0; i < blockNumber; i++) {
      int64_t *addx = (int64_t*)(&(addr[i]));
      // fout << "i: " << i << std::endl;
      // // std::cout << "i: " << i << std::endl;
      for (j = 0; j < BLOCK_DATA_SIZE; ++j) {
        fout << addx[j] << ' ';
        std::cout << addx[j] << ' ';
      }
      if (i % 5 == 0) {
        fout << std::endl;
        std::cout << std::endl;
      }
    }
  } else if (structureSize[structureId] == 8) {
    blockNumber = (int64_t)ceil(1.0*size/(BLOCK_DATA_SIZE/2));
    for (i = 0; i < blockNumber; i++) {
      Bucket_x *addx = (Bucket_x*)(&(addr[i]));
      for (j = 0; j < BLOCK_DATA_SIZE/2; ++j) {
        fout << "(" << addx[j].x << ", " << addx[j].key << ") ";
      }
      if (i % 5 == 0) {
        fout << std::endl;
      }
    }
  }
  printf("\n");
  fout << std::endl;
  // std::cout<< "\nFinished printEnc. \n";
  fout.close();
}

// TODO: change nt types
void test(int64_t **arrayAddr, int64_t structureId, int64_t size) {
  int64_t pass = 1;
  int64_t i;
  // print(structureId);
  if(structureSize[structureId] == 8) {
    for (i = 1; i < size; i++) {
      pass &= ((arrayAddr[structureId])[i-1] <= (arrayAddr[structureId])[i]);
      if (!pass) {
        // // // std::cout<< (arrayAddr[structureId])[i-1] << ' ' << (arrayAddr[structureId])[i];
        break;
      }
    }
  } else {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (i = 1; i < size; i++) {
      pass &= (addr[i-1].x <= addr[i].x);
      if (!pass) {
        // // // std::cout<< addr[i-1].x << ' ' << addr[i].x;
        break;
      }
    }
  }
  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}

void testEnc(int64_t **arrayAddr, int64_t structureId, int64_t size) {
  int64_t pass = 1;
  int64_t i, j, blockNumber;
  EncBlock *addr = (EncBlock*)arrayAddr[structureId];
  if(structureSize[structureId] == 4) {
    // int64_t type
    blockNumber = (int64_t)ceil(1.0*size/BLOCK_DATA_SIZE);
    for (i = 0; i < blockNumber; i++) {
      int64_t *addx = (int64_t*)(&(addr[i]));
      for (j = 0; j < BLOCK_DATA_SIZE-1; ++j) {
        pass &= (addx[j] <= addx[j+1]);
      }
    }
  } else {
    // Bucket Type
    blockNumber = (int64_t)ceil(1.0*size/(BLOCK_DATA_SIZE/2));
    for (i = 0; i < blockNumber; i++) {
      Bucket_x *addx = (Bucket_x*)(&(addr[i]));
      pass &= (addx[0].x <= addx[1].x);
    }
  }
  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}

void testWithDummy(int64_t **arrayAddr, int64_t structureId, int64_t size) {
  int64_t i = 0;
  int64_t j = 0;
  // print(structureId);
  if(structureSize[structureId] == 8) {
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