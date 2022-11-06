#include "common.h"

std::random_device dev2;
std::mt19937 rng2(dev2());
unsigned char key2[16];
mbedtls_aes_context aes2;
mbedtls_ctr_drbg_context ctr_drbg;
mbedtls_entropy_context entropy;


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
  int ret;
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
void cbc_encrypt(EncBlock* buffer, int blockSize) {
  // std::cout << "In cbc_encrypt\n";
  mbedtls_ctr_drbg_random(&ctr_drbg, (uint8_t*)(&(buffer->iv)), 16);
  unsigned char iv[16];
  memcpy(iv, (uint8_t*)(&(buffer->iv)), blockSize);
  mbedtls_aes_crypt_cfb8(&aes2, MBEDTLS_AES_ENCRYPT, blockSize, (uint8_t*)(&(buffer->iv)), (uint8_t*)buffer, (uint8_t*)buffer);
  memcpy((uint8_t*)(&(buffer->iv)), iv, blockSize);
  return;
}

// Assume blockSize = 16 * k
void cbc_decrypt(EncBlock* buffer, int blockSize) {
  // std::cout << "In cbc_decrypt\n";
  // unsigned char iv[16];
  // memcpy(iv, (uint8_t*)(&(buffer->iv)), blockSize);
  mbedtls_aes_crypt_cfb8(&aes2, MBEDTLS_AES_DECRYPT, blockSize, (uint8_t*)(&(buffer->iv)), (uint8_t*)buffer, (uint8_t*)buffer);
  // memcpy(iv, (uint8_t*)(&(buffer->iv)), blockSize);
  // int *a = (int*)buffer;
  // printf("%d, %d, %d, %d\n", a[0], a[1], a[2], a[3]);
  return;
}

// index: EncBlock, offset: int, blockSize: bytes
void OcallRB(int index, int* buffer, int blockSize, int structureId) {
  memcpy(buffer, &(((EncBlock*)(arrayAddr[structureId]))[index]), blockSize); 
  IOcost += 1;
}

// startIdx: index of blocks, offset: data offset in the block, blockSize: bytes of real data
void OcallReadBlock(int startIdx, int offset, int* buffer, int blockSize, int structureId) {
  if (blockSize == 0) {
    return;
  }
  EncBlock *readBuffer = (EncBlock*)malloc(sizeof(EncBlock));
  if (nonEnc) {
    OcallRB(startIdx, (int*)readBuffer, sizeof(EncBlock), structureId);
    int *a = (int*)readBuffer;
    memcpy(buffer, (int*)readBuffer+offset*(structureSize[structureId]/sizeof(int)), blockSize);
  } else {
    OcallRB(startIdx, (int*)readBuffer, sizeof(EncBlock), structureId);
    cbc_decrypt(readBuffer, sizeof(EncBlock)/2);
    memcpy(buffer, (int*)readBuffer+offset*(structureSize[structureId]/sizeof(int)), blockSize);
  }
  free(readBuffer);
}

// index: EncBlock, offset: int, blockSize: bytes(only non-enc mode could use)
void OcallWB(int index, int offset, int* buffer, int blockSize, int structureId) {
  memcpy((int*)(&(((EncBlock*)(arrayAddr[structureId]))[index]))+offset, buffer, blockSize);
  IOcost += 1;
}

// startIdx: index of blocks, offset(int): data offset in the block, blockSize: bytes of real data
void OcallWriteBlock(int startIdx, int offset, int* buffer, int blockSize, int structureId) {
  if (blockSize == 0) {
    return;
  }
  EncBlock* writeBuf = (EncBlock*)malloc(sizeof(EncBlock));
  if (nonEnc) {
    OcallWB(startIdx, offset*(structureSize[structureId]/sizeof(int)), buffer, blockSize, structureId);
  } else {
    if (offset == 0) { // could write the whole block
      memcpy((int*)writeBuf, buffer, blockSize);
      cbc_encrypt(writeBuf, sizeof(EncBlock)/2);
      OcallWB(startIdx, 0, (int*)writeBuf, sizeof(EncBlock), structureId);
    } else { // read&decrypt first, later write&encrypt
      OcallRB(startIdx, (int*)writeBuf, sizeof(EncBlock), structureId);
      cbc_decrypt(writeBuf, sizeof(EncBlock)/2);
      memcpy((int*)writeBuf+offset*(structureSize[structureId]/sizeof(int)), buffer, blockSize);
      cbc_encrypt(writeBuf, sizeof(EncBlock)/2);
      OcallWB(startIdx, 0, (int*)writeBuf, sizeof(EncBlock), structureId);
    }
  }
  free(writeBuf);
}

// TODO: Set this function as OCALL
// allocate encB for outside memory
void freeAllocate(int structureIdM, int structureIdF, int size) {
  // 1. Free arrayAddr[structureId]
  if (arrayAddr[structureIdF]) {
    free(arrayAddr[structureIdF]);
  }
  // 2. malloc new asked size (allocated in outside)
  if (size <= 0) {
    return;
  }
  int *addr = (int*)malloc(size * sizeof(EncBlock));
  memset(addr, DUMMY, size * sizeof(EncBlock));
  // 3. assign malloc address to arrayAddr
  arrayAddr[structureIdM] = addr;
  return ;
}

// index: start index counted by elements (count from 0), blockSize: #elements
void opOneLinearScanBlock(int index, int* block, int blockSize, int structureId, int write, int dummyNum=0) {
  if (blockSize + dummyNum == 0) {
    return ;
  }
  if (dummyNum < 0) {
    printf("Dummy padding error!");
    return ;
  }
  int multi = structureSize[structureId] / sizeof(int);
  int B = BLOCK_DATA_SIZE / multi;
  int startBIdx = index / B;
  int startOffset = index % B;
  int firstSize = B - startOffset;
  int boundary = ceil(1.0 * (blockSize - firstSize) / B) + 1;
  int opStart = 0;
  int Msize, offset = 0;
  if (!write) {
    std::cout << "In reading B: \n";
    for (int i = 0; i < boundary; ++i) {
      if (i != 0 && (i != boundary - 1)) {
        Msize = B;
      } else if (i == 0) {
        Msize = firstSize;
      } else {
        Msize = blockSize - opStart;
      }
      offset = (i == 0) ? startOffset : 0;
      printf("Start Idx: %d, offset: %d, opstart: %d, Msize: %d\n", startBIdx + i, offset, opStart, Msize);
      OcallReadBlock(startBIdx + i, offset, &block[opStart * multi], Msize * structureSize[structureId], structureId);
      opStart += Msize;
    }
  } else {
    for (int i = 0; i < boundary; ++i) {
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
      int *junk = (int*)malloc(dummyNum * multi * sizeof(int));
      for (int j = 0; j < dummyNum * multi; ++j) {
        junk[j] = DUMMY;
      }
      int dummyStart = (index + blockSize) / B;
      int dummyOffset = (index + blockSize) % B;
      int dummyFirstSize = B - dummyOffset;
      int dummyBoundary = ceil(1.0 * (dummyNum - dummyFirstSize) / B) + 1;
      int dummyOpStart = 0;
      for (int j = 0; j < dummyBoundary; ++j) {
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
  int j, t;
  if (n > 1) {
    for (int i = 0; i < n - 1; ++i) {
      // TODO: shuffle error
      j = i + rand() / (INT_MAX / (n - i) + 1);
      t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
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
  std::uniform_int_distribution<int> dist{low, high};
  int randNum = dist(rng2);
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
  std::uniform_int_distribution<int> dist{low, high};
  int randNum = dist(rng2);
  // int randNum = rand() % (high - low + 1) + low;
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
    // printf("In quickSort: %ld, %ld\n", low, high);
    int mid = partition(arr, low, high);
    quickSort(arr, low, mid - 1);
    quickSort(arr, mid + 1, high);
  }
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

// SUPPORT
void print(int* array, int size) {
  int i;
  for (i = 0; i < size; i++) {
    printf("%d ", array[i]);
    if ((i != 0) && (i % 6 == 0)) {
      printf("\n");
    }
  }
  printf("\n");
}

void print(int **arrayAddr, int structureId, int size) {
  int i;
  std::ofstream fout("/home/data/bchenba/output.txt");
  if(structureSize[structureId] == 8) {
    int *addr = (int*)arrayAddr[structureId];
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
void printEnc(int **arrayAddr, int structureId, int size) {
  int i, j, blockNumber;
  std::ofstream fout("/home/data/bchenba/res.txt");
  EncBlock *addr = (EncBlock*)arrayAddr[structureId];
  std::cout << "In printEnc: \n";
  if(structureSize[structureId] == 4) {
    blockNumber = (int)ceil(1.0*size/4);
    for (i = 0; i < blockNumber; i++) {
      int *addx = (int*)(&(addr[i]));
      printf("%d: (", i);
      for (j = 0; j < 4; ++j) {
        fout << addx[j] << ' ';
        std::cout << addx[j] << ' '; 
      }
      std::cout << ")\n";
      std::cout << "(";/*
      for (j = 5; j < 8; ++j) {
        fout << addx[j] << ' ';
        std::cout << addx[j] << ' '; 
      }
      std::cout << ")\n";*/
      if (i % 5 == 0) {
        fout << std::endl;
        std::cout << std::endl;
      }
    }
  } else if (structureSize[structureId] == 8) {
    blockNumber = (int)ceil(1.0*size/2);
    for (i = 0; i < blockNumber; i++) {
      Bucket_x *addx = (Bucket_x*)(&(addr[i]));
      for (j = 0; j < 2; ++j) {
        fout << "(" << addx[j].x << ", " << addx[j].key << ") ";
      }
      if (i % 5 == 0) {
        fout << std::endl;
      }
    }
  }
  // printf("\n");
  fout << std::endl;
  std::cout << "\nFinished printEnc. \n";
  fout.close();
}

// TODO: change nt types
void test(int **arrayAddr, int structureId, int size) {
  int pass = 1;
  int i;
  // print(structureId);
  if(structureSize[structureId] == 8) {
    for (i = 1; i < size; i++) {
      pass &= ((arrayAddr[structureId])[i-1] <= (arrayAddr[structureId])[i]);
      if (!pass) {
        std::cout << (arrayAddr[structureId])[i-1] << ' ' << (arrayAddr[structureId])[i];
        break;
      }
    }
  } else {
    Bucket_x *addr = (Bucket_x*)arrayAddr[structureId];
    for (i = 1; i < size; i++) {
      pass &= (addr[i-1].x <= addr[i].x);
      if (!pass) {
        std::cout << addr[i-1].x << ' ' << addr[i].x;
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