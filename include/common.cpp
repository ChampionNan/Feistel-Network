#include "common.h"

std::random_device dev2;
std::mt19937 rng2(dev2());

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
/* OCall functions */
void ocall_print_string(const char *str) {
  /* Proxy/Bridge will check the length and null-terminate
   * the input string to prevent buffer overflow.
   */
  printf("%s", str);
  fflush(stdout);
}
/*
void OcallReadBlock(int index, int* buffer, size_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(buffer, arrayAddr[structureId] + index, blockSize * structureSize[structureId]);
  memcpy(buffer, arrayAddr[structureId] + index, blockSize);
  IOcost += 1.0*blockSize/structureSize[structureId]/BLOCK_DATA_SIZE;
}

void OcallWriteBlock(int index, int* buffer, size_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(arrayAddr[structureId] + index, buffer, blockSize * structureSize[structureId]);
  memcpy(arrayAddr[structureId] + index, buffer, blockSize);
  IOcost += 1.0*blockSize/structureSize[structureId]/BLOCK_DATA_SIZE;
}*/
void OcallReadBlock(int64_t index, int64_t* buffer, int64_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(buffer, arrayAddr[structureId] + index, blockSize * structureSize[structureId]);
  memcpy(buffer, arrayAddr[structureId] + index, blockSize);
  IOcost += 1;
}

void OcallWriteBlock(int64_t index, int64_t* buffer, int64_t blockSize, int structureId) {
  if (blockSize == 0) {
    // printf("Unknown data size");
    return;
  }
  // memcpy(arrayAddr[structureId] + index, buffer, blockSize * structureSize[structureId]);
  memcpy(arrayAddr[structureId] + index, buffer, blockSize);
  IOcost += 1;
}

// TODO: Set this function as OCALL
void freeAllocate(int structureIdM, int structureIdF, int64_t size) {
  // 1. Free arrayAddr[structureId]
  if (arrayAddr[structureIdF]) {
    free(arrayAddr[structureIdF]);
  }
  // 2. malloc new asked size (allocated in outside)
  if (size <= 0) {
    return;
  }
  int64_t *addr = (int64_t*)malloc(size * sizeof(int64_t));
  memset(addr, DUMMY, size * sizeof(int64_t));
  // 3. assign malloc address to arrayAddr
  arrayAddr[structureIdM] = addr;
  return ;
}

// Functions x crossing the enclave boundary, unit: BLOCK_DATA_SIZE
// TODO: FIX This function ? need checking
/*
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
}*/
void opOneLinearScanBlock(int64_t index, int64_t* block, int64_t blockSize, int structureId, int write, int64_t dummyNum=0) {
  if (blockSize + dummyNum == 0) {
    return ;
  }
  if (dummyNum < 0) {
    printf("Dummy padding error!");
    return ;
  }
  int64_t boundary = (int64_t)((blockSize + BLOCK_DATA_SIZE - 1 )/ BLOCK_DATA_SIZE);
  int64_t Msize, i;
  int multi = structureSize[structureId] / sizeof(int64_t);
  if (!write) {
    // OcallReadBlock(index, block, blockSize * structureSize[structureId], structureId);
    for (i = 0; i < boundary; ++i) {
      Msize = std::min((int64_t)BLOCK_DATA_SIZE, blockSize - i * BLOCK_DATA_SIZE);
      OcallReadBlock(index + multi * i * BLOCK_DATA_SIZE, &block[i * BLOCK_DATA_SIZE * multi], Msize * structureSize[structureId], structureId);
    }
  } else {
    // OcallWriteBlock(index, block, blockSize * structureSize[structureId], structureId);
    for (i = 0; i < boundary; ++i) {
      Msize = std::min((int64_t)BLOCK_DATA_SIZE, blockSize - i * BLOCK_DATA_SIZE);
      OcallWriteBlock(index + multi * i * BLOCK_DATA_SIZE, &block[i * BLOCK_DATA_SIZE * multi], Msize * structureSize[structureId], structureId);
    }
    if (dummyNum > 0) {
      int64_t *junk = (int64_t*)malloc(dummyNum * multi * sizeof(int64_t));
      for (int64_t j = 0; j < dummyNum * multi; ++j) {
        junk[j] = DUMMY;
      }
      int64_t startIdx = index + multi * blockSize;
      boundary = ceil(1.0 * dummyNum / BLOCK_DATA_SIZE);
      for (int64_t j = 0; j < boundary; ++j) {
        Msize = std::min((int64_t)BLOCK_DATA_SIZE, dummyNum - j * BLOCK_DATA_SIZE);
        OcallWriteBlock(startIdx + multi * j * BLOCK_DATA_SIZE, &junk[j * BLOCK_DATA_SIZE * multi], Msize * structureSize[structureId], structureId);
      }
    }
  }
  return;
}

// SUPPORTINGs
int myrand() {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<int> dist(0, INT_MAX);
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
/*
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
}*/

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
  // int randNum = rand() % (high - low + 1) + low;
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
    // printf("In quickSort: %lu, %lu\n", low, high);
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
    printf("%lu", array[i]);
    if ((i != 0) && (i % 5 == 0)) {
      printf("\n");
    }
  }
  printf("\n");
}

void print(int64_t **arrayAddr, int structureId, int64_t size) {
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

// TODO: change nt types
void test(int64_t **arrayAddr, int structureId, int64_t size) {
  int pass = 1;
  int64_t i;
  // print(structureId);
  if(structureSize[structureId] == 8) {
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

void testWithDummy(int64_t **arrayAddr, int structureId, int64_t size) {
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