# include "oq.h"
#include "common.h"
// OQSORT
extern int is_tight;
mbedtls_aes_context aes;
unsigned char key[16];
int base;
int max_num;
int ROUND = 3;
int fail = 0;
std::random_device dev;
std::mt19937 rng(dev()); 

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

int sumArray(bool *M, int left, int right) {
  int total = 0;
  for (int i = left; i < right; ++i) {
    total += M[i];
  }
  return total;
}

void swapArray(int *D, int i, int j) {
  int temp = D[i];
  D[i] = D[j];
  D[j] = temp;
}

void OROffCompact(int *D, bool *M, int left, int right, int z) {
  int n = right - left;
  if (n == 1) return;
  int n2 = n / 2;
  int m = sumArray(M, left, left + n2);
  if (n == 2) {
    bool flag = ((1 - M[left]) * M[left + 1]) ^ (z % 2);
    if (flag) swapArray(D, left, left + 1);
  } else if (n > 2) {
    OROffCompact(D, M, left, left + n2, z % n2);
    OROffCompact(D, M, left + n2, left + n, (z + m % n2) % n2);
    bool s1 = (((z % n2) + m) >= n2) ? 1 : 0;
    bool s2 = (z >= n2) ? 1 : 0;
    bool s = s1 ^ s2;
    for (int i = 0; i < n2; ++i) {
      bool s3 = (i >= ((z + m) % n2)) ? 1 : 0;
      bool b = s ^ s3;
      if (b) swapArray(D, left + i, left + i + n2);
    }
  }
}

void ORCompact(int *D, bool *M, int left, int right) {
  int n = right - left;
  if (n <= 1) return;
  int n1 = pow(2, floor(log2(n)));
  int n2 = n - n1;
  int m = sumArray(M, left, left + n2);
  ORCompact(D, M, left, left + n2);
  OROffCompact(D, M, left + n2, left + n, (n1 - n2 + m) % n1);
  // print(D, n);
  for (int i = 0; i < n2; ++i) { // not 2^k mix
    bool b = (i >= m) ? 1 : 0;
    if (b) swapArray(D, left + i, left + i + n1);
  }
}

// #elem < pivot
int assignM(int *arr, bool *M, int left, int right, int pivot) {
  int total = 0;
  for (int i = left; i < right; ++i) {
    M[i] = (arr[i] < pivot) ? 1 : 0;
    total += M[i];
  }
  return total;
}

void multiLevelPartiton(int *D, bool *M, int low, int high, std::vector<int> pivots,int left, int right, std::vector<int> &partitionIdx) {
  int realNum, dummyNum;
  int pivotIdx, pivot, leftReal, totalLeftNum, leftDummy;
  if (right >= left) {
    // 1. cal pivot
    pivotIdx = (left + right) >> 1;
    pivot = pivots[pivotIdx];
    // 2. assign flag to dummy
    realNum = moveDummy(D, high - low);
    dummyNum = high - low - realNum;
    leftReal = assignM(D, M, low, low + realNum, pivot);
    totalLeftNum = (realNum + dummyNum) >> 1;
    leftDummy = totalLeftNum-leftReal;
    print(&D[low], high-low);
    memset(&M[low+realNum], 1, leftDummy * sizeof(bool));
    memset(&M[low+realNum+leftDummy], 0, (dummyNum-leftDummy) * sizeof(bool));
    // 3. call
    ORCompact(D, M, low, high);
    partitionIdx.push_back(low+leftReal+leftDummy);
    // 4. recusive logP times
    multiLevelPartiton(D, M, low, low+leftReal+leftDummy, pivots, left, pivotIdx-1, partitionIdx);
    multiLevelPartiton(D, M, low+leftReal+leftDummy, high, pivots, pivotIdx+1, right, partitionIdx);
  }
}

// call obliviousPartition
void internalObliviousSort(int *D, int left, int right, int smallM, double beta, double gamma) {
  int M = right - left;
  int hatM = M + ceil((beta + gamma) * M);
  int n = ceil(1.0 * ALPHA * M);
  int r = ceil(log2(1.0 * M / smallM));
  int p = pow(2, r) - 1; // we need 2^r-1 pivots
  // 1. get D's samples
  vector<int> samples;
  floydSampler(M, n, samples);
  for (int i = 0; i < samples.size(); ++i) {
    samples[i] = D[samples[i]]; // replace index to real value
  }
  sort(samples.begin(), samples.end());
  // 2. get pivots
  quantileCal2(samples, 0, samples.size(), p);
  print(samples);
  vector<int> partitionIdx;
  bool *indicator = (bool*)malloc(sizeof(bool) * hatM);
  // copy D to extended memory
  int *extD = (int*)malloc(sizeof(int) * hatM);
  memset(extD, DUMMY, sizeof(int) * hatM);
  memcpy(extD, D, sizeof(int) * M);
  vector<int> size0; // record for level even level
  vector<int> size1; // record for level odd level
  // 3. get partitioned index
  // 3.1 level 0
  int pivot0 = samples[p/2]; // pivot for level 0
  int leftRealNum = assignM(extD, indicator, left, right, pivot0);
  memset(&indicator[right], 1, sizeof(bool) * (hatM/2-leftRealNum));
  memset(&indicator[right+hatM/2-leftRealNum], 0, sizeof(bool) * (hatM-hatM/2-(M-leftRealNum)));
  ORCompact(extD, indicator, 0, hatM);
  size0.push_back(hatM/2);
  size0.push_back(hatM-hatM/2);
  // 3.2 > level 1
  int readStart, secSize, realNum, dummyNum, pivot, rightRealNum;
  for (int i = 1; i < r; ++i) { // index of level
    readStart = 0; // record read start for each section
    vector<int> readSize, writeSize;
    if (i % 2 == 1) { // readSize in size0
      size1.clear();
      readSize = size0;
    } else { // readSize in size1
      size0.clear();
      readSize = size1;
    }
    for (int j = 0; j < pow(2, i); ++j) { // j: #last level sections
      realNum = moveDummy(&extD[readStart], readSize[j]);
      dummyNum = readSize[j] - realNum;
      pivot = samples[(2*j+1)*p/pow(2, i+1)];
      leftRealNum = assignM(extD, indicator, readStart, readStart+realNum, pivot);
      rightRealNum = realNum - leftRealNum;
      memset(&indicator[readStart+realNum], 1, sizeof(bool) * (readSize[j]/2-leftRealNum));
      memset(&indicator[readStart+realNum+(readSize[j]/2-leftRealNum)], 0, sizeof(bool) * (dummyNum-(readSize[j]/2-leftRealNum)));
      ORCompact(extD, indicator, readStart, readStart+readSize[j]);
      writeSize.push_back(readSize[j]/2);
      writeSize.push_back(readSize[j]-readSize[j]/2);
      // each section ending part
      readStart += readSize[j];
    }
    if (i % 2 == 1) { // writeSize in size1
      size1 = writeSize;
    } else { // writeSize in size0
      size0 = writeSize;
    }
  }
  if ((r - 1) % 2 == 1) { // writeSize in size1
    partitionIdx = size1;
  } else { // writeSize in size0
    partitionIdx = size0;
  }
  for (int i = 1; i < partitionIdx.size(); ++i) {
    partitionIdx[i] += partitionIdx[i-1];
  }
  // 4. get the right order of index
  partitionIdx.insert(partitionIdx.begin(), 0);
  // partitionIdx.push_back(hatM);
  // 5. for each section, call bitonic sort, [)
  int writeStart = 0;
  for (int i = 0; i < partitionIdx.size() - 1; ++i) {
    printf("Begin at: %d, number: %d\n", partitionIdx[i], partitionIdx[i+1]-partitionIdx[i]);
    realNum = moveDummy(&extD[partitionIdx[i]], partitionIdx[i+1]-partitionIdx[i]);
    smallBitonicSort(extD, partitionIdx[i], realNum, 0);
    memcpy(&D[writeStart], &extD[partitionIdx[i]], sizeof(int) * realNum);
    writeStart += realNum;
  }
  free(indicator);
  free(extD);
  print(D, M);
}

void floydSampler(int n, int k, std::vector<int> &x) {
  std::unordered_set<int> H;
  for (int i = n - k; i < n; ++i) {
    x.push_back(i);
  }
  std::random_device dev;
  std::mt19937 rng(dev()); 
  int r, j, temp;
  for (int i = 0; i < k; ++i) {
    std::uniform_int_distribution<int> dist{0, n-k+1+i};
    r = dist(rng); // get random numbers with PRNG
    if (H.count(r)) {
      std::uniform_int_distribution<int> dist2{0, i};
      j = dist2(rng);
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

int Sample(int inStructureId, int sampleSize, std::vector<int> &trustedM2, int is_tight, int is_rec=0) {
  printf("In sample\n");
  int N_prime = sampleSize;
  // double alpha = (!is_rec) ? ALPHA : _ALPHA;
  double alpha = ALPHA;
  int n_prime = ceil(1.0 * alpha * N_prime);
  int boundary = ceil(1.0 * N_prime / BLOCK_DATA_SIZE);
  int j = 0, Msize;
  int *trustedM1 = (int*)malloc(BLOCK_DATA_SIZE * sizeof(int));
  nonEnc = 1;
  std::vector<int> sampleIdx;
  floydSampler(N_prime, n_prime, sampleIdx);
  for (int i = 0; i < boundary; ++i) {
    if (is_tight) {
      Msize = std::min((int)BLOCK_DATA_SIZE, N_prime - i * BLOCK_DATA_SIZE);
      opOneLinearScanBlock(i * BLOCK_DATA_SIZE, trustedM1, Msize, inStructureId, 0, 0);
      while ((j < n_prime) && (sampleIdx[j] >= i * BLOCK_DATA_SIZE) && (sampleIdx[j] < (i+1) * BLOCK_DATA_SIZE)) {
        trustedM2.push_back(trustedM1[sampleIdx[j] % BLOCK_DATA_SIZE]);
        j += 1;
      }
    } else if ((!is_tight) && (sampleIdx[j] >= i * BLOCK_DATA_SIZE) && (sampleIdx[j] < (i+1) * BLOCK_DATA_SIZE)) {
      Msize = std::min((int)BLOCK_DATA_SIZE, N_prime - i * BLOCK_DATA_SIZE);
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

std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec) {
  if (inSize <= M) {
    return {inSize, 1};
  }
  printf("In OneLevelPartition\n");
  // double beta = (!is_rec) ? BETA : _BETA;
  double beta = BETA;
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
  int multi = structureSize[outStructureId1]/sizeof(int);
  int totalEncB = ceil(1.0 * boundary1 * smallSectionSize * p0 / (BLOCK_DATA_SIZE / multi));
  freeAllocate(outStructureId1, outStructureId1, totalEncB);
  
  int Msize1, Msize2, index1, index2, writeBackNum;
  int total_blocks = ceil(1.0 * inSize / BLOCK_DATA_SIZE);
  int *trustedM3 = (int*)malloc(sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
  std::vector<int> partitionIdx;
  // OCall
  fyShuffle(inStructureId, inSize, BLOCK_DATA_SIZE);
  // Finish FFSEM implementation in c++
  // pseudo_init(total_blocks);
  // printf("Before partition\n");
  for (int i = 0; i < boundary1; ++i) {
    // Read one M' memory block after fisher-yates shuffle
    printf("OneLevel %d/%d\n", i, boundary1);
    Msize1 = std::min(boundary2 * BLOCK_DATA_SIZE, inSize - i * boundary2 * BLOCK_DATA_SIZE);
    nonEnc = 1;
    opOneLinearScanBlock(i * boundary2 * BLOCK_DATA_SIZE, trustedM3, Msize1, inStructureId, 0, 0);

    int blockNum = moveDummy(trustedM3, dataBoundary);
    quickSortMulti(trustedM3, 0, blockNum-1, samples, 1, p0, partitionIdx);
    sort(partitionIdx.begin(), partitionIdx.end());
    partitionIdx.insert(partitionIdx.begin(), -1);
    for (int j = 0; j < p0; ++j) {
      index1 = partitionIdx[j]+1;
      index2 = partitionIdx[j+1];
      writeBackNum = index2 - index1 + 1;
      if (writeBackNum > smallSectionSize) {
        printf("Overflow in small section M/p0: %d > %d\n", writeBackNum, smallSectionSize);
      }
      nonEnc = 0;
      opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
    }
    memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
    partitionIdx.clear();
  }
  free(trustedM3);
  // mbedtls_aes_free(&aes);
  if (bucketSize0 > M) {
    printf("Each section size is greater than M, adjst parameters: %d, %d", bucketSize0, M);
  }
  return {bucketSize0, p0};
}

int ObliviousTightSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2) {
  int *trustedM;
  printf("In ObliviousTightSort\n");
  aes_init();
  if (inSize <= M) {
    trustedM = (int*)malloc(sizeof(int) * M);
    nonEnc = 1;
    opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
    quickSort(trustedM, 0, inSize - 1);
    int multi = structureSize[outStructureId1]/sizeof(int);
    int totalEncB = ceil(1.0 * inSize / (BLOCK_DATA_SIZE / multi));
    freeAllocate(outStructureId1, outStructureId1, totalEncB);
    nonEnc = 0;
    opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
    free(trustedM);
    return outStructureId1;
  }
  std::vector<int> trustedM2;
  int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
  // printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::pair<int, int> section= OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  // printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  printf("After onelevel\n");
  int sectionSize = section.first;
  int sectionNum = section.second;
  // TODO: IN order to reduce memory, can replace outStructureId2 with inStructureId
  int multi = structureSize[outStructureId2]/sizeof(int);
  int totalEncB = ceil(1.0 * inSize / (BLOCK_DATA_SIZE / multi)); 
  freeAllocate(outStructureId2, outStructureId2, totalEncB);
  trustedM = (int*)malloc(sizeof(int) * M);
  int j = 0, k;
  printf("In final\n");
  for (int i = 0; i < sectionNum; ++i) {
    nonEnc = 0;
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k-1);
    nonEnc = 1;
    opOneLinearScanBlock(j, trustedM, k, outStructureId2, 1, 0);
    j += k;
  }
  // printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  free(trustedM);
  return outStructureId2;
}

std::pair<int, int> ObliviousLooseSort(int inStructureId, int inSize, int outStructureId1, int outStructureId2) {
  printf("In ObliviousLooseSort\n");
  int *trustedM;
  aes_init();
  if (inSize <= M) {
    trustedM = (int*)malloc(sizeof(int) * M);
    nonEnc = 1;
    opOneLinearScanBlock(0, trustedM, inSize, inStructureId, 0, 0);
    quickSort(trustedM, 0, inSize - 1);
    nonEnc = 0;
    opOneLinearScanBlock(0, trustedM, inSize, outStructureId1, 1, 0);
    free(trustedM);
    return {outStructureId1, inSize};
  }
  std::vector<int> trustedM2;
  int realNum = Sample(inStructureId, inSize, trustedM2, is_tight);
  // printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::pair<int, int> section = OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  // printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  int sectionSize = section.first;
  int sectionNum = section.second;
  int totalLevelSize = sectionNum * sectionSize;
  int k;
  int multi = structureSize[outStructureId2]/sizeof(int);
  int totalEncB = ceil(1.0 * totalLevelSize / (BLOCK_DATA_SIZE / multi)); 
  freeAllocate(outStructureId2, outStructureId2, totalEncB);
  trustedM = (int*)malloc(sizeof(int) * M);
  printf("In final\n");
  for (int i = 0; i < sectionNum; ++i) {
    nonEnc = 0;
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId1, 0, 0);
    k = moveDummy(trustedM, sectionSize);
    quickSort(trustedM, 0, k - 1);
    nonEnc = 1;
    opOneLinearScanBlock(i * sectionSize, trustedM, sectionSize, outStructureId2, 1, 0);
  }
  // printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  return {outStructureId2, totalLevelSize};
}

/* Feistel version
std::pair<int, int> OneLevelPartition(int inStructureId, int inSize, std::vector<int> &samples, int sampleSize, int p, int outStructureId1, int is_rec) {
  if (inSize <= M) {
    return {inSize, 1};
  }
  printf("In Onelevel Partition\n");
  // double beta = (!is_rec) ? BETA : _BETA;
  double beta = BETA;
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
  int multi = structureSize[outStructureId1]/sizeof(int);
  int totalEncB = ceil(1.0 * boundary1 * smallSectionSize * p0 / (BLOCK_DATA_SIZE / multi));
  freeAllocate(outStructureId1, outStructureId1, totalEncB);
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
      Msize1 = std::min((int)BLOCK_DATA_SIZE, inSize - read_index * BLOCK_DATA_SIZE);
      nonEnc = 1;
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
        printf("Overflow in small section M/p0: %ld > %ld\n", writeBackNum, smallSectionSize);
      }
      nonEnc = 0;
      opOneLinearScanBlock(j * bucketSize0 + i * smallSectionSize, &trustedM3[index1], writeBackNum, outStructureId1, 1, smallSectionSize - writeBackNum);
    }
    memset(trustedM3, DUMMY, sizeof(int) * boundary2 * BLOCK_DATA_SIZE);
    partitionIdx.clear();
  }
  free(trustedM3);
  free(shuffleB);
  if (bucketSize0 > M) {
    printf("Each section size is greater than M, adjst parameters: %ld, %ld", bucketSize0, M);
  }
  return {bucketSize0, p0};
}*/