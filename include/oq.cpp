# include "oq.h"
#include "common.h"
// OQSORT
extern int is_tight;
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
  printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::pair<int, int> section= OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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
  printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  free(trustedM);
  return outStructureId2;
}

// TODO: TightSort2
int ObliviousTightSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2) {
  std::cout << "In ObliviousTightSort2 && In SampleRec\n";
  std::vector<std::vector<int> > pivots;
  SampleRec(inStructureId, sampleId, sortedSampleId, 1, pivots);
  printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::cout << "In TwoLevelPartition\n";
  std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
  printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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
  printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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
  printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::pair<int, int> section = OneLevelPartition(inStructureId, inSize, trustedM2, realNum, P, outStructureId1, 0);
  printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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
  printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  return {outStructureId2, totalLevelSize};
}

// TODO: Finish
std::pair<int, int> ObliviousLooseSort2(int inStructureId, int inSize, int sampleId, int sortedSampleId, int outStructureId1, int outStructureId2) {
  std::cout << "In ObliviousLooseSort2 && In SasmpleRec\n";
  std::vector<std::vector<int> > pivots;
  SampleRec(inStructureId, sampleId, sortedSampleId, 0, pivots);
  printf("Till Sample IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
  std::cout << "In TwoLevelPartition\n";
  std::pair<int, int> section = TwoLevelPartition(inStructureId, pivots, P, outStructureId1, outStructureId2);
  printf("Till Partition IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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
  printf("Till Final IOcost: %f\n", 1.0*IOcost/N*BLOCK_DATA_SIZE);
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