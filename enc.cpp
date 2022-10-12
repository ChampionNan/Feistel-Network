#include <mbedtls/aes.h>
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"
#include <iostream>
#include <set>
#include <cmath>

mbedtls_aes_context aes;
unsigned char key[16];
int base;
int max_num = 2000000;
int _max_num;

__int128_t prf(__int128_t a) {
  // std::cout << "In prf\n";
  unsigned char input[16] = {0};
  unsigned char encrypt_output[16] = {0};
  for (int i = 0; i < 16; ++i) {
    input[i] = (a >> (120 - i * 8)) & 0xFF;
  }
  // input[15] = tweak & 0xFF;
  mbedtls_aes_crypt_ecb(&aes, MBEDTLS_AES_ENCRYPT, input, encrypt_output);
  __int128_t res = 0;
  // int flag = (1 << base+1) - 1;
  for (int i = 0; i < 16; ++i) {
    res |= encrypt_output[i] << (120 - i * 8);
  }
  return res;
}

// char key[8] = "=-÷×&";
int encrypt(int index, unsigned char *key, int rounds) {
  // std::cout << "In Encrypt\n";
  int l = index / (1 << base);
  int r = index % (1 << base);
  __int128_t e;
  int temp, i = 1;
  while (i <= rounds) {
    e = prf((r << 16 * 8 - base) + rounds);
    temp = r;
    r = l ^ (e >> 16 * 8 - base);
    l = temp;
    i += 1;
  }
  return (l << base) + r;
}

int main(int argc, const char* argv[]) {
  mbedtls_aes_init(&aes);
  // generate key
  mbedtls_ctr_drbg_context ctr_drbg;
  mbedtls_entropy_context entropy;
  char *pers = "aes generate key";
  int ret;
  mbedtls_entropy_init(&entropy);
  mbedtls_ctr_drbg_init(&ctr_drbg);
  if((ret = mbedtls_ctr_drbg_seed(&ctr_drbg, mbedtls_entropy_func, &entropy, (unsigned char *) pers, strlen(pers))) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_init returned -0x%04x\n", -ret);
    return -1;
  }
  if((ret = mbedtls_ctr_drbg_random(&ctr_drbg, key, 16)) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_random returned -0x%04x\n", -ret);
    return -1;
  }
  mbedtls_aes_setkey_enc(&aes, key, 128);
  base = ceil(1.0 * log2(max_num) / 2); //=width/2
  _max_num = 1 << 2 * base;

  std::set<int> hash;
  int i = 0, res;
  std::cout << "Max_num: " << _max_num << std::endl;
  while (hash.size() < _max_num) {
    res = encrypt(i, key, 3);
    if (hash.count(res)) {
      // std::cout << "Duplicate keys!\n";
      continue;
    }
    if (res >= 0 && res < 2000000) {
      hash.insert(res);
    }
    if (i % 1000000 == 0 && i != 0) {
      std::cout << "i: " << i << std::endl;
      std::cout << "Hash size: " << hash.size() << std::endl;
    }
    i += 1;
  }

  mbedtls_aes_free(&aes);
  printf("Encrypt size: %d\n", hash.size());
  printf("Max element: %d\n", *hash.begin());
  printf("Min element: %d\n", *hash.rbegin());
  return 0;
}