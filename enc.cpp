#include <mbedtls/aes.h>
#include "mbedtls/entropy.h"
#include "mbedtls/ctr_drbg.h"
#include <iostream>
#include <set>

mbedtls_aes_context aes;
unsigned char key[32];

void prf(unsigned char *right, char key, int tweak, unsigned char *ret) {
  unsigned char input[16] = {0};
  unsigned char encrypt_output[16] = {0};
  input[0] = right[0];
  input[1] = right[1];
  input[15] = tweak & 0xFF;
  mbedtls_aes_crypt_ecb(&aes, MBEDTLS_AES_ENCRYPT, input, encrypt_output);
  ret[0] = encrypt_output[0];
  ret[1] = encrypt_output[1];
}

void round(unsigned char* data, char key, int tweak, unsigned char* newData) {
  unsigned char leftBits[2] = {data[0], data[1]};
  unsigned char rightBits[2] = {data[2], data[3]};
  newData[0] = rightBits[0];
  newData[1] = rightBits[1];
  unsigned char prfRet[2];
  prf(rightBits, key, tweak, prfRet);
  newData[2] = leftBits[0] ^ prfRet[0];
  newData[3] = leftBits[1] ^ prfRet[1];
}

// char key[8] = "=-÷×&";
int encrypt(int index, unsigned char *key, int rounds) {
  unsigned char bytes[4];
  for (int i = 0; i < 4; ++i) {
    bytes[i] = (index >> (24 - i * 8)) & 0xFF;
  }
  int keyIdx = 0;
  unsigned char newBytes[4];
  while (rounds--) {
    round(bytes, key[keyIdx], rounds, newBytes);
    keyIdx = (keyIdx + 1) % 32;
  }
  int newIndex = 0;
  for (int i = 0; i < 4; ++i) {
    newIndex |= newBytes[i] << (24 - i * 8);
  }
  return newIndex;
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
  if((ret = mbedtls_ctr_drbg_random(&ctr_drbg, key, 32)) != 0) {
    printf(" failed\n ! mbedtls_ctr_drbg_random returned -0x%04x\n", -ret);
    return -1;
  }
  mbedtls_aes_setkey_enc(&aes, key, 256);
  std::set<int> hash;
  int i = 0, res;
  while (hash.size() < 2000000) {
    res = encrypt(i, key, 3);
    if (hash.count(res)) {
      std::cout << "Duplicate!\n";
      return 0;
    }
    if (res >= 0 && res < 2000000 && !hash.count(res)) {
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