#include "bitonic.h"

void smallBitonicMerge(int64_t *a, int64_t start, int64_t size, int64_t flipped) {
  if (size == 1) {
    return;
  } else {
    int64_t swap = 0;
    int64_t mid = greatestPowerOfTwoLessThan((double)size);
    for (int64_t i = 0; i < size - mid; ++i) {
      int64_t num1 = a[start + i];
      int64_t num2 = a[start + mid + i];
      swap = num1 > num2;
      swap = swap ^ flipped;
      a[start + i] = (!swap * num1) + (swap * num2);
      a[start + i + mid] = (swap * num1) + (!swap * num2);
    }
    smallBitonicMerge(a, start, mid, flipped);
    smallBitonicMerge(a, start + mid, size - mid, flipped);
  }
  return;
}

//Correct, after testing
void smallBitonicSort(int64_t *a, int64_t start, int64_t size, int64_t flipped) {
  if (size <= 1) {
    return ;
  } else {
    int64_t mid = greatestPowerOfTwoLessThan((double)size);
    smallBitonicSort(a, start, mid, 1);
    smallBitonicSort(a, start + mid, size - mid, 0);
    smallBitonicMerge(a, start, size, flipped);
  }
  return;
}
