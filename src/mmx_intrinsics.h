/*!
 *  \file mmx_intrinsics.h
 *  \brief      Provides C functions that implement Intel's MMX intrinsic
 *  \details    Based off VOLK kernel functions.
 *  \author    Jake Johnson
 *  \version   4.1a
 *  \date      Jan 23, 2018
 *  \pre       Make sure you have .bin files containing data and lookup tables
 */

#include "mmintrin.h"
#include <stdint.h>
#include <stdio.h>

/*!
 *  \brief     Multiply and accumulate two vectors of int16_t integers with MMX
 * intrinsic functions \details   This class is used to demonstrate a number of
 * section commands. \param[out] int32_t returnValue as result of multiplication
 * and accumulation \param[in] aVector Source vector with factors to multiply
 *  \param[in] bVector Source vector with factors to multiply
 *  \param[in] num_points Number of points to Multiply in the operation
 */
static inline double mmx_mul_and_acc_short(const int16_t *aVector,
                                           const int16_t *bVector,
                                           uint32_t num_points) {

  int32_t returnValue = 0;
  uint32_t number = 0;
  const uint32_t quarterPoints = num_points / 4;

  const int16_t *aPtr = aVector;
  const int16_t *bPtr = bVector;
  int16_t tempBuffer[4];

  __m64 aVal, bVal, cVal;
  __m64 accumulator = _mm_setzero_si64();

  for (; number < quarterPoints; number++) {

    // Load 256-bits of integer data from memory into dst. mem_addr does not
    // need to be aligned on any particular boundary.
    aVal = _m_from_int64(*aPtr);
    bVal = _m_from_int64(*bPtr);

    cVal = _mm_mullo_pi16(aVal, bVal);

    accumulator = _mm_adds_pi16(accumulator, cVal);

    // Increment pointers
    aPtr += 4;
    bPtr += 4;
  }

  //_mm256_storeu_si256((__m256i*)tempBuffer, accumulator);
  *tempBuffer = _m_to_int64(accumulator);

  returnValue = tempBuffer[0];
  returnValue += tempBuffer[1];
  returnValue += tempBuffer[2];
  returnValue += tempBuffer[3];

  // Perform non SIMD leftover operations
  number = quarterPoints * 4;
  for (; number < num_points; number++) {
    returnValue += (*aPtr++) * (*bPtr++);
  }
  return returnValue;
}
