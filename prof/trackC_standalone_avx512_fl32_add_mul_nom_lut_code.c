/*!
 *  \file trackC_standalone_avx512_fl32_add_mul_nom_lut_code.c
 *  \brief      Simulates the tracking stage of a receiver using AVX512
 intrinsics.
 *  \details    Profiles code when using:
 1. Carrier wave generation by means of PLUT method.
 2. Pseudorandom code generation by means of PLUT method.
 3. Down-conversion of the received signal by nominal multiplication.
 4. Multiplication and accumulationof baseband signal with a local replica of
ranging code using AVX512 SIMD intrinsics with fl32 types
 *  \author    Damian Miralles
 *  \author    Jake Johnson
 *  \date      Jan 23, 2018
 *  \pre       Make sure you have .bin files containing data and lookup tables.
 *  \note      Functions in the file must target AVX512 enabled platforms.
 *  \code{.sh}
# Sample compilation script
$ gcc -I ../src/ trackC_standalone_avx512_fl32_add_mul_nom_lut_code.c -g
 -mavx2 -lm -o avx512_fl32_nom_lut_code -O3
 *  \endcode
 */

#include "avx512_intrinsics.h"
#include "read_bin.h"
#include "write_bin.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// Sin and Cos Function approximations
#define gps_sin(x) (((x > 31416) || (x < 0 && x > -31416)) ? -1 : +1)
#define gps_cos(x)                                                             \
  (((x > 15708 && x < 47124) || (x < -15708 && x > -47124)) ? -1 : 1)

int main() {

  // Declarations

  int i, loopcount, blksize;

  int vsmCount, vsmInterval, dataAdaptCoeff;
  double remCodePhase, remCarrPhase, codePhaseStep;
  double earlyLateSpc, seekvalue, samplingFreq, trigarg, carrFreq;
  double I_E, Q_E, I_P, Q_P, I_L, Q_L;
  double carrNco, oldCarrNco, tau1carr, tau2carr, carrError, oldCarrError;
  double PDIcarr, codeNco, oldCodeNco, tau1code, tau2code, codeError,
      oldCodeError, PDIcode;
  double codeFreq, codeFreqBasis, carrFreqBasis, absoluteSample, codeLength;
  double pwr, pwrSum, pwrSqrSum, pwrMean, pwrVar, pwrAvgSqr, pwrAvg, noiseVar,
      CNo, accInt;
  char *rawSignal;
  long int codePeriods;
  const double pi = 3.1415926535;

  FILE *fpdata;

  // Initialization
  remCodePhase = 0;
  remCarrPhase = 0;
  oldCarrNco = 0;
  oldCarrError = 0;
  carrError = 0;
  oldCarrError = 0;
  oldCodeNco = 0;
  oldCodeError = 0;
  absoluteSample = 0;
  vsmCount = 0;
  pwr = 0;
  CNo = 0;
  pwrSum = 0;
  pwrSqrSum = 0;

  // Get all the vectors/integers/strings from .bin files
  float caCode[1025];
  getcaCodeFromFileAsFloat("../data/caCode.bin", caCode);
  blksize = getIntFromFile("../data/blksize.bin");
  codePhaseStep = getDoubleFromFile("../data/codePhaseStep.bin");
  remCodePhase = getDoubleFromFile("../data/remCodePhase.bin");
  earlyLateSpc = getDoubleFromFile("../data/earlyLateSpc.bin");
  samplingFreq = getDoubleFromFile("../data/samplingFreq.bin");
  remCarrPhase = getDoubleFromFile("../data/remCarrPhase.bin");
  carrFreq = getDoubleFromFile("../data/carrFreq.bin");
  char fileid[] = "../data/GPS_and_GIOVE_A-NN-fs16_3676-if4_1304.bin";
  seekvalue = getDoubleFromFile("../data/skipvalue.bin");
  tau1carr = getDoubleFromFile("../data/tau1carr.bin");
  tau2carr = getDoubleFromFile("../data/tau2carr.bin");
  PDIcarr = getDoubleFromFile("../data/PDIcarr.bin");
  carrFreqBasis = getDoubleFromFile("../data/carrFreqBasis.bin");
  tau1code = getDoubleFromFile("../data/tau1code.bin");
  tau2code = getDoubleFromFile("../data/tau2code.bin");
  PDIcode = getDoubleFromFile("../data/PDIcode.bin");
  codeFreq =
      1023002.79220779; // getDoubleFromFile("text_data_files/codeFreq.bin");
  codeFreqBasis =
      1023002.79220779; // getDoubleFromFile("text_data_files/codeFreqBasis.bin");
  codeLength = getDoubleFromFile("../data/codeLength.bin");
  codePeriods = (long int)getIntFromFile("../data/codePeriods.bin");
  dataAdaptCoeff = getIntFromFile("../data/dataAdaptCoeff.bin");
  vsmInterval = getIntFromFile("../data/VSMinterval.bin");
  accInt = getDoubleFromFile("../data/accTime.bin");

  // Declare outputs
  double *carrFreq_output = calloc(codePeriods, sizeof(double));
  double *codeFreq_output = calloc(codePeriods, sizeof(double));
  double *absoluteSample_output = calloc(codePeriods, sizeof(double));
  double *codeError_output = calloc(codePeriods, sizeof(double));
  double *codeNco_output = calloc(codePeriods, sizeof(double));
  double *carrError_output = calloc(codePeriods, sizeof(double));
  double *carrNco_output = calloc(codePeriods, sizeof(double));
  double *I_E_output = calloc(codePeriods, sizeof(double));
  double *I_P_output = calloc(codePeriods, sizeof(double));
  double *I_L_output = calloc(codePeriods, sizeof(double));
  double *Q_E_output = calloc(codePeriods, sizeof(double));
  double *Q_P_output = calloc(codePeriods, sizeof(double));
  double *Q_L_output = calloc(codePeriods, sizeof(double));
  double *VSMIndex = calloc(codePeriods / vsmInterval, sizeof(double));
  double *VSMValue = calloc(codePeriods / vsmInterval, sizeof(double));

  const int lutSize = 256;     // [N=number of bits]
  float sin_LUT_fl32[lutSize]; // our sine wave LUT
  float cos_LUT_fl32[lutSize]; // our sine wave LUT

  // Allocate memory for the signal
  rawSignal = calloc(dataAdaptCoeff * blksize, sizeof(char));

  // Open the file for reading the data and fseek if required
  fpdata = fopen(fileid, "rb");
  fseek(fpdata, dataAdaptCoeff * seekvalue, SEEK_SET);

  // Sine Look-up Table Generation
  for (int i = 0; i < lutSize; ++i) {
    sin_LUT_fl32[i] = (float)(10.0 * sinf(2.0f * pi * (float)i / lutSize));
    cos_LUT_fl32[i] = (float)(10.0 * cosf(2.0f * pi * (float)i / lutSize));
  }

  int sec_count = 0;
  printf("\n*** Running: trackC_standalone_avx512_fl32_add_mul_nom_lut_code "
         "***\n");
  for (loopcount = 0; loopcount < codePeriods; loopcount++) {

    if (loopcount == 1000 * sec_count) {
      printf("  [Completed: %d seconds]\r", sec_count);
      fflush(stdout);
      sec_count += 1;
    }

    I_E = 0;
    Q_E = 0;
    I_P = 0;
    Q_P = 0;
    I_L = 0;
    Q_L = 0;
    i = 0;

    codePhaseStep = codeFreq / samplingFreq;
    blksize = ceil((codeLength - remCodePhase) / codePhaseStep);

    // Create blksize_arr
    double blksize_arr[blksize];
    for (i = 0; i < blksize; i++) {
      blksize_arr[i] = i;
    }

    i = fread(rawSignal, sizeof(char), dataAdaptCoeff * blksize, fpdata);

    // An error check should be added here to see if the required amount of data
    // can be read

    // instantiate vectors

    float mixedcarrSin_vec[blksize];
    float mixedcarrCos_vec[blksize];
    float sin_avx_si32[blksize];
    float cos_avx_si32[blksize];
    float eCode_avx_vec[blksize];
    float lCode_avx_vec[blksize];
    float pCode_avx_vec[blksize];

    // Sine AVX2 NCO Look-up Table Implementation
    avx512_nom_nco_fl32(sin_avx_si32, sin_LUT_fl32, blksize, remCarrPhase,
                        carrFreq, samplingFreq);

    avx512_nom_nco_fl32(cos_avx_si32, cos_LUT_fl32, blksize, remCarrPhase,
                        carrFreq, samplingFreq);

    avx512_nom_code_fl32(eCode_avx_vec, pCode_avx_vec, lCode_avx_vec, caCode,
                         blksize, (float)remCodePhase, (float)codeFreq,
                         (float)samplingFreq);

    // This loop is for parts of code I haven't brought out of loop or haven't
    // figured out how to
    for (i = 0; i < blksize; i++) {
      mixedcarrSin_vec[i] = sin_avx_si32[i] * rawSignal[i];
      mixedcarrCos_vec[i] = cos_avx_si32[i] * rawSignal[i];
    }

    // I_E
    double I_E =
        avx512_mul_and_acc_fl32(eCode_avx_vec, mixedcarrSin_vec, blksize);

    // I_L
    double I_L =
        avx512_mul_and_acc_fl32(lCode_avx_vec, mixedcarrSin_vec, blksize);

    // I_P
    double I_P =
        avx512_mul_and_acc_fl32(pCode_avx_vec, mixedcarrSin_vec, blksize);

    // Q_E
    double Q_E =
        avx512_mul_and_acc_fl32(eCode_avx_vec, mixedcarrCos_vec, blksize);

    // Q_L
    double Q_L =
        avx512_mul_and_acc_fl32(lCode_avx_vec, mixedcarrCos_vec, blksize);

    // Q_P
    double Q_P =
        avx512_mul_and_acc_fl32(pCode_avx_vec, mixedcarrCos_vec, blksize);

    // Compute the VSM C/No
    pwr = I_P * I_P + Q_P * Q_P;
    pwrSum += pwr;
    pwrSqrSum += pwr * pwr;
    vsmCount++;

    if (vsmCount == vsmInterval) {
      pwrMean = pwrSum / vsmInterval;
      pwrVar = pwrSqrSum / vsmInterval - pwrMean * pwrMean;
      pwrAvgSqr = pwrMean * pwrMean - pwrVar;
      pwrAvgSqr = (pwrAvgSqr > 0) ? pwrAvgSqr : -pwrAvgSqr;
      pwrAvg = sqrt(pwrAvgSqr);
      noiseVar = 0.5 * (pwrMean - pwrAvg);
      CNo = (pwrAvg / accInt) / (2 * noiseVar);
      CNo = (CNo > 0) ? CNo : -CNo;
      CNo = 10 * log10(CNo);

      *(VSMIndex + loopcount / vsmInterval) = loopcount + 1;
      *(VSMValue + loopcount / vsmInterval) = CNo;

      vsmCount = 0;
      pwrSum = 0;
      pwrSqrSum = 0;
    }

    remCodePhase = ((remCodePhase) + (blksize)*codePhaseStep - 1023);
    trigarg = (2.0 * pi * carrFreq * (blksize / samplingFreq)) + remCarrPhase;
    remCarrPhase = trigarg - (2 * pi) * ((int)(trigarg / (2 * pi)));

    // Implement carrier loop discriminator (phase detector)
    /* COMMENTING OUT BECAUSE I_P == 0 at loopcount== */
    carrError = atan(Q_P / I_P) / (2.0 * pi);

    // Implement carrier loop filter and generate NCO command
    carrNco = oldCarrNco + (tau2carr / tau1carr) * (carrError - oldCarrError) +
              carrError * (PDIcarr / tau1carr);
    oldCarrNco = carrNco;
    oldCarrError = carrError;

    // Modify carrier freq based on NCO command
    carrFreq = carrFreqBasis + carrNco;

    // Find DLL error and update code NCO -------------------------------------
    codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) /
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

    // Implement code loop filter and generate NCO command
    codeNco = oldCodeNco + (tau2code / tau1code) * (codeError - oldCodeError) +
              codeError * (PDIcode / tau1code);
    oldCodeNco = codeNco;
    oldCodeError = codeError;

    // Modify code freq based on NCO command
    codeFreq = codeFreqBasis - codeNco;
    absoluteSample = ftell(fpdata) / dataAdaptCoeff - remCodePhase;

    // Store values in output arrays         vvv Corresponding variable in
    // Matlab (trackResults)
    carrFreq_output[loopcount] = carrFreq; // codeFreq
    codeFreq_output[loopcount] = codeFreq;
    absoluteSample_output[loopcount] = absoluteSample; // absoluteSample
    codeError_output[loopcount] = codeError;           // dllDiscr
    codeNco_output[loopcount] = codeNco;               // dllDiscrFilt
    carrError_output[loopcount] = carrError;           // pllDiscr
    carrNco_output[loopcount] = carrNco;               // pllDiscrFilt
    I_E_output[loopcount] = I_E;                       // I_E
    I_P_output[loopcount] = I_P;                       // I_P
    I_L_output[loopcount] = I_L;                       // I_L
    Q_E_output[loopcount] = Q_E;                       // Q_E
    Q_P_output[loopcount] = Q_P;                       // Q_P
    Q_L_output[loopcount] = Q_L;                       // Q_L

  } // end for

  // Log tracking results to file
  printf("  [Logging data into the "
         "'plot/data_avx512_fl32_add_mul_nom_lut_code' directory]\n");
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/codeNco_output.bin",
      codeNco_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/codeError_output.bin",
      codeError_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/carrNco_output.bin",
      carrNco_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/carrError_output.bin",
      carrError_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/absoluteSample_output.bin",
      absoluteSample_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/carrFreq_output.bin",
      carrFreq_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/codeFreq_output.bin",
      codeFreq_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/I_E_output.bin",
      I_E_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/I_P_output.bin",
      I_P_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/I_L_output.bin",
      I_L_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/Q_E_output.bin",
      Q_E_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/Q_P_output.bin",
      Q_P_output);
  write_file_fl64(
      "../plot/data_avx512_fl32_add_mul_nom_lut_code/Q_L_output.bin",
      Q_L_output);

  // Clearing up memory
  printf("  [Cleaning up used memory]\n");
  fclose(fpdata);
  free(rawSignal);
  free(carrFreq_output);
  free(codeFreq_output);
  free(absoluteSample_output);
  free(codeError_output);
  free(codeNco_output);
  free(carrError_output);
  free(carrNco_output);
  free(I_E_output);
  free(I_P_output);
  free(I_L_output);
  free(Q_E_output);
  free(Q_P_output);
  free(Q_L_output);
  free(VSMIndex);
  free(VSMValue);

  // Exit the program
  printf("*** Job Completed Succesfully! ***\n\n");
  return EXIT_SUCCESS;
}
