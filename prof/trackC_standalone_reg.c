/*!
 *  \file trackC_standalone_reg.c
 *  \brief      Simulates the tracking stage of a receiver using nominal C code
 *  \details    Profiles code when using:
 1. Carrier wave generation by means of DLUT method.
 2. Pseudorandom code generation by means of DLUT method.
 3. Down-conversion of the received signal by nominal multiplication.
 4. Multiplication of baseband signal with a local replica of
ranging code using with nominal C operations
 5. Accumulation to generate the correlation value with nominal C operations
 *  \author    Damian Miralles
 *  \author    Jake Johnson
 *  \date      Jan 23, 2018
 *  \pre       Make sure you have .bin files containing data and lookup tables.
 *  \code{.sh}
# Sample compilation script
$ gcc -I ../src/ trackC_standalone_reg.c -g
 -mavx2 -lm -o reg_standalone -O3
 *  \endcode
 */

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
  int i, loopcount, angle, blksize, pCode, eCode, lCode;
  int vsmCount, vsmInterval, dataAdaptCoeff;
  double remCodePhase, remCarrPhase, codePhaseStep;
  double earlyLateSpc, seekvalue, samplingFreq, trigarg, carrCos, carrSin,
      carrFreq;
  double I_E, Q_E, I_P, Q_P, I_L, Q_L, baseCode, mixedcarrSin = 0.0,
                                                 mixedcarrCos = 0.0;
  double carrNco, oldCarrNco, tau1carr, tau2carr, carrError, oldCarrError;
  double PDIcarr, codeNco, oldCodeNco, tau1code, tau2code, codeError,
      oldCodeError, PDIcode;
  double codeFreq, codeFreqBasis, carrFreqBasis, absoluteSample, codeLength;
  double pwr, pwrSum, pwrSqrSum, pwrMean, pwrVar, pwrAvgSqr, pwrAvg, noiseVar,
      CNo, accInt;
  char *rawSignal, *rawSignalI, *rawSignalQ;
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

  double caCode[1025]; //**
  getcaCodeFromFile("../data/caCode.bin", caCode);
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
  codeFreq = 1023002.79220779;
  codeFreqBasis = 1023002.79220779;
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

  // Allocate memory for the signal
  rawSignal = calloc(dataAdaptCoeff * blksize, sizeof(char));
  rawSignalI = calloc(blksize, sizeof(char));
  rawSignalQ = calloc(blksize, sizeof(char));

  // Open the file for reading the data and fseek if required
  fpdata = fopen(fileid, "rb");
  fseek(fpdata, dataAdaptCoeff * seekvalue, SEEK_SET);

  int sec_count = 0;
  printf("\n*** Running: trackC_standalone_reg ***\n");
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
    // blksize is supposed to be ceiling'd but for some reason compiler won't
    // find <math.h> header... >:(
    blksize = ceil((codeLength - remCodePhase) / codePhaseStep);

    i = fread(rawSignal, sizeof(char), dataAdaptCoeff * blksize, fpdata);

    // An error check should be added here to see if the required amount of data
    // can be read

    if (dataAdaptCoeff == 2) {
      for (i = 0; i < blksize; i++) {
        *(rawSignalI + i) = *(rawSignal + 2 * i);
        *(rawSignalQ + i) = *(rawSignal + 2 * i + 1);
      }
    }

    // Set up all the code phase tracking information
    // Define index into early,late and prompt codes

    for (i = 0; i < blksize; i++) {

      // Generate the carrier frequency to mix the signal to baseband
      trigarg = (2.0 * pi * carrFreq * (i / samplingFreq)) + remCarrPhase;
      angle = (int)(trigarg * 10000) % 62832;
      carrCos = (short)(round(8 * (gps_cos(angle))));
      carrSin = (short)(round(8 * (gps_sin(angle))));

      // First mix to baseband
      if (dataAdaptCoeff == 1) {
        mixedcarrSin = (carrSin * (*(rawSignal + i)));
        mixedcarrCos = (carrCos * (*(rawSignal + i)));
      } else if (dataAdaptCoeff == 2) {
        mixedcarrSin =
            (*(rawSignalI + i) * carrSin + *(rawSignalQ + i) * carrCos);
        mixedcarrCos =
            (*(rawSignalI + i) * carrCos - *(rawSignalQ + i) * carrSin);
      }

      baseCode = (i * codePhaseStep + remCodePhase);

      pCode = (int)(baseCode) < baseCode ? (baseCode + 1) : baseCode;
      eCode = (int)(baseCode - earlyLateSpc) < (baseCode - earlyLateSpc)
                  ? (baseCode - earlyLateSpc + 1)
                  : (baseCode - earlyLateSpc);
      lCode = (int)(baseCode + earlyLateSpc) < (baseCode + earlyLateSpc)
                  ? (baseCode + earlyLateSpc + 1)
                  : (baseCode + earlyLateSpc);

      // Now get early, late, and prompt values for each
      I_E += (*(caCode + eCode) * mixedcarrSin);
      Q_E += (*(caCode + eCode) * mixedcarrCos);
      I_P += (*(caCode + pCode) * mixedcarrSin);
      Q_P += (*(caCode + pCode) * mixedcarrCos);
      I_L += (*(caCode + lCode) * mixedcarrSin);
      Q_L += (*(caCode + lCode) * mixedcarrCos);
    }

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
  }

  // Log tracking results to file
  printf("  [Logging data into the 'plot/data_reg' directory]\n");
  write_file_fl64("../plot/data_reg/codeNco_output.bin", codeNco_output);
  write_file_fl64("../plot/data_reg/codeError_output.bin", codeError_output);
  write_file_fl64("../plot/data_reg/carrNco_output.bin", carrNco_output);
  write_file_fl64("../plot/data_reg/carrError_output.bin", carrError_output);
  write_file_fl64("../plot/data_reg/absoluteSample_output.bin",
                  absoluteSample_output);
  write_file_fl64("../plot/data_reg/carrFreq_output.bin", carrFreq_output);
  write_file_fl64("../plot/data_reg/codeFreq_output.bin", codeFreq_output);
  write_file_fl64("../plot/data_reg/I_E_output.bin", I_E_output);
  write_file_fl64("../plot/data_reg/I_P_output.bin", I_P_output);
  write_file_fl64("../plot/data_reg/I_L_output.bin", I_L_output);
  write_file_fl64("../plot/data_reg/Q_E_output.bin", Q_E_output);
  write_file_fl64("../plot/data_reg/Q_P_output.bin", Q_P_output);
  write_file_fl64("../plot/data_reg/Q_L_output.bin", Q_L_output);

  // Clearing up memory
  printf("  [Cleaning up used memory]\n");
  fclose(fpdata);
  free(rawSignal);
  free(rawSignalI);
  free(rawSignalQ);
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
