//VectorMath.hpp

//Miscellaneous convenience functions.

#ifndef __VECTOR_MATH_HPP__
#define __VECTOR_MATH_HPP__

#include "TypeAbbreviations.hpp"
#include "FFT.hpp"


//#define VECTOR_MATH_HW_ACCELERATION_ENABLED



// Clear / Fill / Cast



//vm_clrv()
//Set elements of single-precision float array to zero
//++++++++++++++++++++++++++
void vm_clrv(float data[],
             const uif32 length);

void vm_fallback_clrv(float data[],
                      const uif32 length);
//--------------------------



//vm_fillv()
//Set elements of single-precision float array to given value
//++++++++++++++++++++++++++
void vm_fillv(const float* value,
              float data[],
              const uif32 length);

void vm_fallback_fillv(const float* value,
                       float data[],
                       const uif32 length);
//--------------------------



//vm_copy_nth()
//Set elements of single-precision float array to given value
//++++++++++++++++++++++++++
void vm_hopcopy(const float source[],
                float dest[],
                const uif32 hopLength,
                const uif32 sourceLength);

void vm_fallback_hopcopy(const float source[],
                         float dest[],
                         const uif32 hopLength,
                         const uif32 sourceLength);
//--------------------------



//vm_si16v_to_f32v()
//Converts an array of signed 16-bit integers
//to single-precision floating-point values
//++++++++++++++++++++++++++
void vm_si16v_to_f32v(const si16 source[],
                      float dest[],
                      const uif32 length);

void vm_fallback_si16v_to_f32v(const short source[],
                               float dest[],
                               const uif32 length);
//--------------------------



// Min / Max / Mean / Abs / Sum / MaxMag



//vm_minv()
//Max value for single-precision float array
//++++++++++++++++++++++++++
void vm_minv(const float data[],
             float* result,
             const uif32 length);

void vm_fallback_minv(const float data[],
                      float* result,
                      const uif32 length);
//--------------------------



//vm_minvi()
//Max value and its index, over a single-precision float array
//++++++++++++++++++++++++++
void vm_minvi(const float data[],
              float* result,
              uif32* resultIndex,
              const uif32 length);

void vm_fallback_minvi(const float data[],
                       float* result,
                       uif32* resultIndex,
                       uif32 length);
//--------------------------



//vm_maxv()
//Max value for single-precision float array
//++++++++++++++++++++++++++
void vm_maxv(const float data[],
             float* result,
             const uif32 length);

void vm_fallback_maxv(const float data[],
                      float* result,
                      const uif32 length);
//--------------------------



//vm_maxvi()
//Max value and its index, over a single-precision float array
//++++++++++++++++++++++++++
void vm_maxvi(const float data[],
              float* result,
              uif32* resultIndex,
              const uif32 length);

void vm_fallback_maxvi(const float data[],
                       float* result,
                       uif32* resultIndex,
                       uif32 length);
//--------------------------



//vm_meanv()
//Mean value for single-precision float array
//++++++++++++++++++++++++++
void vm_meanv(const float data[],
              float* result,
              const uif32 length);

void vm_fallback_meanv(const float data[],
                       float* result,
                       const uif32 length);
//--------------------------



//vm_absv()
//Element-wise absolute value for single-precision float array
//++++++++++++++++++++++++++
void vm_absv(const float source[],
             float dest[],
             const uif32 length);

void vm_fallback_absv(const float source[],
                      float dest[],
                      const uif32 length);
//--------------------------



//vm_sumv()
//Sum over array A
//++++++++++++++++++++++++++
void vm_sumv(const float data[],
             float* result,
             const uif32 length);

void vm_fallback_sumv(const float data[],
                      float* result,
                      const uif32 length);
//--------------------------



//vm_maxmgvi()
//Max magnitude and index over array A
//++++++++++++++++++++++++++
void vm_maxmgvi(const float data[],
                float* result,
                uif32* resultIndex,
                const uif32 length);

void vm_fallback_maxmgvi(const float data[],
                         float* result,
                         uif32* resultIndex,
                         const uif32 length);
//--------------------------



// Add / Subtract / Multiply / Divide



//vm_addv()
//Element-wise A + B
//++++++++++++++++++++++++++
void vm_addv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length);

void vm_fallback_addv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length);
//--------------------------


//vm_addvs()
//Add scalar value B to all elements of a A
//++++++++++++++++++++++++++
void vm_addvs(const float dataA[],
              const float* B,
              float* dataC,
              const uif32 length);

void vm_fallback_addvs(const float dataA[],
                       const float* B,
                       float* dataC,
                       const uif32 length);
//--------------------------



//vm_subv()
//Element-wise A - B
//++++++++++++++++++++++++++
void vm_subv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length);

void vm_fallback_subv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length);
//--------------------------



//vm_mulv()
//Multiply two single-precision float arrays, element by element
//++++++++++++++++++++++++++
void vm_mulv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length);

void vm_fallback_mulv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length);
//--------------------------



//vm_mulcv()
//Multiply two interleaved complex single-precision float arrays
//of the format [R0,I0, R1,I1, ...] element-pair by element-pair,
//where the complex length is 1/2 the number of floats in each
//vector. //Can be in place with respect to A or B.
//++++++++++++++++++++++++++
void vm_mulcv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 complexLength);

void vm_fallback_mulcv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 complexLength);
//--------------------------



//vm_cvmags()
//Calculate the absolute square magnitude of a complex single-precision
//float arrays of the format [R0,I0, R1,I1, ...], yielding a real result
//in complex format.
//* Can be in-place with respect to input A.
//* Result is real-only, but returned in complex format
//++++++++++++++++++++++++++
void vm_cvmags(const float dataA[],
              float dataC[],
              const uif32 complexLength);

void vm_fallback_cvmags(const float dataA[],
                        float dataC[],
                        const uif32 complexLength);
//--------------------------



//vm_mulvs()
//Multiply all elements of a single-precision float array by a
//scalar value
//++++++++++++++++++++++++++
void vm_mulvs(const float dataA[],
              const float* B,
              float dataC[],
              const uif32 length);

void vm_fallback_mulvs(const float dataA[],
                       const float* B,
                       float dataC[],
                       const uif32 length);
//--------------------------



//Corrrelation



//vm_acf1() - used by Yin
//
//Calculates Type 1 autocorrelation function (ACF):
//
//  for (tau = 0, tau < tauRange; tau++) {
//      sum = 0.0f;
//      for (j = 0, j < W, j++) {
//          sum += signal[j] * signal[j + tau];
//      }
//      acf[tau] = runSum;
//  }
//
//  If the input signal is NOT zero-padded outside (after) integration window W,
//  the function will NOT taper with tau. This is generally desireable for pitch
//  detection, as high-lag peaks are detectable; the Type 1 ACF is used by Yin.
//
//  If the input signal IS zero-padded outside (after) integration window W,
//  or doesn't exist after W, the function WILL taper with tau; in this case,
//  Type 1 autocorrelation gives the same (tapered) result as Type 2
//  autocorrelation. The tapering is also known as "biasing". It can be
//  compensated for by multiplying ACF terms by 1/(W-j), but this can highlight
//  numerical errors at high lag values due to the integration covering
//  smaller and smaller portions of the input signal as lag values increase.
//
//  Normalization by division ensures that the maximum (first) value is
//  always 1, and that the function ranges from -1 to 1.
//
//  For some discussion of Type 1 vs Type 2, see:
//    * "A Smarter Way to Find Pitch", by Philip McLeod
//      http://www.cs.otago.ac.nz/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf
//    * "YIN, a fundamental frequency estimator for speech and music" by
//      Alain de Cheveigne, JASA:
//      http://audition.ens.fr/adc/pdf/2002_JASA_YIN.pdf
//
//++++++++++++++++++++++++++
void vm_acf1(const float inputFrame[],
             float processingFrame[], //Space >= inputFrameLength
             float resultFrame[], //Space >= inputFrameLength, result = inputFrameLengthDiv2
             const uif32 inputFrameLength, //Power of 2
             const FFT* fftContextOfInputFrameLength,
             const bool normalized);

void vm_fallback_acf1(const float inputFrame[],
                      float processingFrame[], //Space >= inputFrameLength
                      float resultFrame[], //Space >= inputFrameLength, result = inputFrameLengthDiv2
                      const uif32 inputFrameLength, //Power of 2
                      const FFT* fftContextOfInputFrameLength,
                      bool normalized);
//--------------------------



//vm_acf2() - used by MPM
//
//Calculates Type 2 autocorrelation function (ACF):
//
//  for (tau = 0, tau < tauRange; tau++) {
//      sum = 0.0f;
//      for (j = 0, j < W - tau, j++) {  //<--- "- tau" not present for Type 1.
//          sum += signal[j] * signal[j + tau];
//      }
//      acf[tau] = runSum;
//  }
//
//  This autocorrelation tapers with increasing lag tau.
//
//  The tapering is also known as "biasing". It can be compensated for by
//  multiplying ACF terms by 1/(W-j), but this can highlight numerical errors
//  at high lag values due to the integration covering smaller and smaller
//  portions of the input signal as lag values increase.
//
//  Normalization can ensures that the maximum (first) value is always 1,
//  and that the function ranges from -1 to 1. (**MAY NOT BE TRUE IF
//  UNBIASED, DEPENDING ON NUMERICAL ERRORS AT HIGH LAG VALUES**...)
//
//  For some discussion of Type 1 vs Type 2, see:
//    * "A Smarter Way to Find Pitch", by Philip McLeod
//      http://www.cs.otago.ac.nz/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf
//    * "YIN, a fundamental frequency estimator for speech and music" by
//      Alain de Cheveigne, JASA:
//      http://audition.ens.fr/adc/pdf/2002_JASA_YIN.pdf
//
//++++++++++++++++++++++++++
void vm_acf2(const float inputFrame[],
             float resultFrame[], //Space >= 2*inputFrameLength; end result = inputFrameLength
             const uif32 inputFrameLength, //Power of 2
             const FFT* fftContextOfInputFrameLengthX2,
             const bool unbiased,
             const bool normalized);

void vm_fallback_acf2(const float inputFrameSamples[],
                      float resultFrame[], //Space >= 2*inputFrameLength; end result = inputFrameLength
                      const uif32 inputFrameLength, //Power of 2
                      const FFT* fftContextOfInputFrameLengthX2,
                      const bool unbiased,
                      bool normalized);
//--------------------------



// Window functions: blackman / hamming / hanning / gaussian



//vm_blackmanv()
//Calculate the coefficients of a blackman window
//++++++++++++++++++++++++++
void vm_blackmanv(float coefficients[],
                  const uif32 length);

void vm_fallback_blackmanv(float coefficients[],
                           const uif32 length);
//--------------------------



//vm_hammingv()
//Calculate the coefficients of a hamming window
//++++++++++++++++++++++++++
void vm_hammingv(float coefficients[],
                 const uif32 length);

void vm_fallback_hammingv(float coefficients[],
                          const uif32 length);
//--------------------------



//vm_hanningv()
//Calculate the coefficients of a hanning window
//++++++++++++++++++++++++++
void vm_hanningv(float coefficients[],
                 const uif32 length);

void vm_fallback_hanningv(float coefficients[],
                          const uif32 length);
//--------------------------



//vm_gaussianv()
//Calculate the coefficients of a gaussian window
//++++++++++++++++++++++++++
void vm_gaussianv(const float sigma,
                  float coefficients[],
                  const uif32 length);

void vm_fallback_gaussianv(const float sigma,
                           float coefficients[],
                           const uif32 length);
//--------------------------



// Print / Test

void printVector(const float data[],
                 const uif32 length);


//vm_test()
void vm_test();



#endif //__VECTOR_MATH_HPP__
