//VectorMathApple.cpp

#include "VectorMath.hpp"

#include "TypeAbbreviations.hpp"

#include <Accelerate/Accelerate.h>


#ifndef __APPLE__
    #error VectorMathApple included in a non-Apple build.
#endif



// Clear / Fill / Cast



void vm_clrv(float data[],
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vclr(data, 1,
              length);
#else
    vm_fallback_clrv(data,
                     length);
#endif
}



void vm_fillv(const float* value,
              float data[],
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vfill(value,
               data, 1,
               length);
#else
    vm_fallback_fillv(value,
                      data,
                      length);
#endif
}



void vm_hopcopy(const float source[],
                float dest[],
                const uif32 hopLength,
                const uif32 sourceLength) {
        
    vm_fallback_hopcopy(source,
                        dest,
                        hopLength,
                        sourceLength);
}



void vm_si16v_to_f32v(const si16 source[],
                      float dest[],
                      const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vflt16(source, 1,
                dest, 1,
                length);
#else
    vm_fallback_si16v_to_f32v(source,
                             dest,
                             length);
#endif
}



// Min / Max / Mean / Abs / Sum / MaxMag



void vm_minv(const float data[],
             float* result,
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_minv(data, 1,
              result,
              length);
#else
    vm_fallback_minv(data,
                     result,
                     length);
#endif
}



void vm_minvi(const float data[],
              float* result,
              uif32* resultIndex,
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_minvi(data, 1,
               result,
               (vDSP_Length *)resultIndex,
               length);
#else
    vm_fallback_minvi(data,
                      result,
                      resultIndex,
                      length);
#endif
}



void vm_maxv(const float data[],
             float* result,
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_maxv(data, 1,
              result,
              length);
#else
    vm_fallback_maxv(data,
                     result,
                     length);
#endif
}



void vm_maxvi(const float data[],
              float* result,
              uif32* resultIndex,
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_maxvi(data, 1,
               result,
               (vDSP_Length *)resultIndex,
               length);
#else
    vm_fallback_maxvi(data,
                      result,
                      resultIndex,
                      length);
#endif
}



void vm_meanv(const float data[],
              float* result,
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_meanv(data, 1,
               result,
               length);
#else
    vm_fallback_meanv(data,
                      result,
                      length);
#endif
}



void vm_absv(const float source[],
             float dest[],
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vabs(source, 1,
              dest, 1,
              length);
#else
    vm_fallback_absv(source,
                     dest,
                     length);
#endif
}



void vm_sumv(const float data[],
             float* result,
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_sve(data, 1,
             result,
             length);
#else
    vm_fallback_sumv(data,
                     result,
                     length);
#endif
}



void vm_maxmgvi(const float data[],
                float* result,
                uif32* resultIndex,
                const uif32 length) {
#ifdef VEC_MATH_HW_ACCELERATION_ENABLED
    vDSP_maxmgvi(data, 1,
                 result,
                 (vDSP_Length*) resultIndex,
                 (vDSP_Length) length);
#else
    vm_fallback_maxmgvi(data,
                        result,
                        resultIndex,
                        length);
#endif
}



// Add / Subtract / Multiply / Divide



void vm_addv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vadd(dataA, 1,
              dataB, 1,
              dataC, 1,
              length);
#else
    vm_fallback_addv(dataA,
                     dataB,
                     dataC,
                     length);
#endif
}



void vm_addvs(const float dataA[],
              const float* B,
              float dataC[],
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vsadd(dataA, 1,
               B,
               dataC, 1,
               length);
#else
    vm_fallback_addvs(dataA,
                      B,
                      dataC,
                      length);
#endif
}



void vm_subv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vsub(dataB, 1, //NOTE: Swapping order here to make more intuitive
              dataA, 1,
              dataC, 1,
              length);
#else
    vm_fallback_subv(dataA,
                     dataB,
                     dataC,
                     length);
#endif
}



void vm_mulv(const float dataA[],
             const float dataB[],
             float dataC[],
             const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vmul(dataA, 1,
              dataB, 1,
              dataC, 1,
              length);
#else
    vm_fallback_mulv(dataA,
                     dataB,
                     dataC,
                     length);
#endif
}



void vm_mulcv(const float dataA[],
              const float dataB[],
              float dataC[],
              const uif32 complexLength) {
    
    //Not using vDSP here, since it provides no means for multiplying
    //interleaved (DSPComplex) format directly, and it seems intensive
    //to convert to everything to DSPSplitComplex format and back again.
    
    vm_fallback_mulcv(dataA,
                      dataB,
                      dataC,
                      complexLength);
}



void vm_cvmags(const float dataA[],
               float dataC[],
               const uif32 complexLength) {
    
    //Not using vDSP here, since it provides no means for multiplying
    //interleaved (DSPComplex) format directly, and it seems intensive
    //to convert to everything to DSPSplitComplex format and back again.
    
    vm_fallback_cvmags(dataA,
                       dataC,
                       complexLength);
}



void vm_mulvs(const float dataA[],
              const float* B,
              float dataC[],
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_vsmul(dataA, 1,
               B,
               dataC, 1,
               length);
#else
    vm_fallback_mulvs(dataA,
                      B,
                      dataC,
                      length);
#endif
}


//Correlation



//See vm_fallback_acf1()
void vm_acf1(const float inputFrame[],
             float processingFrame[],  //Space >= inputFrameLength
             float resultFrame[],  //Space >= inputFrameLength, result = inputFrameLength
             const uif32 inputFrameLength, //Power of 2
             const FFT* fftContextOfInputFrameLength,
             const bool normalized) {
    
    vm_fallback_acf1(inputFrame,
                     processingFrame,  //Space >= inputFrameLength
                     resultFrame,  //Space >= inputFrameLength
                     inputFrameLength, //Power of 2
                     fftContextOfInputFrameLength,
                     normalized);
}

//See vm_fallback_acf2()
void vm_acf2(const float inputFrameSamples[],
            float resultFrame[],
            const uif32 inputFrameLength,
            const FFT* fftContextOfInputFrameLengthX2,
            const bool unbiased,
            const bool normalized) {
    
    vm_fallback_acf2(inputFrameSamples,
                     resultFrame,
                     inputFrameLength,
                     fftContextOfInputFrameLengthX2,
                     unbiased,
                     normalized);
}



// Window functions: blackman / hamming / hanning / gaussian



void vm_blackmanv(float data[],
                  const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_blkman_window(data,
                       length,
                       0); //0 for full-sized window
#else
    vm_fallback_blackmanv(data,
                          length);
#endif
}



void vm_hammingv(float data[],
                 const uif32 length) {
#ifdef VEC_MATH_HW_ACCELERATION_ENABLED
    vDSP_hamm_window(data,
                     length,
                     0); //0 for full-sized window
#else
    vm_fallback_hammingv(data,
                         length);
#endif
}



void vm_hanningv(float data[],
              const uif32 length) {
#ifdef VECTOR_MATH_HW_ACCELERATION_ENABLED
    vDSP_hann_window(data,
                     length,
                     0); //0 for full-sized window
#else
    vm_fallback_hanningv(data,
                         length);
#endif
}



void vm_gaussianv(const float sigma,
                  float coefficients[],
                  const uif32 length) {

    vm_fallback_gaussianv(sigma,
                          coefficients,
                          length);
}



