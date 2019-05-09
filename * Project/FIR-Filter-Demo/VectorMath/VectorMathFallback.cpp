//VectorMath.cpp

#include "VectorMath.hpp"

#include "Logger.hpp"
#include "TypeAbbreviations.hpp"

#include <cmath>
#include <cassert>



// Clear / Fill / Cast



void vm_fallback_clrv(float data[],
                      const uif32 length) {
    //memset() works here, since
    //IEEE standard representation for 0.0 is all zero bits.
    //See:
    //http://stackoverflow.com/a/877401
    std::memset(data, 0, sizeof(float) * length);
}



void vm_fallback_fillv(const float* value,
                       float data[],
                       const uif32 length) {
    
    std::fill_n(data, length, *value);
}



void vm_fallback_hopcopy(const float source[],
                         float dest[],
                         const uif32 hopLength,
                         const uif32 sourceLength) {
    uif32 j = 0;
    for (uif32 i = 0; i < sourceLength; i += hopLength) {
        dest[j++] = source[i];
    }
}



void vm_fallback_si16v_to_f32v(const si16 source[],
                               float dest[],
                               const uif32 length) {
    
    for (uif32 i = 0; i < length; ++i) {
        dest[i] = (float) source[i];
    }
}



// Min / Max / Mean / Abs / Sum / MaxMag



void vm_fallback_minv(const float data[],
                      float* result,
                      const uif32 length) {
    
    float min = std::numeric_limits<float>::max();
    
    for (uif32 i = 0; i < length; ++i) {
        if (min > data[i]) {
            min = data[i];
        }
    }
    
    *result = min;
}



void vm_fallback_minvi(const float data[],
                       float* result,
                       uif32* resultIndex,
                       const uif32 length) {
    
    assert(length > 0);
    
    float min = std::numeric_limits<float>::max();
    uif32 index = 0; //Should be undefined
    
    for (uif32 i = 0; i < length; ++i) {
        if (min > data[i]) {
            min = data[i];
            index = i;
        }
    }
    
    *result = min;
    *resultIndex = index;
}



void vm_fallback_maxv(const float data[],
                      float* result,
                      const uif32 length) {
    
    float max = std::numeric_limits<float>::lowest();
    
    for (uif32 i = 0; i < length; ++i) {
        if (max < data[i]) {
            max = data[i];
        }
    }
    
    *result = max;
}



void vm_fallback_maxvi(const float data[],
                       float* result,
                       uif32* resultIndex,
                       const uif32 length) {
    
    assert(length > 0);
    
    float max = std::numeric_limits<float>::lowest(); //use lowest or min?
    uif32 index = 0; //Should be undefined
    
    for (uif32 i = 0; i < length; ++i) {
        if (max < data[i]) {
            max = data[i];
            index = i;
        }
    }
    
    *result = max;
    *resultIndex = index;
}



void vm_fallback_meanv(const float data[],
                       float* result,
                       const uif32 length) {
    
    float sum = 0.0f;
    
    for (uif32 i = 0; i < length; ++i) {
        sum += data[i];
    }
    
    *result = sum /= length;
}



void vm_fallback_absv(const float source[],
                      float dest[],
                      const uif32 length) {

    for (uif32 i = 0; i < length; ++i) {
        dest[i] = std::fabsf(source[i]);
    }
}



void vm_fallback_sumv(const float data[],
                      float* result,
                      const uif32 length) {
    
    float sum = 0.0f;
    
    for (uif32 i = 0; i < length; i++) {
        sum += data[i];
    }
    
    *result = sum;
}



void vm_fallback_maxmgvi(const float data[],
                         float* result,
                         uif32* resultIndex,
                         const uif32 length) {
    float maxMag = 0.0f;
    uif32 index = 0; //Should be undefined
    
    for (uif32 i = 0; i < length; ++i) {
        float mag = std::fabsf(data[i]);
        if (maxMag < mag) {
            maxMag = mag;
            index = i;
        }
    }
    
    *result = maxMag;
    *resultIndex = index;
}



// Add / Subtract / Multiply / Divide



void vm_fallback_addv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length) {
    for (uif32 i = 0; i < length; ++i) {
        dataC[i] = dataA[i] + dataB[i];
    }
}



void vm_fallback_addvs(const float dataA[],
                       const float* B,
                       float dataC[],
                       const uif32 length) {
    float b = *B;
    for (uif32 i = 0; i < length; ++i) {
        dataC[i] = dataA[i] + b;
    }
}



void vm_fallback_subv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length) {
    
    for (uif32 i = 0; i < length; ++i) {
        dataC[i] = dataA[i] - dataB[i];
    }
}



void vm_fallback_mulv(const float dataA[],
                      const float dataB[],
                      float dataC[],
                      const uif32 length) {
    
    for (uif32 i = 0; i < length; ++i) {
        dataC[i] = dataA[i] * dataB[i];
    }
}



void vm_fallback_mulcv(const float dataA[],
                       const float dataB[],
                       float dataC[],
                       const uif32 complexLength) {

    for (uif32 i = 0; i < complexLength; ++i) {
        
        uif32 realIndex = i * 2;
        uif32 imagIndex = realIndex + 1;
        
        //Re[(a+bi)(c+di)] = (ac - bd)
        //Im[(a+bi)(c+di)] = (a + b)(c + d) - ac - bd
        
        float ac = dataA[realIndex] * dataB[realIndex];
        float bd = dataA[imagIndex] * dataB[imagIndex];
        
        //Imaginary
        dataC[imagIndex] =
        ((dataA[realIndex] + dataA[imagIndex]) *
         (dataB[realIndex] + dataB[imagIndex])) - ac - bd;
        
        //Real
        dataC[realIndex] = ac - bd;
    }
}



void vm_fallback_cvmags(const float dataA[],
                        float dataC[],
                        const uif32 complexLength) {

    //Complex square magnitude = |x+iy|^2 = x^2 + y^2
    
    for (uif32 i = 0; i < complexLength; ++i) {
        
        uif32 realIndex = i * 2;
        uif32 imagIndex = realIndex + 1;
        
        float re = dataA[realIndex];
        float im = dataA[imagIndex];
        
        dataC[realIndex] = re * re  +  im * im;
        dataC[imagIndex] = 0.0f;
    }
}



void vm_fallback_mulvs(const float dataA[],
                       const float* B,
                       float dataC[],
                       const uif32 length) {
    float b = *B;
    for (uif32 i = 0; i < length; ++i) {
        dataC[i] = dataA[i] * b;
    }
}



// Correlation

void vm_fallback_acf1(const float inputFrame[],
                      float processingFrame[], //Space >= inputFrameLength
                      float resultFrame[], //Space >= inputFrameLength, result = inputFrameLength
                      const uif32 inputFrameLength, //Power of 2
                      const FFT* fftContextOfInputFrameLength,
                      bool normalized) {
    
    //See notes/references on autocorrelation in:
    // * VectorMath.h/vm_fallback_acf2()
    // * F0Detector.cpp
    
    uif32 inputFrameLengthDiv2 = inputFrameLength / 2;
    
    //Fill resultFrame with input samples
    memcpy(resultFrame,
           inputFrame,
           sizeof(float) * inputFrameLength);
    
    //Forward FFT of audio data
    fftContextOfInputFrameLength->forward(resultFrame);
    
    //Fill processingFrame with first half of input samples,
    //followed by zero padding. (Reversal in the time
    //domain equates to complex conjugation in the frequency
    //domain).
    for (uif32 i = 0; i < inputFrameLengthDiv2; ++i) {
        processingFrame[i] = inputFrame[inputFrameLengthDiv2-1-i];
    }
    vm_clrv(processingFrame + inputFrameLengthDiv2,
            inputFrameLengthDiv2);
    
    //Forward FFT of 1/2 audio data, reversed
    fftContextOfInputFrameLength->forward(processingFrame);
    
    //Set 0 frequency bins (which are real-only) to 0
    resultFrame[0] = 0.0f;
    processingFrame[0] = 0.0f;
    
    //Stash (real) nyquist components, oddly stored by FFT,
    //in preparation for obtaining square magnitudes
    float nyquistComponentA = resultFrame[1];
    resultFrame[1] = 0.0f;
    float nyquistComponentB = processingFrame[1];
    processingFrame[1] = 0.0f;
    
    //Complex multiply
    //(Complex multiplication in the frequency domain equates to
    //real convolution in the time domain.)
    vm_mulcv(processingFrame,
            resultFrame,
            resultFrame, //result
            inputFrameLengthDiv2);
    
    //Restore nyquist components, multiplied, to result
    resultFrame[1] = nyquistComponentA * nyquistComponentB;
    
    //Perform Inverse FFT
    fftContextOfInputFrameLength->backward(resultFrame);
    
    //If normalizing, divide all elements by the first
    if (normalized) {
        float divisor = 1.0f / resultFrame[0];
        vm_mulvs(resultFrame,
                 &divisor,
                 resultFrame,
                 inputFrameLength);
    }
}



void vm_fallback_acf2(const float inputFrame[],
                      float resultFrame[], //Space >= 2*inputFrameLength
                      const uif32 inputFrameLength, //Power of 2
                      const FFT* fftContextOfInputFrameLengthX2, //2XinputFrameLength
                      bool unbiased,
                      bool normalized) {
    
    //See notes/references on autocorrelation in:
    // * VectorMath.h/vm_fallback_acf2()
    // * F0Detector.cpp
    
    //Fill resultFrame with input samples and zero-pading, of equal length
    vm_clrv(resultFrame + inputFrameLength, inputFrameLength);
    memcpy(resultFrame,
           inputFrame,
           sizeof(float) * inputFrameLength);
    
    //Forward FFT
    fftContextOfInputFrameLengthX2->forward(resultFrame);
    
    //Set 0 frequency bin (which is real-only) to 0
    resultFrame[0] = 0.0f;
    
    //Stash (real) nyquist component, oddly stored by FFT,
    //in preparation for obtaining square magnitudes
    float nyquistComponent = resultFrame[1];
    resultFrame[1] = 0.0f;
    
    //Get square magnitudes.
    //Result is a real vector 1/2 zpadded length
    vm_cvmags(resultFrame, //input
              resultFrame, //output
              inputFrameLength); // 1/2 zpadded float vector lengths
    
    //Restore nyquist component (squared)
    resultFrame[1] = nyquistComponent * nyquistComponent;
    
    //Perform Inverse FFT
    fftContextOfInputFrameLengthX2->backward(resultFrame);
    
    //If unbiasing, divide by smaller and smaller amounts
    if (unbiased) {
        for (uif32 i = 0; i < inputFrameLength; ++i) {
            resultFrame[i] /= (inputFrameLength - i);
        }
    }
    
    //If normalizing, divide all elements by the first
    if (normalized) {
        float divisor = 1.0f / resultFrame[0];
        vm_mulvs(resultFrame,
                 &divisor,
                 resultFrame,
                 inputFrameLength);
    }
}



// Window functions: blackman / hamming / hanning / gaussian



void vm_fallback_blackmanv(float coefficients[],
                           const uif32 length) {
    
    for (uif32 i = 0; i < length; ++i) {
        coefficients[i] = 0.42 - (0.5f * cos(  2.0f * M_PI * i / length ) ) +
            (0.08 * cos( 4 * M_PI * i / length) );
    }
}



void vm_fallback_hammingv(float coefficients[],
                          const uif32 length) {
    for (uif32 i = 0; i < length; ++i) {
        coefficients[i] = 0.54f - (0.46f * cos( (2.0f * M_PI * i) / length ) );
    }
}



void vm_fallback_hanningv(float coefficients[],
                          const uif32 length) {

    //From matlab:
    //http://stackoverflow.com/a/24947868
    
    uif32 half, idx;
    
    if( length % 2 == 0) {
        
        half = length / 2;
        
        for (uif32 i = 0; i < half; ++i) {
            coefficients[i] = 0.5f * (1.0f - cos( 2.0f * M_PI * (i + 1) / (length + 1)));
        }
        
        idx = half - 1;
        for(uif32 i = half; i < length; ++i) {
            coefficients[i] = coefficients[idx];
            --idx;
        }
    }
    else {
        
        half = (length + 1) / 2;
        
        for(uif32 i = 0; i < half; ++i) {
            coefficients[i] = 0.5f * (1.0f - cos(2.0f * M_PI * (i + 1) / (length + 1)));
        }
        
        idx = half - 2;
        for(uif32 i = half; i < length; ++i) {
            coefficients[i] = coefficients[idx];
            --idx;
        }
    }
}



void vm_fallback_gaussianv(const float sigma,
                           float coefficients[],
                           const uif32 length) {
    
    //Based on: https://ccrma.stanford.edu/~jos/sasp/Matlab_Gaussian_Window.html
    
    float n[length];
    n[0] = -((float)length - 1.0f) / 2.0f;
    for (uif32 i = 1; i < length; ++i) {
        n[i] = n[i-1] + 1;
    }
    
    float x[length];
    for (uif32 i = 0; i < length; ++i) {
        x[i] = -n[i];
        x[i] *= n[i];
        x[i] /= 2.0f * sigma * sigma;
        coefficients[i] = std::expf(x[i]);
    }
}



// Print / Test


void printVector(const float data[],
                 const uif32 length) {
    for (uif32 i=0; i < length; i++) {
        iLog("[{}]: {:4.2f}", i, data[i]);
    }
    iLog(" ");
}





void vm_test() {
    
    float result;
    uif32 resultIndex;
    float value;
    
    // Clear / Fill / Cast
    
    iLog("Testing vm_clrv()...");
    
    float dataA[] = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    uif32 dataLength = sizeof(dataA) / sizeof(dataA[0]);
    float dataB[dataLength];
    float dataC[dataLength];
    
    iLog("dataA, length:{}", dataLength);
    printVector(dataA, dataLength);
    
    iLog("Called clear.");
    vm_clrv(dataA, dataLength);
    
    iLog("dataA:");
    printVector(dataA, dataLength);
    
    iLog("Testing vm_fillv()...");
    
    iLog("Set every element to 5.");
    value = 5;
    vm_fillv(&value, dataA, dataLength);
    
    iLog("dataA:");
    printVector(dataA, dataLength);

    iLog("Testing vm_si16v_to_f32v()...");
    
    si16 dataS16[dataLength];
    iLog("dataS16:");
    for (si16 i = 0; i < dataLength; ++i) {
        dataS16[i] = i - 2;
        iLog("[{}]: {}", i, dataS16[i]);
    }
    
    iLog("Converting si16 to float.");
    vm_si16v_to_f32v(dataS16,
                     dataA,
                     dataLength);
    
    iLog("dataA:");
    printVector(dataA, dataLength);
    
    iLog("Testing min/max value/index functions...");
    
    dataA[0] = -1.0f;
    dataA[1] = -2.0f;
    dataA[2] = -3.0f;
    dataA[3] = -4.0f;
    dataA[4] = -5.0f;
    
    iLog("dataA:");
    printVector(dataA, dataLength);
    
    vm_minv(dataA, &result, dataLength);
    iLog("Minimum: {}", result);
    
    vm_minvi(dataA, &result, &resultIndex, dataLength);
    iLog("Min: {} Index: {}", result, resultIndex);
    
    vm_maxv(dataA, &result, dataLength);
    iLog("Maximum: {}", result);
    
    vm_maxvi(dataA, &result, &resultIndex, dataLength);
    iLog("Max: {} Index: {}", result, resultIndex);
    
    vm_meanv(dataA, &result, dataLength);
    iLog("Mean: {}", result);
    
    vm_absv(dataA, dataB, dataLength);
    iLog("Abs:");
    printVector(dataB, dataLength);
    
    vm_sumv(dataA, &result, dataLength);
    iLog("Sum: {}", result);
    
    vm_maxmgvi(dataA, &result, &resultIndex, dataLength);
    iLog("MaxMag: {} Index: {}", result, resultIndex);

    dataA[0] = 5.0f;
    dataA[1] = 4.0f;
    dataA[2] = 3.0f;
    dataA[3] = 2.0f;
    dataA[4] = 1.0f;
    
    iLog("dataA:");
    printVector(dataA, dataLength);

    vm_minv(dataA, &result, dataLength);
    iLog("Minimum: {}", result);
    
    vm_minvi(dataA, &result, &resultIndex, dataLength);
    iLog("Min: {} Index: {}", result, resultIndex);
    
    vm_maxv(dataA, &result, dataLength);
    iLog("Maximum: {}", result);
    
    vm_maxvi(dataA, &result, &resultIndex, dataLength);
    iLog("Max: {} Index: {}", result, resultIndex);
    
    vm_meanv(dataA, &result, dataLength);
    iLog("Mean: {}", result);
    
    vm_absv(dataA, dataB, dataLength);
    iLog("Abs:");
    printVector(dataB, dataLength);
    
    vm_sumv(dataA, &result, dataLength);
    iLog("Sum: {}", result);
    
    vm_maxmgvi(dataA, &result, &resultIndex, dataLength);
    iLog("MaxMag: {} Index: {}", result, resultIndex);
    
    dataA[0] = 5.0f;
    dataA[1] = 4.0f;
    dataA[2] = 3.0f;
    dataA[3] = 4.0f;
    dataA[4] = 5.0f;
    
    iLog("dataA:");
    printVector(dataA, dataLength);
    
    vm_minv(dataA, &result, dataLength);
    iLog("Minimum: {}", result);
    
    vm_minvi(dataA, &result, &resultIndex, dataLength);
    iLog("Min: {} Index: {}", result, resultIndex);
    
    vm_maxv(dataA, &result, dataLength);
    iLog("Maximum: {}", result);
    
    vm_maxvi(dataA, &result, &resultIndex, dataLength);
    iLog("Max: {} Index: {}", result, resultIndex);
    
    vm_meanv(dataA, &result, dataLength);
    iLog("Mean: {}", result);
    
    vm_absv(dataA, dataB, dataLength);
    iLog("Abs:");
    printVector(dataB, dataLength);
    
    vm_sumv(dataA, &result, dataLength);
    iLog("Sum: {}", result);
    
    vm_maxmgvi(dataA, &result, &resultIndex, dataLength);
    iLog("MaxMag: {} Index: {}", result, resultIndex);

    
    iLog("Testing addition/subtraction/multiplication/division functions...");
    
    dataA[0] = 5.0f;
    dataA[1] = 4.0f;
    dataA[2] = 3.0f;
    dataA[3] = 4.0f;
    dataA[4] = 5.0f;
    
    dataB[0] = 1.0f;
    dataB[1] = 2.0f;
    dataB[2] = 3.0f;
    dataB[3] = 4.0f;
    dataB[4] = 10.0f;
    
    iLog("dataA:");
    printVector(dataA, dataLength);
    
    iLog("dataB:");
    printVector(dataB, dataLength);
    
    iLog("Adding dataA and dataB, storing in dataC.");
    vm_addv(dataA, dataB, dataC, dataLength);
    
    iLog("dataC:");
    printVector(dataC, dataLength);
    
    iLog("Adding 5 to dataA, storing in dataC.");
    value = 5;
    vm_addvs(dataA, &value, dataC, dataLength);
    
    iLog("dataC:");
    printVector(dataC, dataLength);
    
    iLog("Subtracting dataB from dataA, storing in dataC.");
    vm_subv(dataA, dataB, dataC, dataLength);
    
    iLog("dataC:");
    printVector(dataC, dataLength);
    
    iLog("Multiplying dataA by dataB, storing in dataC");
    
    vm_mulv(dataA, dataB, dataC, dataLength);
    
    iLog("dataC:");
    printVector(dataC, dataLength);
    
    iLog("Multiplying dataA by 10, storing in dataC");
    
    value = 10;
    vm_mulvs(dataA, &value, dataC, dataLength);
    
    iLog("dataC:");
    printVector(dataC, dataLength);

    const uif32 coefficientsLength = 50;
    float coefficients[coefficientsLength];
    
    vm_blackmanv(coefficients, coefficientsLength);
    iLog("Blackman coefficients:");
    for (uif32 i = 0; i < coefficientsLength; ++i) {
        printf("%f\n", coefficients[i]);
    }
    printf("\n");
    
    vm_hammingv(coefficients, coefficientsLength);
    iLog("Hamming coefficients:");
    for (uif32 i = 0; i < coefficientsLength; ++i) {
        printf("%f\n", coefficients[i]);
    }
    printf("\n");
    
    vm_hanningv(coefficients, coefficientsLength);
    iLog("Hanning coefficients:");
    for (uif32 i = 0; i < coefficientsLength; ++i) {
        printf("%f\n", coefficients[i]);
    }
    printf("\n");
    
    vm_gaussianv((float)coefficientsLength/6.0f, coefficients, coefficientsLength);
    iLog("Gaussian coefficients, r6:");
    for (uif32 i = 0; i < coefficientsLength; ++i) {
        printf("%f\n", coefficients[i]);
    }
    printf("\n");
     
    
    vm_gaussianv((float)coefficientsLength/8.0f, coefficients, coefficientsLength);
    iLog("Gaussian coefficients, r8:");
    for (uif32 i = 0; i < coefficientsLength; ++i) {
        printf("%f\n", coefficients[i]);
    }
    printf("\n");
     
    
}



