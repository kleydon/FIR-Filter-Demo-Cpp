//FFT.cpp
//
//Wrapper around Ooura's FFTsg implementation; see OouraFFTsg.cpp

#include "FFT.hpp"

#include "Logger.hpp"
#include "OouraFFTsg.hpp"
#include "Utility.hpp"
#include "VectorMath.hpp"

#include <cassert>

//#include "fftsg.hpp"



FFT::FFT(const uif32 length) {
    
    //Check arguments
    assert(isPowerOf2(length));
    this->length = length;
    
    //Oora-specific setup
    uif32 sqrtLength = floorf(sqrtf(length));
    ip = new int[sqrtLength + 2];
    w = new float[length * 5 / 4];
    ip[0] = 0; //Necessary before first transform
}



FFT::~FFT() {
    //Delete any dynamic allocations
    safeDeleteArray(ip);
    safeDeleteArray(w);
}



void FFT::forward(float data[]) const {
    rdft(length,
         +1,
         data,
         ip,
         w);
}



void FFT::backward(float data[]) const {
    
    rdft(length,
         -1,
         data,
         ip,
         w);
    
    //Normalize
    float scaleFactor = 2.0f / length;
    vm_mulvs(data,
             &scaleFactor,
             data,
             length);
}



void FFT::complexForward(float data[]) const {
    rdft(length,
         +1,
         data,
         ip,
         w);
}



void FFT::complexBackward(float data[]) const {
    
    cdft(length,
         -1,
         data,
         ip,
         w);
    
    //Normalize
    float scaleFactor = 1.0f / length;
    vm_mulvs(data,
             &scaleFactor,
             data,
             length);
}



uif32 FFT::getLength() const {
    return length;
}



uif32 FFT::calcComplexLength() const {
    return (length / 2) + 1;
}



void FFT::test() {
    
    const float duration = 1.0f; //sec
    const float frequency = 10.0f;
    const uif32 numSamples = 256;
    const uif32 fftLength = numSamples;
    
    FFT* fft256 = new FFT(fftLength);

    float data[numSamples];
    
    printf("\nFFT Test\n");
    printf("\t Samples: %u\n", numSamples);
    printf("\t Duration: %0.2fsec\n", duration);
    printf("\t Cosine Frequency: %0.2fHz\n", frequency);
    printf("\t FFT Size: %u\n", numSamples);
    printf("\t Bin Width: %0.2fHz\n", (numSamples/duration) / fftLength);
    printf("\n");
    
    //Calculate input waveform
    for (uif32 i = 0; i < numSamples; ++i) {
        float t = i * (duration / numSamples);
        //printf("t:%f\n", t[i]);
        data[i] = cos(2.0f * M_PI * frequency * t);
    }
    
    printf("Input waveform:\n");
    printf("\t");
    for (uif32 i = 0; i < numSamples; ++i) {
        printf("%f ", data[i]);
    }
    printf("\n\n");
    
    //Forward
    fft256->forward(data);
    
    printf("Forward Transform Output:\n");
    printf("\t");
    for (uif32 i = 0; i < numSamples/2; ++i) {
        printf("%f ", data[i*2]);
    }
    printf("\n\n");
    
    fft256->backward(data);
    
    printf("Backward FFT Output:\n");
    printf("\t");
    for (uif32 i = 0; i < numSamples; ++i) {
        printf("%f ", data[i]);
    }
    printf("\n\n");
}


