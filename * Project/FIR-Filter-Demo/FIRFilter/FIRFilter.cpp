//FIRFilter.cpp


#include "FIRFilter.hpp"

#include "Logger.hpp"
#include "TypeAbbreviations.hpp"
#include "Utility.hpp"
#include "VectorMath.hpp"

#include <Accelerate/Accelerate.h>
#include "assert.h"
#include <math.h>



const float FIRFilter::VALUE_NOT_AVAILABLE = 0.0f; //**SIMILAR NAME IN F0DETECTOR****



FIRFilter::FIRFilter() {
    //Begin with all pointers set to nullptr
    filterCoefficients = nullptr;
    filterPadded = nullptr;
    signalPadded = nullptr;
    signalPaddedXfilterPadded = nullptr;
    overlapBuffer = nullptr;
    tempResult = nullptr;
    tempResultB = nullptr;
    
    fft = nullptr;
};



FIRFilter::~FIRFilter() {
    //Delete any dynamic allocations
    safeDeleteArray(filterCoefficients);
    safeDeleteArray(filterPadded);
    safeDeleteArray(signalPadded);
    safeDeleteArray(signalPaddedXfilterPadded);
    safeDeleteArray(overlapBuffer);
    safeDeleteArray(tempResult);
    safeDeleteArray(tempResultB);
    
    safeDelete(fft);
};



bool FIRFilter::initializeAsLowPass(const float sampleRate,
                                    const uif32 frameLength,
                                    const float cutoffFreq,
                                    const float transBw) {
    return configure(sampleRate,
                     frameLength,
                     VALUE_NOT_AVAILABLE,
                     VALUE_NOT_AVAILABLE,
                     cutoffFreq,
                     transBw);
}



bool FIRFilter::initializeAsHighPass(const float sampleRate,
                                     const uif32 frameLength,
                                     const float cutoffFreq,
                                     const float transBw) {
    return configure(sampleRate,
                     frameLength,
                     cutoffFreq,
                     transBw,
                     VALUE_NOT_AVAILABLE,
                     VALUE_NOT_AVAILABLE);
}



bool FIRFilter::initializeAsBandPass(const float sampleRate,
                                     const uif32 frameLength,
                                     const float lowCutoffFreq,
                                     const float lowTransBw,
                                     const float highCutoffFreq,
                                     const float highTransBw) {
    return configure(sampleRate,
                     frameLength,
                     lowCutoffFreq,
                     lowTransBw,
                     highCutoffFreq,
                     highTransBw);
}



bool FIRFilter::configure(const float sampleRate,
                          const uif32 frameLength, //even power of two
                          const float lowCutoffFreq,
                          const float initialLowTransBw,
                          const float highCutoffFreq,
                          const float initialHighTransBw) {
    
    this->sampleRate = sampleRate;
    
    this->frameLength = frameLength;
    if (this->frameLength & 1) { //Length must be even
        this->frameLength++;
    }
    
    this->lowCutoffFreq = lowCutoffFreq;
    this->highCutoffFreq = highCutoffFreq;

    float provisionalLowTransBw = initialLowTransBw;
    float provisionalHighTransBw = initialHighTransBw;
    
    
    //Set:
    // * filterCoefficients
    // * filterLength - odd
    // * overlapLength - even
    //+++++++++
    //Start by setting filter coefficients and length from given paramters
    setFilterCoefficientsAndLength(sampleRate,
                                   lowCutoffFreq,
                                   provisionalLowTransBw,
                                   highCutoffFreq,
                                   provisionalHighTransBw);

    //Set overlapLength (overlap across filter calls) - even
    uif32 initialOverlapLength = filterLength - 1;
    overlapLength = initialOverlapLength;
    
    //If required, modify transition bandwidths and recalculate, until
    //overlapLength <= frameLength (to ensure overlap-and-add method
    //will be accurate).
    while (overlapLength > frameLength) {
        
        //Increase transition bandwidths while keeping ratio the same
        provisionalLowTransBw *= 1.2f;
        provisionalHighTransBw *= 1.2f;
            
        setFilterCoefficientsAndLength(sampleRate,
                                       lowCutoffFreq,
                                       provisionalLowTransBw,
                                       highCutoffFreq,
                                       provisionalHighTransBw);
        
        overlapLength = filterLength - 1;
    }
    
    this->lowTransBw = provisionalLowTransBw;
    this->highTransBw = provisionalHighTransBw;

    //If transition bandwidths had to be broadened...
    if (provisionalLowTransBw != initialLowTransBw ||
        provisionalHighTransBw != initialHighTransBw) {
        
        iLog("FIRFilter::initialize() - ***WARNING***");
        
        iLog("\t Initial transition bandwidths ({}, {}) forced to ({}, {})",
             initialLowTransBw,
             initialHighTransBw,
             this->lowTransBw,
             this->highTransBw);
        
        iLog("\t to ensure final overlapLength: ({}) <= frameLength: {}.",
             overlapLength,
             frameLength);
        
        iLog("\t (With original transition bandwidth parameters, overlapLength: {}).",
             initialOverlapLength);
    }
    //---------

    //Calculate array lengths
    convolutionResultLength = frameLength + filterLength - 1; //Even
    fftLength = nextPow2(convolutionResultLength);
    fftLengthDiv2 = fftLength / 2;

    /*
    iLog("FIRFilter::initialize() - Parameters:");
    iLog("\t sampleRate: {}", sampleRate);
    iLog("\t frameLength: {}", frameLength);
    iLog("\t lowCutoffFreq: {}", lowCutoffFreq);
    iLog("\t lowTransBw: {}", lowTransBw);
    iLog("\t highCutoffFreq: {}", highCutoffFreq);
    iLog("\t highTransBw: {}", highTransBw);
    iLog("\t overlapLength: {}", overlapLength);
    iLog("\t convolutionResultLength: {}", convolutionResultLength);
    iLog("\t fftLength: {}", fftLength);
    */
    
    //Allocate and initialize array memory to default (zero) values
    //+++++++++
    safeDeleteArray(filterPadded);
    filterPadded = new float[fftLength]();
    
    safeDeleteArray(signalPadded);
    signalPadded = new float[fftLength]();
    
    safeDeleteArray(signalPaddedXfilterPadded);
    signalPaddedXfilterPadded = new float[fftLength]();

    safeDeleteArray(overlapBuffer);
    overlapBuffer = new float[overlapLength]();
    
    
    safeDeleteArray(tempResult);
    tempResult = new float[fftLength]();
    
    safeDeleteArray(tempResultB);
    tempResultB = new float[fftLength]();

    
    safeDelete(fft);
    fft = new FFT(fftLength);
    //---------
    
    //Set up FFT parameters
    
    //Calculate FFT of (padded) filter kernel
    //+++++++++
    //Copy filter coefficients to array for transformation
    memcpy(filterPadded, //Destination
           filterCoefficients, //Source
           sizeof(float) * filterLength); //number of bytes
    
    //printArrayRow("filterPadded Before FFT", filterPadded, fftLength);

    //Calculate FFT of (padded) filter kernel, in place
    fft->forward(filterPadded);
    
    //printArray("filterPadded after FFT", filterPadded, 2, fftLength, true);
    
    //Get the nyquist real value
    //After a forward transform, the nyquist and 0 frequency bins are both
    //real-valued (no imaginary components). To save space (and keep to powers
    //of 2), the FFT format stores the real nyquist value in the (unused)
    //place of 0Hz frequency bin; (i.e. at index: 1).
    //To simplify future vector multiplication, we store the value and zero
    //the array component here:
    filterNyquistValue = filterPadded[1];
    filterPadded[1] = 0.0f; //For ease of subsequent complex multiplication
    //---------
    
    //Initialize prev/current initial transient sample counters.
    //Set prev to value != current, to avoid requiring an additional
    //conditional in applyFilter()
    initialTransientSampleCount = 0;
    prevInitialTransientSampleCount = 1;
    
    //Initialize insertion index for delay line filter
    delayLineInsertionIndex = 0;
    
    /*
    iLog("\n");
    iLog("Initial Filter:");
    
    if (lowCutoffFreq == VALUE_NOT_AVAILABLE) {
        iLog("\t Low Pass (Cutoff: {:4.2f}Hz, TransBW: {:4.2f}Hz)",
            highCutoffFreq,
            highTransBw);
    }
    else if (highCutoffFreq == VALUE_NOT_AVAILABLE) {
        iLog("\t High Pass (Cutoff: {:4.2f}Hz, TransBW: {:4.2f}Hz)",
            lowCutoffFreq,
            lowTransBw);
    }
    else {
        iLog("\t Band Pass: (CutoffLow:{:4.2f}Hz, TransBWLow:{:4.2f}Hz; CutoffHigh:{:4.2f}Hz, TransBWHigh:{:4.2f}Hz)",
            lowCutoffFreq,
            lowTransBw,
            highCutoffFreq,
            highTransBw);
    }
    */
    
    iLog("\t Filter Length: {} taps", filterLength);
    float filterDelay = float((float)filterLength / (float)sampleRate);
    iLog("\t Filter Delay: {0:.4f} sec", filterDelay);
    iLog("\n");
    
    if (overlapLength > frameLength) {
        iLog("**** WARNING: overlapLength ({}) > frameLength {}", overlapLength, frameLength);
        return false;
    }
    
    return true;
}



void FIRFilter::reset() {
    
    //Clear arrays that might need to be cleared
    std::fill_n(signalPadded, fftLength, 0.0f);
    std::fill_n(overlapBuffer, overlapLength, 0.0f);
    
    //Initialize prev/current initial transient sample counters.
    //Set prev to value != current, to avoid requiring an additional
    //conditional in applyFilter()
    initialTransientSampleCount = 0;
    prevInitialTransientSampleCount = 1;
    
    //Initialize insertion index for delay line filter
    delayLineInsertionIndex = 0;
}



//applyFilter()
//Can be in-place wrt input and output
//
void FIRFilter::applyFilter(const float inputFrameSamples[],
                            float outputFrameSamples[]) {
            
    applyFFTConvolutionFilter(inputFrameSamples, outputFrameSamples);
    
    //applyDelayLineFilter(inputFrameSamples, outputFrameSamples);
    //applyInputSideConvolutionFilter(inputFrameSamples, outputFrameSamples);
    //applyOutputSideConvolutionFilter(inputFrameSamples, outputFrameSamples);
    
    //printArrayRow("After applying filter", outputFrameSamples, frameLength);
    
    //Determine if output is past the initial transient
    if (prevInitialTransientSampleCount != initialTransientSampleCount) {
        
        prevInitialTransientSampleCount = initialTransientSampleCount;
        
        if (initialTransientSampleCount < frameLength) {
            initialTransientSampleCount += frameLength;
        }
    }
}



bool FIRFilter::pastInitialTransient() { //Called AFTER applyFilter call
    return prevInitialTransientSampleCount >= frameLength;
}



void FIRFilter::test(const float frequency,
                     const float duration) {
    
    const float sampleRate = 44100.0f;
    const uif32 numSamples = (uif32) ceil(sampleRate * duration);
    
    float signal[numSamples];
    float t[numSamples];
    
    iLog("Initial Filter Test ({0:4.5f}Hz Signal, {0:4.5f}sec long)",
         frequency,
         duration);
    
    //Calculate input waveform
    for (uif32 i = 0; i < numSamples; ++i) {
        t[i] = i * (duration / (float) numSamples);
        signal[i] = cos(2 * M_PI * frequency * t[i]);
    }
    
    //Print original waveform
    printv("Original Signal", signal, numSamples);
    
    //Conduct filtering, frame by frame
    uif32 filteredSampleCount = 0;
    for (uif32 i = 0; i + frameLength < numSamples; i += frameLength) {
        
        //Get pointer to start of frame
        float* signalFrame = signal + i;
        
        //Print before filter
        //printv("Before Filter", signalFrame, frameLength);
        
        //iLog("FRAME {} AFTER IFFT(FFT(s)):", i / frameLength);

        //Filter the frame
        applyFilter(signalFrame, signalFrame);
        
        //Print after filter
        //printv("After Filter", signalFrame, frameLength);
        
        filteredSampleCount += frameLength;
    }

    //Print final result
    printv("Final Result", signal, filteredSampleCount);
}



//applyDelayLineFilter() - O(N^2)
//
// Can be in-place wrt input and output
//
// O(N^2)
//
//Based on: http://www.eas.uccs.edu/~mwickert/ece5655/lecture_notes/ece5655_chap7.pdf
//p7-40
//
void FIRFilter::applyDelayLineFilter(const float inputFrameSamples[],
                                     float outputFrameSamples[]) {
    
    const uif32 filterLengthMinus1 = filterLength - 1;
    
    for (uif32 i = 0; i < frameLength; ++i) {
        
        //Update insertion index
        ++delayLineInsertionIndex;
        if (delayLineInsertionIndex >= filterLength) {
            delayLineInsertionIndex = 0;
        }
        
        //Insert newest entry
        tempResult[delayLineInsertionIndex] = inputFrameSamples[i];
        
        float s = 0.0f;
        sif32 tempI = delayLineInsertionIndex;
        for (uif32 j = 0; j < filterLength; ++j) {
            
            s += filterCoefficients[j] * tempResult[tempI];
            
            --tempI;
            if (tempI < 0) {
                tempI = filterLengthMinus1;
            }
        }
        
        outputFrameSamples[i] = s;
    }
}



// applyOutputSideConvolutionFilter()
//
// Can be in-place wrt input and output
//
// Makes use of kernel symetry for some performance optimization
//
// Filtered content has a more abrupt beginning, compared to
// input-side convolution, as described.
//
// (But how does performance compare?)
//
// O(N^2)
//
//https://christianfloisand.wordpress.com/tag/convolution/
//
void FIRFilter::applyOutputSideConvolutionFilter(const float inputFrameSamples[],
                                                 float outputFrameSamples[]) {
    
    uif32 i, j, k;
    
    const uif32 halfOrder = (filterLength-1) / 2;
    
    for (i = 0; i < frameLength; ++i) {
        
        tempResultB[0] = inputFrameSamples[i];
        tempResult[i] = 0.0f;
        
        for (j=0; j < halfOrder; ++j) {
            tempResult[i] += filterCoefficients[j] * (tempResultB[j] + tempResultB[filterLength-j-1]);
        }
        
        tempResult[i] += (filterCoefficients[j] * tempResultB[j]);
        
        for (k = filterLength - 1; k > 0; --k) {
            tempResultB[k] = tempResultB[k - 1];
        }
    }
    
    //Note comparison here between temp_result and convolution result from another
    //algorithm doesn't work; this is different...
    
    //Copy to output
    for (k = 0; k < frameLength; ++k) {
        outputFrameSamples[k] = tempResult[k];
    }
}



// applyFFTFilter()
//
// O(NlogN)
//
// Adapted from:
//   http://hamiltonkibbe.com/finite-impulse-response-filters-using-apples-accelerate-framework-part-iii/
//
void FIRFilter::applyFFTConvolutionFilter(const float inputFrameSamples[],
                                          float outputFrameSamples[]) {
    
    //Copy audio signal into beginning of a padded signal buffer
    memcpy(signalPadded, //Destination
           inputFrameSamples, //Source
           sizeof(float) * frameLength); //number of bytes
    
    //Fill the rest of the padded signal buffer with zeros
    vm_clrv(signalPadded + frameLength,
            fftLength - frameLength);
    
    
    //Calculate FFT of input signal
    //+++++++++
    fft->forward(signalPadded);
    //---------
    
    
    //Multiply signal and filter in the frequency domain
    //+++++++++
    //After a forward transform, the nyquist and 0 frequency bins are both
    //real-valued (no imaginary components). To save space (and keep to powers
    //of 2), the FFT stores the real nyquist value in the (unused) place of 0Hz
    //frequency bin; (i.e. at index: 1); which is a nuissance, because we now
    //need to move it in order facilitatate complex multiplication.
    //(The filter's nyquist component has already been moved and stored in
    //filterNyquistValue).
    float nyquistMultiplied = filterNyquistValue * signalPadded[1];
    
    //Clear nyquist component, to facilitate complex multiplication
    signalPadded[1] = 0.0f;
    
    //Complex multiplication of interleaved [R0,I0, R1,I1, ...] vectors
    vm_mulcv(signalPadded,
             filterPadded,
             signalPaddedXfilterPadded,
             fftLengthDiv2);

    //Restore nyquist component
    signalPaddedXfilterPadded[1] = nyquistMultiplied;
    //---------

    
    //Calculate inverse FFT of result
    //+++++++++
    fft->backward(signalPaddedXfilterPadded);
    //---------
    
    
    //Add in old overlap from previous function call
    vm_addv(signalPaddedXfilterPadded,
            overlapBuffer,
            signalPaddedXfilterPadded,
            overlapLength);
    
    
    //Copy new overlap into overlap buffer
    memcpy(overlapBuffer, //Destination
           signalPaddedXfilterPadded + frameLength, //Source
           sizeof(float) * overlapLength); //number of bytes
    
    
    //Copy final result to result audio frame
    memcpy(outputFrameSamples, //Destination
           signalPaddedXfilterPadded, //Source
           sizeof(float) * frameLength); //number of bytes
}



//setFilterCoefficientsAndLength
//
//For explanation and a simulation tool, see:
//https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
//https://tomroelandts.com/articles/how-to-create-simple-band-pass-and-band-reject-filters
//Additional good explanation here:
//http://www.analog.com/media/en/technical-documentation/dsp-book/dsp_book_Ch16.pdf
//http://www.physics.queensu.ca/~phys352/lect10.pdf
//
void FIRFilter::setFilterCoefficientsAndLength(const float sampleRate,
                                               const float lowCutoffFreq,
                                               const float lowTransBw,
                                               const float highCutoffFreq,
                                               const float highTransBw) {
    
    //Assert that the filter has at least one cutoff frequency - high or low
    assert((highCutoffFreq != VALUE_NOT_AVAILABLE) || (lowCutoffFreq != VALUE_NOT_AVAILABLE));
    
    //While the provided cutoff frequencies are where the signal *starts* to dip,
    //FIR filter design method calls for high/low frequencies at the half-attenuation point.
    //Calculate the half-attenuation points.
    float fhaHigh = highCutoffFreq + highTransBw / 2.0f;
    float fhaLow = lowCutoffFreq - lowTransBw / 2.0f;
    
    //High frequency for specifying low-pass filter
    float fhaHighDivFs = fhaHigh / sampleRate;  //Fraction of the sampling rate (in (0, 0.5))
    if (fhaHighDivFs > 0.5f) {
        eLog("calcFilterCoefficients - ERROR: High cutoff frequency is too high.");
    }
    
    //Low frequency for specifying high-pass filter
    float fhaLowDivFs = fhaLow / sampleRate;  //Fraction of the sampling rate (in (0, 0.5))
    if (fhaLowDivFs > 0.5f) {
        eLog("calcFilterCoefficients - ERROR: Low cutoff frequency is too high.");
    }
    
    //Low cutoff transition band-width, for specifying high-pass filter
    float ltbwDivFs = lowTransBw / sampleRate;  //Fraction of the sampling rate (in (0, 0.5)).
    if (ltbwDivFs > 0.5f) {
        eLog("calcFilterCoefficients - ERROR: Low transition bandwidth is too high.");
    }
    
    //High cutoff transition band-width, for specifying low-pass filter
    float htbwDivFs = highTransBw / sampleRate;  //Fraction of the sampling rate (in (0, 0.5)).
    if (htbwDivFs > 0.5f) {
        eLog("calcFilterCoefficients - ERROR: High transition bandwidth is too high.");
    }
    
    //Determine low-pass filter length
    uif32 lpLen = (uif32) ceil(4.0f / htbwDivFs);
    if ((lpLen & 1) == 0) { //Length must be odd
        lpLen++;
    }
    
    //Determine high-pass filter length
    uif32 hpLen = (uif32) ceil(4.0f / ltbwDivFs);
    if ((hpLen & 1) == 0) { //Length must be odd
        hpLen++;
    }
    
    //Construct the low-pass filter
    float lp[lpLen];
    if (highCutoffFreq != VALUE_NOT_AVAILABLE) {
        
        calcLowPassFilterCoefficients(fhaHighDivFs, lp, lpLen);
        
        //If filter is low-pass ONLY, finish up and return
        if (lowCutoffFreq == VALUE_NOT_AVAILABLE) {
            
            //Finalize filter length
            filterLength = lpLen;
            
            //Finalize filter coefficients
            safeDeleteArray(filterCoefficients);
            filterCoefficients = new float[filterLength];
            for (uif32 i = 0; i < filterLength; ++i) {
                filterCoefficients[i] = lp[i];
            }
            
            //printFilterCoefficients();
            
            return;
        }
    }
    
    //Construct high-pass filter
    float hp[hpLen];
    if (lowCutoffFreq != VALUE_NOT_AVAILABLE) {
        
        calcHighPassFilterCoefficients(fhaLowDivFs, hp, hpLen);
        
        //If filter is high-pass ONLY, finish up and return
        if (highCutoffFreq == VALUE_NOT_AVAILABLE) {
            
            //Finalize filter length
            filterLength = hpLen;
            
            //Finalize filter coefficients
            safeDeleteArray(filterCoefficients);
            filterCoefficients = new float[filterLength];
            for (uif32 i = 0; i < filterLength; ++i) {
                filterCoefficients[i] = hp[i];
            }
            
            //printFilterCoefficients();
            
            return;
        }
    }
    
    //Convolve high and low-pass filters to get band-pass filter
    //+++++++++++
    //Finalize filter length
    filterLength = lpLen + hpLen - 1;
    
    //Finalize filter coefficients
    safeDeleteArray(filterCoefficients);
    filterCoefficients = new float[filterLength];
    convolve(hp, hpLen, lp, lpLen, filterCoefficients);
    //-----------
    
    //Print coefficients
    //printFilterCoefficients();
}



void FIRFilter::calcLowPassFilterCoefficients(const float fhaHighDivFs, //Half-attenuation frequency / sampling frequency
                                              float coefficients[],
                                              const uif32 coefficientsLength) {
    
    const uif32 length = coefficientsLength;
    uif32 i;
    
    //Compute sinc
    float sl[length];
    for (i = 0; i < length; ++i) {
        sl[i] = sinc(2 * fhaHighDivFs * (i - (length - 1) / 2.0f));
    }
    
    //Compute blackman window
    float bl[coefficientsLength];
    for (i = 0; i < length; ++i) {
        bl[i] = 0.42 - 0.5 * cos(2.0f * M_PI * i / (length - 1)) + 0.08f * cos(4.0f * M_PI * i / (length - 1));
    }
    
    //Multiply sinc by blackman window
    for (i = 0; i < length; ++i) {
        coefficients[i] = sl[i] * bl[i];
    }
    
    //Get sum of filter coefficients
    float sum = 0.0;
    for (i = 0; i < length; ++i) {
        sum += coefficients[i];
    }
    
    //Normalize for unity gain
    for (i = 0; i < length; ++i) {
        coefficients[i] = coefficients[i] / sum;
    }
}



void FIRFilter::calcHighPassFilterCoefficients(const float fhaLowDivFs, //Half-attenuation frequency / sampling frequency
                                               float coefficients[],
                                               const uif32 coefficientsLength) {
    
    const uif32 length = coefficientsLength;
    uif32 i;
    
    //Compute sinc
    float s[length];
    for (i = 0; i < length; ++i) {
        s[i] = sinc(2 * fhaLowDivFs * (i - (length - 1) / 2.0));
    }
    
    //Compute blackman window
    float b[length];
    for (i = 0; i < length; ++i) {
        b[i] = 0.42f - 0.5f * cos(2.0f * M_PI * i / (length - 1)) + 0.08f * cos(4.0f * M_PI * i / (length - 1));
    }
    
    //Multiply sinc by blackman window
    for (i = 0; i < length; ++i) {
        coefficients[i] = s[i] * b[i];
    }
    
    //Get sum of filter coefficients
    float sum = 0.0f;
    for (i = 0; i < length; ++i) {
        sum += coefficients[i];
    }
    
    //Normalize for unity gain
    for (i = 0; i < length; ++i) {
        coefficients[i] = coefficients[i] / sum;
    }
    
    //Spectrally invert (negate samples, add 1 to center sample)
    for (i = 0; i < length; ++i) {
        coefficients[i] = -coefficients[i];
    }
    coefficients[(length - 1) / 2] += 1;
}



void FIRFilter::printFilterCoefficients() {
    
    printf("Initial Filter Kernel:\n");
    
    printf("\tLength:%u Fcl:%4.2f Tbwl:%4.2f Fch:%4.2f Tbwh:%4.2f\n",
           filterLength,
           lowCutoffFreq,
           lowTransBw,
           highCutoffFreq,
           highTransBw);
    printf("--------------\n");
    
    printv("Coefficients", filterCoefficients, filterLength);
    
    printf("\n--------------\n");
}
