//Utility.cpp

#include "Utility.hpp"

#include "Logger.hpp"
#include "TypeAbbreviations.hpp"

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>



//Timing-related
//++++++++++++++++++++
void waitMs(uif32 durationMs) {
    ChronoSteadyClockTimePoint start = ChronoSteadyClock::now();
    while (true) {
        ChronoDurationFloatMs duration( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
        if (duration.count() > durationMs) {
            break;
        }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(durationMs));
}

void sleepMs(uif32 durationMs) {
    std::this_thread::sleep_for(std::chrono::milliseconds(durationMs));
}

void operationTimingTest() {
    
    printf("Operation Timing Test:\n\n");
    
    uif32 reps = 1000;
    
    printf("\t Iterations for each operation: %u\n", reps);
    
    ChronoSteadyClockTimePoint start = ChronoSteadyClock::now();
    for (float i = 0; i < reps; ++i) {
        float x = i + (i-1);
        x = 0.0f;
    }
    ChronoDurationFloatMs addsDurationMs( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
    float addsDurationCountMs = addsDurationMs.count();
    printf("\t %f ms / addition\n", addsDurationCountMs/reps);
    
    start = ChronoSteadyClock::now();
    for (float i = 0; i < reps; ++i) {
        float x = i * (i-1);
        x = 0.0f;
    }
    ChronoDurationFloatMs multsDurationMs( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
    float multsDurationCountMs = multsDurationMs.count();
    printf("\t %f ms / multiplication\n", multsDurationCountMs/reps);
    
    start = ChronoSteadyClock::now();
    for (float i = 0; i < 1000; ++i) {
        float x = expf(i);
        x = 0.0f;
    }
    ChronoDurationFloatMs expsDurationMs( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
    float expsDurationCountMs = expsDurationMs.count();
    printf("\t %f ms / expf()\n", expsDurationCountMs/reps);
    
    start = ChronoSteadyClock::now();
    for (float i = 0; i < 1000; ++i) {
        float x = powf(i, i);
        x = 0.0f;
    }
    ChronoDurationFloatMs powsDurationMs( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
    float powsDurationCountMs = powsDurationMs.count();
    printf("\t %f ms / powf()\n", powsDurationCountMs/reps);
    
    start = ChronoSteadyClock::now();
    for (float i = 0; i < 1000; ++i) {
        float x = logf(i);
        x = 0.0f;
    }
    ChronoDurationFloatMs logsDurationMs( ChronoDurationCast<ChronoDurationFloatMs>(ChronoSteadyClock::now() - start) );
    float logsDurationCountMs = logsDurationMs.count();
    printf("\t %f ms / logf()\n", logsDurationCountMs/reps);
    
    printf("\n");
    
    printf("\t => %f multiplications()s / addition\n", multsDurationCountMs / addsDurationCountMs);
    
    printf("\t => %f expf()s / multiplication\n", expsDurationCountMs / multsDurationCountMs);
    
    printf("\t => %f powf()s / multiplication\n", powsDurationCountMs / multsDurationCountMs);
    
    printf("\t => %f logf()s / multiplication\n", logsDurationCountMs / multsDurationCountMs);
}
//--------------------



//Waveform-related
//++++++++++++++++++++
void createCosineWave(const float frequency,
                      const float sampleRate,
                      float result[],
                      const uif32 resultLength) {
    
    for (uif32 i = 0; i < resultLength; ++i) {
        result[i] = cos(2 * M_PI * frequency * i / sampleRate);
    }
}
//--------------------

    
    
//Power/Log
//++++++++++++++++++++
float sqr(const float x) {
	return x * x;
}



bool isPowerOf2(const ui32 x) {
    return x > 0 && !(x & (x - 1));
}



uif32 nextPow2(uif32 x) {
    
    --x;
    
    for (uif32 i = 1; i < sizeof(x) * 8; i *= 2) {
        x |= x >> i;
    }
    
    ++x;
    
    return x;
}



uif32 logBase2(const uif32 x) {
    
    switch (x) {
            
        case 32:
            return 5;
            break;
            
        case 64:
            return 6;
            break;
            
        case 128:
            return 7;
            break;
            
        case 256:
            return 8;
            break;
            
        case 512:
            return 9;
            break;
            
        case 1024:
            return 10;
            break;
            
        case 2048:
            return 11;
            break;
            
        case 4096:
            return 12;
            break;
            
        case 8192:
            return 13;
            break;
            
        case 16384:
            return 14;
            break;
    }
    
    eLog("Error in logBase2() - Invalid argument given: {}", x);
    
    //Invalid argument; return 0
    return 0;
}
//--------------------




//Interpolation
//++++++++++++++++++++

//parabolicLogMaxInterp()
//Interpolates the maximum and fractional index of the maximum of a parabolic fit through three points
//NOTE1: Interpolation takes place on a LOG scale, for a more exact parabolic fit.
//NOTE2: Window audio data using a GAUSSIAN function with steep fall-off for greatest accuracy.
//Math based on:
// https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_qint_m.html
// https://gist.github.com/endolith/255291
// Google "Improving FFT Frequency Measurement Resolution by Parabolic and Gaussian Spectrum Interpolation" paper
//Structure based on:
// "fvec_quadratic_peak_pos()" function in "aubio" code library (mathutils.c)
// Available here http://dev.aubio.org
//Also see:
// "parabolicInterpolation()" function in "TarsoDSP" code library (Yin.java)
//Available here: https://github.com/JorenSix/TarsosDSP
//
void parabolicLogMaxInterp(const float data[],
                           const uif32 closestWholeIndex, //Closest whole index
                           const uif32 length,
                           float* interpolatedMaximum,
                           float* interpolatedMaximumFractionalIndex) {
    
    //If the selected index is an end-point...
    //(Not making extrapalatory assumptions)
    if ( closestWholeIndex == 0 || closestWholeIndex == (length - 1) ) {
        *interpolatedMaximum = data[closestWholeIndex];
        *interpolatedMaximumFractionalIndex = closestWholeIndex;
        return;
    }
    
    //Determine previous and next whole indices, being safe about array bounds
    uif32 prevWholeIndex = (closestWholeIndex < 1) ? closestWholeIndex : closestWholeIndex - 1;
    uif32 nextWholeIndex = (closestWholeIndex + 1 < length) ? closestWholeIndex + 1 : closestWholeIndex;
    
    if (prevWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[nextWholeIndex]) {
            *interpolatedMaximum = data[nextWholeIndex];
            *interpolatedMaximumFractionalIndex = nextWholeIndex;
            return;
        }
        else {
            *interpolatedMaximum = data[closestWholeIndex];
            *interpolatedMaximumFractionalIndex = closestWholeIndex;
            return;
        }
    }
    
    if (nextWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[prevWholeIndex]) {
            *interpolatedMaximum = data[prevWholeIndex];
            *interpolatedMaximumFractionalIndex = prevWholeIndex;
            return;
        }
        else {
            *interpolatedMaximum = data[closestWholeIndex];
            *interpolatedMaximumFractionalIndex = closestWholeIndex;
            return;
        }
    }
    
    //Get samples on log scale
    float s0 = logf(data[prevWholeIndex]);
    float s1 = logf(data[closestWholeIndex]);
    float s2 = logf(data[nextWholeIndex]);
    
    //Equations from: https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_qint_m.html
    float deltaFromClosestWholeIndex = (s2 - s0) / (2 * (2 * s1 - s2 - s0));
    *interpolatedMaximum = s1 - (0.25 * (s0 - s2) * deltaFromClosestWholeIndex);
    //Convert back to linear scale
    *interpolatedMaximum = expf(*interpolatedMaximum);
    *interpolatedMaximumFractionalIndex = (float)closestWholeIndex + deltaFromClosestWholeIndex;
}


//parabolicMaxInterp()
//Interpolates the maximum and fractional index of the maximum of a parabolic fit through three points
//Math based on:
// https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_qint_m.html
//Structure based on:
// "fvec_quadratic_peak_pos()" function in "aubio" code library (mathutils.c)
// Available here http://dev.aubio.org
//Also see:
// "parabolicInterpolation()" function in "TarsoDSP" code library (Yin.java)
//Available here: https://github.com/JorenSix/TarsosDSP
//
void parabolicMaxInterp(const float data[],
                        const uif32 closestWholeIndex,
                        const uif32 length,
                        float* interpolatedMaximum,
                        float* interpolatedMaximumFractionalIndex) {
    
    //If the selected index is an end-point...
    if ( closestWholeIndex == 0 || closestWholeIndex == length - 1 ) {
        *interpolatedMaximum = data[closestWholeIndex];
        *interpolatedMaximumFractionalIndex = closestWholeIndex;
        return;
    }
    
    //Determine prev and next whole indices, being safe about array bounds
    uif32 prevWholeIndex = (closestWholeIndex < 1) ? closestWholeIndex : closestWholeIndex - 1;
    uif32 nextWholeIndex = (closestWholeIndex + 1 < length) ? closestWholeIndex + 1 : closestWholeIndex;
    
    if (prevWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[nextWholeIndex]) {
            *interpolatedMaximum = data[nextWholeIndex];
            *interpolatedMaximumFractionalIndex = nextWholeIndex;
            return;
        }
        else {
            *interpolatedMaximum = data[closestWholeIndex];
            *interpolatedMaximumFractionalIndex = closestWholeIndex;
            return;
        }
    }
    
    if (nextWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[prevWholeIndex]) {
            *interpolatedMaximum = data[prevWholeIndex];
            *interpolatedMaximumFractionalIndex = prevWholeIndex;
            return;
        }
        else {
            *interpolatedMaximum = data[closestWholeIndex];
            *interpolatedMaximumFractionalIndex = closestWholeIndex;
            return;
        }
    }
    
    float s0 = data[prevWholeIndex];
    float s1 = data[closestWholeIndex];
    float s2 = data[nextWholeIndex];
    
    //Equations from: https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_qint_m.html
    float deltaFromClosestWholeIndex = (s2 - s0) / (2 * (2 * s1 - s2 - s0));
    *interpolatedMaximum = s1 - 0.25 * (s0 - s2) * deltaFromClosestWholeIndex;
    *interpolatedMaximumFractionalIndex = (float)closestWholeIndex + deltaFromClosestWholeIndex;
}


//parabolicMinInterp()
//Interpolates the minimum and fractional index of the maximum of a parabolic fit through three points
//Math based on:
// https://ccrma.stanford.edu/~jos/sasp/Matlab_listing_qint_m.html
//Structure based on:
// "fvec_quadratic_peak_pos()" function in "aubio" code library (mathutils.c)
// Available here http://dev.aubio.org
//Also see:
// "parabolicInterpolation()" function in "TarsoDSP" code library (Yin.java)
//Available here: https://github.com/JorenSix/TarsosDSP
//
void parabolicMinInterp(const float data[],
                        const uif32 closestWholeIndex,
                        const uif32 vlen,
                        float* interpolatedMinimum,
                        float* interpolatedMinimumFractionalIndex) {
    
    
    //If the selected index is an end-point...
    if ( closestWholeIndex == 0 || closestWholeIndex == vlen - 1 ) {
        *interpolatedMinimum = data[closestWholeIndex];
        *interpolatedMinimumFractionalIndex = closestWholeIndex;
        return;
    }
    
    //Determine prev and next whole indices, being safe about array bounds
    uif32 prevWholeIndex = (closestWholeIndex < 1) ? closestWholeIndex : closestWholeIndex - 1;
    uif32 nextWholeIndex = (closestWholeIndex + 1 < vlen) ? closestWholeIndex + 1 : closestWholeIndex;
    

    if (prevWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[nextWholeIndex]) {
            *interpolatedMinimum = data[closestWholeIndex];
            *interpolatedMinimumFractionalIndex = closestWholeIndex;
            return;
        }
        else {
            *interpolatedMinimum = data[nextWholeIndex];
            *interpolatedMinimumFractionalIndex = nextWholeIndex;
            return;
        }
    }
    
    if (nextWholeIndex == closestWholeIndex) {
        if (data[closestWholeIndex] <= data[prevWholeIndex]) {
            *interpolatedMinimum = data[closestWholeIndex];
            *interpolatedMinimumFractionalIndex = closestWholeIndex;
            return;
        }
        else {
            *interpolatedMinimum = data[prevWholeIndex];
            *interpolatedMinimumFractionalIndex = prevWholeIndex;
            return;
        }
    }
    
    float s0 = data[prevWholeIndex];
    float s1 = data[closestWholeIndex];
    float s2 = data[nextWholeIndex];
    
    float deltaFromClosestWholeIndex = (s2 - s0) / (2 * (2 * s1 - s2 - s0) );
    
    *interpolatedMinimum = s1 - 0.25f * (s0-s2) * deltaFromClosestWholeIndex;
    
    *interpolatedMinimumFractionalIndex = (float)closestWholeIndex + deltaFromClosestWholeIndex;
}



float linearInterp(const float data[],
                   const uif32 length,
                   const float fractionalIndex) {
    
    long indexBelow = floor(fractionalIndex);
    long indexAbove = ceil(fractionalIndex);

    //If interpolation would lead out of bounds...
    if (indexBelow < 0) {
        return data[0];
    }
    else if (indexAbove > length - 1) {
        return data[length - 1];
    }
    
    //If interpolation is unnecessary...
    if (fractionalIndex == indexBelow) {
        return data[indexBelow];
    }
    else if (fractionalIndex == indexAbove) {
        return data[indexAbove];
    }
    
    //Interpolate
    return data[indexBelow] + (data[indexAbove] - data[indexBelow]) * ((fractionalIndex - indexBelow) / (indexAbove - indexBelow));
}



//sincInterp()
//Interpolates the value of the function given in "input",
//at the index "fractionalIndex". See resampleSinc() below.
float sinc(float x) {
    x *= M_PI;
    return (x == 0.0f) ? 1.0f : (sin(x) / x);
}
float sincInterp(const float data[],
                 uif32 dataLength,
                 float fractionalIndex) {
    float result = 0.0f;
    for (uif32 i = 0; i < dataLength; ++i) {
        result += data[i] * sinc(fractionalIndex - i);
    }
    return result;
}
//--------------------



//Resampling / Remapping
//++++++++++++++++++++
//resampleSinc()
//Resample (upsample or downsample) using sinc interpolation.
//Sinc interpolation takes more time than linear interpolation, but:
// * Is smooth (good for audio)
// * Ensures the resampled function passes through the original points
//See:
//  http://stackoverflow.com/a/31878492
//Other resources:
//  https://ccrma.stanford.edu/~jos/resample/Free_Resampling_Software.html
//  http://avisynth.nl/index.php/Resampling
//  http://stackoverflow.com/questions/31836598/subsampling-an-array-of-numbers
//  https://ccrma.stanford.edu/~jos/resample/What_Bandlimited_Interpolation.html
//
void resampleSinc(const float input[],
                  uif32 inputLength,
                  float output[],
                  uif32 outputLength) {
    
    float dx = (float) (inputLength - 1) / (outputLength - 1);
    
    for (uif32 i = 0; i < outputLength; ++i) {
        output[i] = sincInterp(input, inputLength, i * dx);
    }
}

void resampleSincTest() {
    
    //Input
    float input[] = {
        0.0,
        1.0,
        0.5,
        0.2,
        0.1,
        0.0,
    };
    const uif32 inputLength = sizeof(input) / sizeof(input[0]);
    
    //Output
    const uif32 outputLength = 20;
    float output[outputLength];
    
    //Print input
    printf("inputIndices = [");
    for (uif32 i = 0; i < inputLength; ++i) {
        printf("%u ", i);
    }
    printf("];\n\n");
    printf("inputValues = [");
    for (uif32 i = 0; i < inputLength; ++i) {
        printf("%.6f ", input[i]);
    }
    printf("];\n\n");
    
    //Resample
    resampleSinc(input,
                 inputLength,
                 output,
                 outputLength);
    
    //Print output
    printf("outputIndices = [");
    for (uif32 i = 0; i < outputLength; ++i) {
        float dx = (float) (inputLength - 1) / (outputLength - 1);
        printf("%.6f ", (float) i * dx);
    }
    printf("];\n\n");
    printf("outputValues = [");
    for (uif32 i = 0; i < outputLength; ++i) {
        printf("%.6f ", output[i]);
    }
    printf("];\n\n");
    
    printf("plot(inputIndices, inputValues, \'r:\', "
           "outputIndices, outputValues, \'b:\');\n\n");
    
}



//See https://www.arduino.cc/en/Reference/Map
float map(const float inValue,
          const float inRangeMinimum,
          const float inRangeMaximum,
          const float outRangeMinimum,
          const float outRangeMaximum) {
    
    return (inValue - inRangeMinimum) * (outRangeMaximum - outRangeMinimum) / (inRangeMaximum - inRangeMinimum) + outRangeMinimum;
}
//--------------------



//Sorting
//++++++++++++++++++++

//QuickSort for float array
//Adapted from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
void quickSortFloatArray(float data[],
                         const uif64 length) {
    
	if (length < 2) {
		return;
	}
    
	float pivot = data[length / 2];
	float *left = data;
	float *right = data + length - 1;
    
	while (left <= right) {
		if (*left < pivot) {
			++left;
		}
		else if (*right > pivot) {
			--right;
		}
		else {
			float temp = *left;
			*left = *right;
			*right = temp;
			++left;
			--right;
		}
	}
    
	quickSortFloatArray(data, right - data + 1);
	quickSortFloatArray(left, data + length - left);
}



//QuickSort for uif32 array
//Adapted from http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#C
void quickSortUIF32Array(uif32 data[],
                         const uif64 length) {
    
	if (length < 2) {
		return;
	}
    
	uif32 p = data[length / 2];
	uif32 *left = data;
	uif32 *right = data + length - 1;
	while (left <= right) {
		if (*left < p) {
			++left;
		}
		else if (*right > p) {
			--right;
		}
		else {
			uif32 temp = *left;
			*left = *right;
			*right = temp;
			++left;
			--right;
		}
	}
    
	quickSortUIF32Array(data, right - data + 1);
	quickSortUIF32Array(left, data + length - left);
}


//InsertionSort for float array
//Sorts in place
//Adapted from https://stat.ethz.ch/pipermail/r-devel/attachments/20100223/077b46e0/attachment.c
//Coded by David Simcha
uif32 insertionSortFloatArray(float data[],
                              const uif32 length) {
	
	uif32 i, jMax;
	uif32 swapCount = 0;
	
	if (length < 2) {
		return 0;
	}
	
	jMax = length - 1;
	for (i = length - 2; i < length; --i) {
        
		uif32 j = i;
		float val = data[i];
		
		for (; j < jMax && data[j + 1] < val; ++j) {
			data[j] = data[j + 1];
		}
		
		data[j] = val;
		swapCount += (j - i);
	}
	
	return swapCount;
}


//Sorts indices of a float array, such that indices reference ascending
//values of the float array. Length of indices can be <= length of data
//array.
void insertionSortIndicesOfFloatArray(const float data[],
                                      uif32 indices[],
                                      const uif32 indicesLength) {
    uif32 i, jMax;
    
    if (indicesLength < 2) {
        return;
    }
    
    jMax = indicesLength - 1;
    for (i = indicesLength - 2; i < indicesLength; --i) {

        uif32 j = i;

        float val = data[indices[i]];
        float valIndex = indices[i];
        
        for (; j < jMax && data[indices[j + 1]] < val; ++j) {
            indices[j] = indices[j + 1];
        }
        
        indices[j] = valIndex;
    }
}



//Merge sort, for float array
//Sorts in place
//Adapted from https://stat.ethz.ch/pipermail/r-devel/attachments/20100223/077b46e0/attachment.c
//Coded by David Simcha
//+++++++++++++++++
static unsigned long mergeFloatArrays(float* from,
                                      float* to,
                                      const uif32 middle,
                                      const uif32 length) {
	
	uif32 bufferIndex, leftLength, rightLength;
	uif32 swaps;
	float* left;
	float* right;
	
	bufferIndex = 0;
	swaps = 0;
	
	left = from;
	right = from + middle;
	rightLength = length - middle;
	leftLength = middle;
	
	while (leftLength && rightLength) {
		
		if (right[0] < left[0]) {
			to[bufferIndex] = right[0];
			swaps += leftLength;
			--rightLength;
			++right;
		}
		else {
			to[bufferIndex] = left[0];
			--leftLength;
			++left;
		}
		++bufferIndex;
	}
	
	if (leftLength) {
		memcpy(to + bufferIndex, left, leftLength * sizeof(float));
	}
	else if (rightLength) {
		memcpy(to + bufferIndex, right, rightLength * sizeof(float));
	}
	
	return swaps;
}



uif32 mergeSortFloatArray(float data[],
                          float buffer[],
                          const uif32 length) {
	
	if (length < 10) {
		return insertionSortFloatArray(data, length);
	}
    
	uif32 lengthDiv2 = length / 2;
    
	uif32 swaps = mergeSortFloatArray(data,
                                      buffer,
                                      lengthDiv2);
    
	swaps += mergeSortFloatArray(data + lengthDiv2,
                                 buffer + lengthDiv2,
                                 length - lengthDiv2);
    
	swaps += mergeFloatArrays(data,
                              buffer,
                              lengthDiv2,
                              length);
	
	memcpy(data, buffer, length * sizeof(float));
	
	return swaps;
}
//--------------------



//Filtering
//++++++++++++++++++++



//Median Filter
//See: http://www.librow.com/articles/article-1

void medianFilterInner(const float extendedSignal[],
                       float result[],
                       const uif32 extendedSignalLength,
                       const uif32 windowLength) {
    
    const uif32 windowLengthDiv2 = windowLength / 2;
    
    //Move window through all elements of the signal
    for (uif32 i = windowLengthDiv2;
            i < extendedSignalLength - windowLengthDiv2; ++i) {
        
        //Pick up window elements
        float window[windowLength];
        for (uif32 j = 0; j < windowLength; ++j) {
            window[j] = extendedSignal[i - windowLengthDiv2 + j];
        }
        
        //Find and store median result
        result[i - windowLengthDiv2] = fastMedian(window, windowLength);
    }
}



void medianFilter(const float inputData[],
                  float outputData[],
                  const uif32 dataLength,
                  const uif32 windowLength) {
    
    //Check arguments
    if (!inputData || dataLength < 1) {
        return;
    }
    
    //Special case: windowLength <= 1
    if (windowLength <= 1) {
        memcpy(outputData,
               inputData,
               sizeof(float) * dataLength);
        return;
    }
    
    //Special case: n == 1
    if (dataLength == 1) {
        if (outputData) {
            outputData[0] = inputData[0];
        }
        return;
    }
    
    //Allocate "extended" input: input, with windowLength/2 on both ends
    const uif32 windowLengthDiv2 = windowLength / 2;
    const uif32 extendedInputDataLength = dataLength + (2 * windowLengthDiv2);
    float* extendedInputData = new float[extendedInputDataLength];
    
    //Copy signal to middle of extended signal
    memcpy(extendedInputData + windowLengthDiv2,
           inputData,
           dataLength * sizeof(float));
    
    //Copy first and last windowLengthDiv2 samples of signal, reflected,
    //to start and end of signalExtended
    for (uif32 i = 0; i < windowLengthDiv2; ++i) {
        //Start
        uif32 j = (windowLengthDiv2 - 1 - i) % dataLength;
        extendedInputData[i] = inputData[j];
        //End
        uif32 k = (dataLength - 1 - i) % dataLength;
        extendedInputData[dataLength + windowLengthDiv2 + i] = inputData[k];
    }
    
    //Call median filter implementation
    medianFilterInner(extendedInputData,
                      outputData,
                      extendedInputDataLength,
                      windowLength);
    
    //Free memory
    delete[] extendedInputData;
}
//--------------------



//Comparison
//++++++++++++++++++++
//NOT FOOL-PROOF; see: http://stackoverflow.com/a/15012792
bool approxEqualAbsTolerance(const float x,
                             const float y) {
    
    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon();
}



bool approxEqualRelTolerance(const float x,
                             const float y) {
    
    float maxXY = std::max( std::fabs(x) , std::fabs(y) );
    
    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon()*maxXY ;
}



bool approxEqualWithinDiffFraction(const float x,
                                   const float y,
                                   const float diffFrac) {
    
    return diffFraction(x, y) <= diffFrac;
}
//--------------------



//Stats
//++++++++++++++++++++

//Fast median routine known as "QuickSelect", used median filtering.
//See: http://ndevilla.free.fr/median/median.pdf
//
//* Does not preserve input array
//* If input array is even, result is the lower of the two middle values
//
#define SWAP(a, b) {float t=(a); (a) = (b); (b) = t;}
float fastMedian(float data[],
                 const uif32 dataLength) {
    
    uif32 low, high, median; //indices
    uif32 middle, ll, hh; //indices
    
    low = 0;
    high = dataLength - 1;
    median = (low + high) / 2;
    
    for (;;) {
        
        // One element only
        if (high <= low) {
            return data[median];
        }
        
        //Two elements only
        if (high == low + 1) {
            if (data[low] > data[high]) {
                SWAP(data[low], data[high]);
            }
            return data[median];
        }
        
        //Find median of low, middle and high items; swap into position low
        middle = (low + high) / 2;
        if (data[middle] > data[high]) {
            SWAP(data[middle], data[high]);
        }
        if (data[low] > data[high]) {
            SWAP(data[low], data[high]);
        }
        if (data[middle] > data[low]) {
            SWAP(data[middle], data[low]);
        }
        
        //Swap low item (now in position middle) into position (low + 1)
        SWAP(data[middle], data[low+1]) ;
        
        ///Nibble from each end towards middle, swapping items when stuck
        ll = low + 1;
        hh = high;
        for (;;) {
            
            do {
                ++ll;
            } while (data[low] > data[ll]);
            
            do {
                --hh;
            } while (data[hh] > data[low]);
            
            if (hh < ll) {
                break;
            }
            
            SWAP(data[ll], data[hh]);
        }
        
        //Swap middle item (in position low) back into correct position
        SWAP(data[low], data[hh]);
        
        //Re-set active partition
        if (hh <= median) {
            low = ll;
        }
        
        if (hh >= median) {
            high = hh - 1;
        }
    }
}
#undef SWAP



//Median Routines from http://rosettacode.org/wiki/Averages/Median
//++++++++
float getMedianWithoutPresevingInputArray(float data[],
                                          const uif32 length) {
	
	quickSortFloatArray(data, length);
	
    if (length % 2 == 0) {
        return (data[length/2] + data[length/2 - 1]) / 2.0f;
    }
    
	return data[length/2];
}



float getMedianWhilePresevingInputArray(const float data[],
                                        const uif32 length) {
	
	float a[length];
	for (uif32 i = 0; i < length; ++i) {
		a[i] = data[i];
	}
	
	quickSortFloatArray(a, length);
	
    if (length % 2 == 0) {
        return (a[length/2] + a[length/2 - 1]) / 2.0f;
    }
    
    return a[length/2];
}
//--------



//linearRegression()
//
//See: http://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c
//
bool linearRegression(const float dataX[],
                      const float dataY[],
                      const uif32 length,
                      float* slope,
                      float* yIntercept,
                      float* correlationCoefficient) { //Optional
	
	float sumX = 0.0f; //Sum of dataX
	float sumX2 = 0.0f; //Sum of X^2
	float sumXY = 0.0f; //Sum of X * Y
	float sumY = 0.0f; //Sum of y
	float sumY2 = 0.0f; //Sum of Y^2
	
	for (uif32 i = 0; i < length; ++i) {
		sumX  += dataX[i];
		sumX2 += sqr(dataX[i]);
		sumXY += dataX[i] * dataY[i];
		sumY  += dataY[i];
		sumY2 += sqr(dataY[i]);
	}
	
	float denominator = (length * sumX2 - sqr(sumX));
	if (denominator == 0) {
		//Singular matrix; can't solve the problem...
		*slope = 0.0f;
		*yIntercept = 0.0f;
		*correlationCoefficient = 0.0f;
		return false;
	}
	
	*slope = (length * sumXY  -  sumX * sumY) / denominator;
	*yIntercept = (sumY * sumX2  -  sumX * sumXY) / denominator;
	if (correlationCoefficient != nullptr) {
		*correlationCoefficient = (sumXY - sumX * sumY / length) /
            sqrt( (sumX2 - sqr(sumX) / length) * (sumY2 - sqr(sumY) / length) );
	}
	
	return true;
}



//erfcc()
//From "Numerical Recipies, The Art of Scientific Computing, Second Edition", p221
//http://www2.units.it/ipl/students_area/imm2/files/
//Returns the complementary error function erfc(x) with fractional error everywhere less than 1.2 × 10−7.
//This function is used in kendl1()
float erfcc(const float x) {
	float t, z, ans;
	z = fabs(x);
	t = 1.0f / (1.0f + (0.5f * z));
    ans = t * exp(-z * z - 1.26551223f +
                  t * (1.00002368f +
                       t * (0.37409196f +
                            t * (0.09678418f +
                                 t * (-0.18628806f +
                                      t * (0.27886807f +
                                           t * (-1.13520398f +
                                                t * (1.48851587f +
                                                     t * (-0.82215223f +
                                                          t * 0.17087277f)))))))));
    return x >= 0.0f ? ans : 2.0f - ans;
}



//kendallN2()
//Calculate Kendall's Tau - O(n^2)
//From "Numerical Recipies, The Art of Scientific Computing, Second Edition", p667
//http://www2.units.it/ipl/students_area/imm2/files/
//
//No sorting of input data required before-hand
//
//Given data arrays data1[1..n] and data2[1..n], returns: Kendall’s tau, its number of standard deviations from zero
//as z, and its two-sided significance level as prob.
//
//Small values of prob indicate a significant correlation (tau positive) or anticorrelation (tau negative).
//p-values:
//http://en.wikipedia.org/wiki/P-value
void kendallN2(float dataA[],
               float dataB[],
               unsigned long length,
               float* tau,
               float* z,
               float* p) {
    
    //NOTE: Book's listing had the wrong indices. Fixed.
    
    unsigned long lastIndex = length - 1;
    unsigned long n2 = 0, n1 = 0, k, j;
    long is = 0;
    float svar;
    float aa,a2,a1;
    
    //Loop over first member of pair, and second member.
    for (j = 0; j < lastIndex; ++j) {
        for (k = (j+1); k <= lastIndex; ++k) {
            a1 = dataA[j] - dataA[k];
            a2 = dataB[j] - dataB[k];
            aa= a1 * a2;
            if (aa) { //Neither array has a tie.
                ++n1;
                ++n2;
                (aa > 0.0f) ? ++is : --is;
            }
            else { //One or both arrays have ties.
                if (a1) {
                    ++n1;
                }
                if (a2) {
                    ++n2;
                }
            }
        }
    }
    
    *tau = is/(sqrt((double) n1)*sqrt((double) n2)); //Equation (14.6.8)
    svar = (4.0f * length + 10.0f) / (9.0f * length * (length - 1.0f)); //Equation (14.6.9)
    *z = (*tau) / sqrt(svar);
    *p = erfcc(fabs(*z) / 1.4142136f);  //Significance
}



//kendallNlogN(NlogN)
//Assumes arr1 already sorted and arr2 already re-ordered in lock-step
//Adapted from https://stat.ethz.ch/pipermail/r-devel/attachments/20100223/077b46e0/attachment.c
//Coded by David Simcha
//From paper: A Computer Method for Calculating Kendall's Tau with Ungrouped Data
// William R. Knight Journal of the American Statistical Association, Vol. 61,
// No. 314, Part 1 (Jun., 1966), pp. 436-439
void kendallNlogN(float dataA[],
                  float dataB[],
                  uif32 length,
                  float* tau,
                  float* z,
                  float* p) {
	
    unsigned long  m1 = 0;
    unsigned long m2 = 0;
	long s;
	uif32 i;
    uif32 swapCount;
	
	uif32 nPair = length * (length - 1) / 2;
	s = nPair;
	
	uif32 tieCount = 0;
	for (i = 1; i < length; ++i) {
		if (dataA[i - 1] == dataA[i]) {
			tieCount++;
		}
		else if (tieCount > 0) {
			insertionSortFloatArray(dataB + i - tieCount - 1, tieCount + 1);
			m1 += tieCount * (tieCount + 1) / 2;
			s += getMs(dataB + i - tieCount - 1, tieCount + 1);
			tieCount++;
			tieCount = 0;
		}
	}
	if (tieCount > 0) {
		insertionSortFloatArray(dataB + i - tieCount - 1, tieCount + 1);
		m1 += tieCount * (tieCount + 1) / 2;
		s += getMs(dataB + i - tieCount - 1, tieCount + 1);
		tieCount++;
	}
	
	swapCount = mergeSortFloatArray(dataB, dataA, length);
	
	m2 = getMs(dataB, length);
	s -= (m1 + m2) + 2 * swapCount;
	
	float denominator1 = nPair - m1;
	float denominator2 = nPair - m2;
    
	*tau = s / sqrt(denominator1) / sqrt(denominator2);
    
	float svar = (4.0f * length + 10.0f) / (9.0f * length * (length - 1.0f)); //Equation (14.6.9), from "Numerical Recipies, The Art of Scientific Computing, Second Edition", p667
	
	*z= (*tau) / sqrt(svar);
	
	*p = erfcc(fabs(*z) / 1.4142136f);  //Significance
}


//Helper function for kendallNlogN()
//Assumes data paramter is sorted
unsigned long getMs(const float data[],
                    const uif32 length) {
	
    uif32 ms = 0;
    uif32 tieCount = 0;
	
	for (uif32 i = 1; i < length; i++) {
		
		if (data[i] == data[i-1]) {
			++tieCount;
		}
		else if (tieCount) {
			ms += (tieCount * (tieCount + 1)) / 2;
			++tieCount;
			tieCount = 0;
		}
	}
	
	if (tieCount) {
		ms += (tieCount * (tieCount + 1)) / 2;
		++tieCount;
	}
	
	return ms;
}



void printHistogram(const float data[],
                    const ui32 length,
                    const uif32 numBins,
                    const uif32 maxNumCharsForBinBar) {
	
	if (length <= 0) {
		return;
	}
	
	//Find maximum/minimum values
	float maxVal = data[0];
	float minVal = data[0];
	for (uif32 i = 0; i < length; ++i) {
		if (maxVal < data[i]) {
			maxVal = data[i];
		}
		if (minVal > data[i]) {
			minVal = data[i];
		}
	}
    if (isnan(minVal)) {
        iLog("Could not print histogram; min_val == NaN");
        return;
    }
	
	//Find the bin size
	float binWidth = ((maxVal - minVal) / numBins);
	if (binWidth == 0) {
        iLog("Could not print histogram; binWidth == 0");
		return;
	}
	
	//Initialize bins
	uif32 bins[numBins];
	for (uif32 i = 0; i < numBins; ++i) {
		bins[i] = 0;
	}

	//Fill bins
	for (uif32 i = 0; i < length; ++i) {
        if (isnan(data[i])) {
            iLog("Could not print histogram; data value == NaN");
            return;
        }
		uif32 binIndex = floor( (data[i] - minVal) / binWidth );
		//Ensure bin_i isn't too high
		if (binIndex > (numBins - 1)) {
			binIndex = numBins - 1;
		}
		++(bins[binIndex]);
	}
	
	//Find the longest bin
	uif32 maxBin = bins[0];
	for (uif32 i = 0; i < numBins; ++i) {
		if (maxBin < bins[i]) {
			maxBin = bins[i];
		}
	}

	
	//Scale, if required
	float binScaleFactor = ((float)maxNumCharsForBinBar) / ((float)maxBin);
	if (binScaleFactor < 1.0f) {
		for (uif32 i = 0; i < numBins; ++i) {
			bins[i] = ceil(bins[i] * binScaleFactor);
		}
	}
	
	//Print
	printf("Histogram:\n");
	for (uif32 i = 0; i < numBins; ++i) {
		printf("\t[%09.3f - %09.3f]\t\t", (i * binWidth) + minVal, ((i+1) * binWidth) + minVal);
		for (uif32 j = 0; j < bins[i]; ++j) {
			printf("*");
		}
		printf("\n");
	}
}
//------------------



//Misc Math
//++++++++++++++++++++
bool oppositeSigns(const float x,
                   const float y) {
    return (x < 0)? (y >= 0): (y < 0);
}



void convolve(const float* dataA,
              const ui32 aLength,
              const float* dataB,
              const ui32 bLength,
              float result[]) { //result_len is aLength + bLength - 1
    
    ui32 kmin, kmax, k;
    
    for (ui32 n = 0; n < aLength + bLength - 1; ++n) {
        
        kmin = (n >= bLength - 1) ? n - (bLength - 1) : 0;
        kmax = (n < aLength - 1) ? n : aLength - 1;
        
        result[n] = 0.0f;
        for (k = kmin; k <= kmax; ++k) {
            result[n] += dataA[k] * dataB[n - k];
        }
    }
}



float diffFraction(const float x,
                   const float y) {
    
    float absDiff = std::abs(x - y);
    float avg = (x + y) / 2;
    
    return (absDiff / avg);
}



float errorFraction(const float observed,
                    const float expected) {
    
    float absDiff = std::abs(observed - expected);
    float result = absDiff / expected;
    //Clipping to keep result [0 1]
    if (result > 1.0) {
        result = 1.0f;
    }
    return result;
}
//--------------------



//Printing
//++++++++++++++++++++

//uint/int vectors
//++++++
void printv(const char name[],
            const uif32 data[],
            const uif32 length) {
    printv(name,
           (sif32*)data,
           1,
           length,
           false);
}



void printv(const char name[],
            const sif32 data[],
            const uif32 length) {
    printv(name,
           data,
           1,
           length,
           false);
}



void printv(const char name[],
            const sif32 data[],
            const uif32 stride,
            const uif32 length,
            const bool asColumn) {
    printf("\n%s:\n", (name) ? name : "Vector:");
    for (uif32 i = 0; i < length; i += stride) {
        printf("%d %s", data[i], (asColumn) ? "\n" : "");
    }
    printf("\n\n");
}
//------



//float vectors
//++++++
void printv(const char name[],
            const float data[],
            const uif32 length) {
    printv(name,
           data,
           1,
           length,
           false);
}



void printv(const char name[],
            const float data[],
            const uif32 length,
            const bool asColumn) {
    printv(name,
           data,
           1,
           length,
           asColumn);
}



void printv(const char name[],
            const float data[],
            const uif32 stride,
            const uif32 length) {
    printv(name,
           data,
           stride,
           length,
           false);
}



void printv(const char name[],
            const float data[],
            const uif32 stride,
            const uif32 length,
            const bool asColumn) {
    printf("\n%s:\n", (name) ? name : "Vector:");
    for (uif32 i = 0; i < length; i += stride) {
        printf("%f %s", data[i], (asColumn) ? "\n" : "");
    }
    printf("\n\n");
}
//------
//--------------------

