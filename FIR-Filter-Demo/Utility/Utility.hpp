//Utility.hpp

//Miscellaneous convenience functions.

#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__


#include "TypeAbbreviations.hpp"

#include <chrono>
#include <thread>



//Clock/Timer-related
//++++++++++++++++++++
//SteadyClock is monotonic
typedef std::chrono::steady_clock ChronoSteadyClock; //Monotonic
typedef std::chrono::steady_clock::time_point ChronoSteadyClockTimePoint;
typedef std::chrono::duration<float, std::milli> ChronoDurationFloatMs;

template <typename T, typename U>
auto ChronoDurationCast(U const& u) -> decltype(std::chrono::duration_cast<T>(u)) {
    return std::chrono::duration_cast<T>(u);
}
//--------------------


//Float infinity values - for finding min/max values
//++++++++++++++++++++
static const float QX_FLOAT_NEG_INF = std::numeric_limits<float>::lowest();
static const float QX_FLOAT_POS_INF = std::numeric_limits<float>::max();
//--------------------



//Timing-related
//++++++++++++++++++++
void waitMs(uif32 durationMs);

void sleepMs(uif32 durationMs);

void operationTimingTest();
//--------------------



//Memory-related
//++++++++++++++++++++
//Based on: http://stackoverflow.com/a/1265681
//
//NOTE: Use of these deletion functions depend on pointers being initialized to nullptr
//when first created; otherwise, there is the possibility of attempting to delete things
//that were never initially created.

template<class T> void safeDelete(T*& pVal) {
    if (pVal != nullptr) {
        delete pVal;
        pVal = nullptr;
    }
}

template<class T> void safeDeleteArray(T*& pVal) {
    if (pVal != nullptr) {
        delete [] pVal;
        pVal = nullptr;
    }
}
//--------------------




//Waveform-related
//++++++++++++++++++++
void createCosineWave(const float frequency,
                      const float sampleRate,
                      float result[],
                      const uif32 resultLength);
//--------------------




//Power/Log
//++++++++++++++++++++
float sqr(const float x);

bool isPowerOf2(const uif32 x);

uif32 nextPow2(uif32 x);

uif32 logBase2(const uif32 x);
//--------------------



//Interpolation / Mapping
//++++++++++++++++++++
void parabolicLogMaxInterp(const float data[],
                           const uif32 i1,
                           const uif32 length,
                           float* interpolatedMaximum,
                           float* interpolatedMaximumFractionalIndex);

void parabolicMaxInterp(const float data[],
                        const uif32 i1,
                        const uif32 length,
                        float* interpolatedMaximum,
                        float* interpolatedMaximumFractionalIndex);

void parabolicMinInterp(const float data[],
                        const uif32 i1,
                        const uif32 length,
                        float* interpolatedMinimum,
                        float* interpolatedMinimumFractionalIndex);

float linearInterp(const float data[],
                   const uif32 length,
                   const float fractionalIndex);

float sinc(const float x);

float sincInterp(const float data[],
                 uif32 dataLength,
                 float fractionalIndex);
//--------------------



//Resampling / Remapping
//++++++++++++++++++++
void resampleSinc(const float input[],
                  uif32 inputLength,
                  float output[],
                  uif32 outputLength);

void resampleSincTest();

float map(const float inValue,
          const float inRangeMinimum,
          const float inRangeMaximum,
          const float outRangeMinimum,
          const float outRangeMaximum);
//--------------------



//Sorting
//++++++++++++++++++++
void quickSortFloatArray(float data[],
                         const uif64 length);
void quickSortUIF32Array(uif32 data[],
                         const uif64 length);

uif32 insertionSortFloatArray(float data[],
                              const uif32 length);

void insertionSortIndicesOfFloatArray(const float data[],
                                      uif32 indices[],
                                      const uif32 indicesLength);

uif32 mergeSortFloatArray(float data[],
                          float buffer[],
                          const uif32 length); //Sorts in place
//--------------------


//Filtering
//++++++++++++++++++++
void medianFilter(const float inputData[],
                  float outputData[],
                  const uif32 dataLength,
                  const uif32 windowLength);
//--------------------


//Comparison
//++++++++++++++++++++
bool approxEqualAbsTolerance(const float x,
                             const float y);

bool approxEqualRelTolerance(const float x,
                             const float y);

bool approxEqualWithinDiffFraction(const float x,
                                   const float y,
                                   const float diffFraction); //factor: [0.0 1.0]
//--------------------



//Stats
//++++++++++++++++++++
float fastMedian(float data[],
                 const uif32 n);

float getMedianWithoutPresevingInputArray(float data[],
                                          const uif32 length);
float getMedianWhilePresevingInputArray(const float data[],
                                        const uif32 length);

bool linearRegression(const float dataX[],
                      const float dataY[],
                      uif32 length,
                      float* m, //Slope
                      float* b, //Y intercept
                      float* r); //Output correlation coefficient

float erfcc(const float x);

void kendallN2(float dataA[],
               float dataB[],
               uif32 length,
               float* tau,
               float* z,
               float* p);

void kendallNlogN(float dataA[],
                  float dataB[],
                  uif32 length,
                  float* tau,
                  float* z,
                  float* p);

unsigned long  getMs(const float data[],
                     const uif32 length);

void printHistogram(const float data[],
                    const uif32 length,
                    const uif32 numBins,
                    const uif32 maxNumCharsForBinBar);
//--------------------



//Misc Math
//++++++++++++++++++++
bool oppositeSigns(const float x,
                   const float y);

void convolve(const float dataA[],
              const uif32 aLength,
              const float dataB[],
              const uif32 bLength,
              float* result); //resultLength = aLength + bLength - 1

float diffFraction(const float x, //effectively, percent difference
                   const float y);

float errorFraction(const float observed,
                    const float expected);
//--------------------



//Printing
//++++++++++++++++++++

//uint/sint vector
//++++++
void printv(const char name[],
            const uif32 data[],
            const uif32 length);

void printv(const char name[],
            const sif32 data[],
            const uif32 length);

void printv(const char name[],
            const sif32 data[],
            const uif32 stride,
            const uif32 length,
            const bool asColumn);
//------

//float vector
//++++++
void printv(const char name[],
            const float data[],
            const uif32 length);

void printv(const char name[],
            const float data[],
            const uif32 length,
            const bool asColumn);

void printv(const char name[],
            const float data[],
            const uif32 stride,
            const uif32 length);

void printv(const char* name,
            const float data[],
            const uif32 stride,
            const uif32 length,
            const bool asColumn);
//------

//--------------------



#endif //__UTILITY_HPP__
