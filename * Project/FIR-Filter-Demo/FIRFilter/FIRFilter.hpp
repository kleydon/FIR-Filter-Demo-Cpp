//FIRFilter.hpp

//A filter for audio samples
//Includes low, high and band-pass variants, and "output good" indicatation.
//
//Specifically, a windowed sinc FIR filter, employing the overlap-and-add method.
//
// See:
//
// * Awsome practical overview of windowed sinc filters:
//      http://www.analog.com/media/en/technical-documentation/dsp-book/dsp_book_Ch16.pdf
//
// * Filter design method and tools from:
//      https://tomroelandts.com/articles/how-to-create-a-simple-low-pass-filter
//      https://tomroelandts.com/tools/filter-designer

// * Helpful text/examples of overlap-and-add method:
//      http://www.comm.utoronto.ca/~dkundur/course_info/real-time-DSP/implementation/Kundur_Lab3_FFT_Convolution_6437.pdf
//      https://www.dsprelated.com/freebooks/sasp/Overlap_Add_OLA_STFT_Processing.html
//      https://www.dsprelated.com/freebooks/sasp/Example_Overlap_Add_Convolution.html
//
// * Helfpul notes on dealing with vDSP implementation, and packing wierdness:
//      http://www.kvraudio.com/forum/viewtopic.php?t=373118
// * Hamilton Kibbe's discussion/code on OLA with vDSP:
//      http://hamiltonkibbe.com/finite-impulse-response-filters-using-apples-accelerate-framework-part-i/
//      http://hamiltonkibbe.com/finite-impulse-response-filters-using-apples-accelerate-framework-part-2/
//      http://hamiltonkibbe.com/finite-impulse-response-filters-using-apples-accelerate-framework-part-iii/
// * Christian Flois' discussion/code on input-side / output-side variations and symmetry-based optimization:
//      https://christianfloisand.wordpress.com/tag/convolution/
// * Iowa Hills' delay-line implementation (first works; second doesnt):
//      http://www.iowahills.com/Example%20Code/FIRFloatingPtImplementation.txt
// * MWickert's ISR on page 7-40 provides a delay-line implementation with rotating indices:
//      http://www.eas.uccs.edu/~mwickert/ece5655/lecture_notes/ece5655_chap7.pdf
// * Helpful note on FFT symmetry properties and real / imaginary results:
//      http://blogs.mathworks.com/steve/2010/07/16/complex-surprises-from-fft/


#ifndef __FIR_FILTER_HPP__
#define __FIR_FILTER_HPP__


#include "FFT.hpp"
#include "TypeAbbreviations.hpp"



class FIRFilter {
    
    public:
    
        static const float VALUE_NOT_AVAILABLE;
    
        static FIRFilter& getSingleton() {
            static FIRFilter singleton;
            return singleton;
        }
    
        bool initializeAsLowPass(const float sampleRate,
                                 const uif32 frameLength, //even
                                 const float cutoffFreq,
                                 const float transBw);
    
        bool initializeAsHighPass(const float sampleRate,
                                  const uif32 frameLength, //even
                                  const float cutoffFreq,
                                  const float transBw);
    
        bool initializeAsBandPass(const float sampleRate,
                                  const uif32 frameLength, //even
                                  const float lowCutoffFreq,
                                  const float initialLowTransBw,
                                  const float highCutoffFreq,
                                  const float initialHighTransBw);
    
        bool configure(const float sampleRate,
                        const uif32 frameLength, //even
                        const float lowCutoffFreq,
                        const float lowTransBw,
                        const float highCutoffFreq,
                        const float highTransBw);
    
        void reset();
    
        void applyFilter(const float inputFrameSamples[],
                         float outputFrameSamples[]);
    
        bool pastInitialTransient();
    
        void test(const float frequency,
                  const float duration); //in seconds
    
    
    private:
    
        FIRFilter();
        ~FIRFilter();
        FIRFilter(FIRFilter const&) = delete;
        void operator = (FIRFilter const&) = delete;

        //Setting up filter coefficients
        //++++++    
        void setFilterCoefficientsAndLength(const float sampleRate,
                                            const float lowCutoff,
                                            const float lowTransBw,
                                            const float highCutoff,
                                            const float highTransBw);
    
        void calcLowPassFilterCoefficients(const float fhaHighDivFs,
                                           float coefficients[],
                                           const uif32 length);
    
        void calcHighPassFilterCoefficients(const float fhaLowDivFs,
                                            float coefficients[],
                                            const uif32 coefficientsLength);
    
        void printFilterCoefficients();
        //------

    
        //Filter variants
        //++++++
        void applyDelayLineFilter(const float inputFrameSamples[],
                                  float outputFrameSamples[]);
    
        void applyInputSideConvolutionFilter(const float inputFrameSamples[],
                                             float outputFrameSamples[]);
    
        void applyOutputSideConvolutionFilter(const float inputFrameSamples[],
                                              float outputFrameSamples[]);
    
        void applyFFTConvolutionFilter(const float inputFrameSamples[],
                                       float outputFrameSamples[]);
        //------

    
        //Data
        //++++++
        //Audio
        float sampleRate;
        float lowCutoffFreq; //Frequency at which signal *starts* to dip - not half-amplitude frequency
        float lowTransBw;
        float highCutoffFreq; //Frequency at which signal *starts* to dip - not half-amplitude frequency
        float highTransBw;
    
        //Lengths
        uif32 frameLength; //Even
        uif32 filterLength; //Odd
        uif32 fftLength; //Even power of 2
        uif32 fftLengthDiv2; //Even power of 2
        uif32 convolutionResultLength; //Even
        uif32 overlapLength; //Overlap across filter calls. Even.
    
        //Main arrays
        float* filterCoefficients;
        float* filterPadded;
        float* signalPadded;
        float* signalPaddedXfilterPadded;
        float* overlapBuffer;

        float* tempResult;
        float* tempResultB;
    
        float filterNyquistValue; //Real-only at nyquist frequency
        FFT* fft;
    
        //Delay-line specific
        uif32 delayLineInsertionIndex;
    
        //Transient
        uif32 initialTransientSampleCount;
        uif32 prevInitialTransientSampleCount;
        //------
};



#endif // __FIR_FILTER_HPP__
