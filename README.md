# FIR-Filter, C++

Easy-to-use DSP filter code, written in C++, with low, high, and band-pass variants.
Implements a windowed sync FIR using the overlap-and-add method.

## Usage:

```
FIRFilter* filter = &FIRFilter::getSingleton();
    
//Low Pass
filter->initializeAsLowPass(44100.0f, //sampleRate
                            256, //frameLength
                            4000.0f, //cutoffFreq
                            1000.0f); //transBw
    
//..or High Pass
filter->initializeAsHighPass(44100.0, //sampleRate
                             256, //frameLength
                             1000.0f, //cutoffFreq
                             1000.0f); //transBw
    

//...or Band Pass
filter->initializeAsBandPass(44100.0, //sampleRate
                             256, //frameLength
                             500.0f, //lowCutoffFreq
                             200.0f, //lowTransBw
                             4000.0f, //highCutoffFreq
                             200.0f); //highTransBw
    
float frequency = 4500.0f;
float duration = 0.010f;
filter->test(frequency, duration);

//filter->test() calls:
applyFilter(const float inputFrameSamples[],
            float outputFrameSamples[]) 
//...which can be in-place w.r.t. input and output.
```
