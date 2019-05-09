# FIR-Filter, C++

Easy-to-use DSP filter code, written in C++, with low, high, and band-pass variants.
Implements a windowed sync FIR using the overlap-and-add method.

## Usage:

```
FIRFilter* filter = &FIRFilter::getSingleton();
    
//Low Pass:
filter->initializeAsLowPass(44100.0f, //sample rate
                            256, //frame length
                            4000.0f, //cutoff frequency
                            1000.0f); //transition bandwidth
    
//High Pass:
filter->initializeAsHighPass(44100.0, //sample rate
                             256, //frame length
                             1000.0f, //cutoff frequency
                             1000.0f); //transition bandwidth
    

//Band Pass:
filter->initializeAsBandPass(44100.0, //sample rate
                             256, //frame length
                             500.0f, //low cutoff frequency
                             200.0f, //low transition bandwidth
                             4000.0f, //high cutoff frequency
                             200.0f); //high transition bandwidth

//Demo/Test:
float frequency = 4500.0f;
float duration = 0.010f;
filter->test(frequency, duration);

//filter->test() calls:
applyFilter(const float inputFrameSamples[],
            float outputFrameSamples[]) 
//...which can be in-place w.r.t. input and output.
```
