//main.cpp
//
//FIR Filter Demo
//
//Classes to implement a FIR filter, with low/high/band-pass variants and "output good indication".
//Specifically, its a windowed sinc FIR filter, employing the overlap-and-add method.

#include <iostream>
#include "FIRFilter.hpp"

int main(int argc, const char * argv[]) {
    
    FIRFilter* filter = &FIRFilter::getSingleton();
    
    //Low pass
    filter->initializeAsLowPass(44100.0f, //sampleRate
                                       256, //frameLength
                                       4000.0f, //cutoffFreq
                                       1000.0f); //transBw
    
    
    //High pass
    //initialFilter->initializeAsHighPass(44100.0, //sampleRate
    //                                    256, //frameLength
    //                                    1000.0f, //cutoffFreq
    //                                    1000.0f); //transBw
    
    
    
    //Band pass
    //initialFilter->initializeAsBandPass(44100.0, //sampleRate
    //                                    256, //frameLength
    //                                   500.0f, //lowCutoffFreq
    //                                    200.0f, //lowTransBw
    //                                    4000.0f, //highCutoffFreq
    //                                    200.0f); //highTransBw
    
    float frequency = 4500.0f;
    float duration = 0.010f;
    filter->test(frequency,
                 duration);
    
    
    

    return 0;
}
