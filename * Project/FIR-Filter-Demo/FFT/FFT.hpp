//FFT.hpp
//
//Computes forward (real->complex) and backward (complex->real)
//real in-place fast fourier transforms.

//  * Complex data stored in interleaved format; [real0,imag0, real1,imag1, ...]
//
//  * The forward transform takes N real numbers as input, and produces
//    N/2+1 complex numbers [0-N/2] as output.
//    (Redundant conjugate pairs are omitted.)
//
//  * The backward transform takes N/2+1 complex numbers [0-N/2] as input,
//    and produces N real numbers as output.
//
//  * After a forward transform, the nyquist and 0 frequency bins are both
//    real-valued (no imaginary components). To save space (and keep to powers
//    of 2), the real nyquist value is stored in the (unused) place of first
//    frequency bin; (i.e. at index: 1). The same technique is used by vDSP.



#ifndef __QX_FFT_HPP__
#define __QX_FFT_HPP__


#include "TypeAbbreviations.hpp"


class FFT {
    
    public:
    
        FFT(const uif32 length);
    
        ~FFT();
    
        void forward(float data[]) const;
    
        void backward(float data[]) const;
    
        void complexForward(float data[]) const;
    
        void complexBackward(float data[]) const;
    
        uif32 getLength() const;
    
        uif32 calcComplexLength() const;
    
        void test();

    
    private:
    
        FFT(FFT const&) = delete;
        void operator = (FFT const&) = delete;
        
        uif32 length;
    
        //Variables specific to Oora's implementation;
        //see notes in OouraFFTsg.cpp
        int* ip;
        float* w;
};



#endif //__QX_FFT_HPP__
