//OouraFFTsg.h
//
//Split-Radix Real<->Complex FFT

#ifndef __OOURA_FFTSG__
#define __OOURA_FFTSG__


#include "TypeAbbreviations.hpp"


//cdft: Complex Discrete Fourier Transform
void cdft(uif32 n,
          sif32 dir,
          float *a,
          int *ip,
          float *w);


//rdft: Real Discrete Fourier Transform
void rdft(uif32 n, //length of data a) to transform (forward), or b) transformed (backwards)
          sif32 dir,
          float a[], //data to transform
          int ip[], //temp storage
          float w[]); //temp storage



#endif //__OOURA_FFTSG__

