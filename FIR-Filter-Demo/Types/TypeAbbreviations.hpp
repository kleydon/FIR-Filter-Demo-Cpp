//TypeAbbreviations.hpp

//Convenience abbreviations for commonly-used standard types

#ifndef __TYPE_ABBREVIATIONS_HPP__
#define __TYPE_ABBREVIATIONS_HPP__

#include <cstddef> //For ptrdiff_t
#include <cstdint>
#include <vector>



//Exact, Signed
typedef int8_t si8;
typedef int16_t si16;
typedef int32_t si32;
typedef int64_t si64;

//Exact, Unsigned
typedef uint8_t ui8;
typedef uint16_t ui16;
typedef uint32_t ui32;
typedef uint64_t ui64;

//Fast, Signed
typedef int_fast8_t sif8;
typedef int_fast16_t sif16;
typedef int_fast32_t sif32;
typedef int_fast64_t sif64;

//Fast, Unsigned
typedef uint_fast8_t uif8;
typedef uint_fast16_t uif16;
typedef uint_fast32_t uif32;
typedef uint_fast64_t uif64;

//Pointer
typedef ptrdiff_t ptrdiff;
typedef intptr_t siptr;
typedef uintptr_t uiptr;

//Vectors
typedef std::vector<si16> si16vec;
typedef std::vector<float> fvec;



#endif //__TYPE_ABBREVIATIONS_HPP__