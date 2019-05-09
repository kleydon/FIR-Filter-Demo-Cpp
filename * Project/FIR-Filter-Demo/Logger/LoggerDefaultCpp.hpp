//LoggerDefaultCpp.hpp
//
//See: Logger.hpp

#ifndef __LOGGER_DEFAULT_CPP_HPP__
#define __LOGGER_DEFAULT_CPP_HPP__


#if defined(__LOGGER_ANDROID__) || defined(__LOGGER_IOS__)
    #error LoggerDefaultCpp.hpp included from an Android or an IOS build
#endif



#include "LoggerSharedFunctions.hpp"



#ifndef NDEBUG //If not **not** debugging, i.e. if debugging...

    #define iLog(format, ...) logLineWithDefaultCppWrapperMacro(stdout, QX_LOG_LEVEL_INFO, format, ##__VA_ARGS__)

#else

    #define iLog(...) (QX_LOGGER_NOOP) //Silenced

#endif


#define eLog(format, ...) logLineWithDefaultCppWrapperMacro(stderr, QX_LOG_LEVEL_ERROR, format, ##__VA_ARGS__)



#endif  //__LOGGER_DEFAULT_CPP_HPP__
