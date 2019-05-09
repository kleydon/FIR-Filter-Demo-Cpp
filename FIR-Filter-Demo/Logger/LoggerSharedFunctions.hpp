//LoggerSharedFunctions.hpp
//
//See: Logger.hpp

#ifndef __LOGGER_SHARED_FUNCTIONS_HPP__
#define __LOGGER_SHARED_FUNCTIONS_HPP__


#include "TypeAbbreviations.hpp"

#include "format.h"

#include <chrono>
#include <ctime>
#include <inttypes.h>
#include <iostream>


using namespace std::chrono;



//Log Levels:
#define QX_LOG_LEVEL_INFO ("I")
#define QX_LOG_LEVEL_ERROR ("E")

#define QX_LOGGER_NOOP ((void)0)



namespace QXLogger {

    
    
    static const char* CONTEXT_FORMAT_STR = "{} {} {}:{} ({}):\t";
    
    

    #define logLineWithDefaultCppWrapperMacro(stream, level, format, ...) \
        QXLogger::logLineWithDefaultCppWrapper( \
            stream, \
            level, \
            __FILE__, \
            sizeof(__FILE__) - 1, \
            __LINE__, \
            __func__, \
            format, \
            ##__VA_ARGS__)
    
    
    
    inline const char* callersFileName(const char callersFilePath[], uif32 length) {
        for (ptrdiff i = length - 1; i > 0; --i) {
            if (callersFilePath[i] == '/' || callersFilePath[i] == '\\') {
                return &callersFilePath[i+1];
            }
        }
        return callersFilePath;
    }
    
    

    static std::atomic<bool> inConstructDateTimeString(false);
    inline const char* constructDateTimeString(char* dateTimeString,
                                               const uif32 maxLength) {
        
        while(inConstructDateTimeString.load());
        inConstructDateTimeString.store(true);
        
        system_clock::time_point tp = system_clock::now();
        
        milliseconds ms = duration_cast<milliseconds>(tp.time_since_epoch());
        seconds s = duration_cast<seconds>(ms);
        
        std::time_t t = (long)s.count();
        uif32 fractionalSeconds = ms.count() % 1000;
        
        std::strftime(dateTimeString, maxLength, "%Y-%m-%d %H:%M:%S:", localtime(&t));
        
        sprintf(dateTimeString + strlen(dateTimeString), "%03" PRIu32, fractionalSeconds);
        
        inConstructDateTimeString.store(false);
        
        return dateTimeString;
    }
    
    
    
    static std::atomic<bool> inLogLineWithDefaultCppWrapper(false);
    template <typename... Args> void logLineWithDefaultCppWrapper(std::FILE* stream,
                                                                  const char level[],
                                                                  const char filePathStr[],
                                                                  uif32 filePathStrLength,
                                                                  uif32 lineNumber,
                                                                  const char functionName[],
                                                                  const char format[],
                                                                  const Args&... args) {
        
        while(inLogLineWithDefaultCppWrapper.load(std::memory_order_acquire));
        inLogLineWithDefaultCppWrapper.store(true, std::memory_order_release);
        
        uif32 maxDateTimeLength = 100;
        char dateTimeStr[maxDateTimeLength];
        constructDateTimeString(dateTimeStr, maxDateTimeLength);
        
        fmt::memory_buffer line;
        //Context
        format_to(line,
                  CONTEXT_FORMAT_STR,
                  &dateTimeStr[0],
                  level,
                  callersFileName(filePathStr, filePathStrLength),
                  lineNumber,
                  functionName);
        //Message
        format_to(line,
                  format,
                  (args)...);
        
        //Add placeholder to ensure line.data will be large enough to hold terminator...
        //(line.data() is a char array, not a string)
        // *** There is certainly a faster way...
        format_to(line, "*");
        *(line.end() - 1) = (char)'\0'; //Add terminator
        
        printf("%s\n", line.data());
        
        inLogLineWithDefaultCppWrapper.store(false, std::memory_order_release);
    }

    
} //namespace QXLogger



#endif  //__LOGGER_SHARED_FUNCTIONS_HPP__
