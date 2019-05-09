//LoggerAndroid.hpp
//
//See: Logger.hpp

#ifndef __LOGGER_ANDROID_HPP__
#define __LOGGER_ANDROID_HPP__


#include "LoggerSharedFunctions.hpp"



#ifndef QX_ANDROID
    #error LoggerAndroid.hpp included in a non-Android build
#endif



#include <android/log.h>



#define LOG_TAG "QX"



namespace QXLogger {
    
    static std::atomic<bool> inLogLineAndroidWrapper(false);
    template <typename... Args> void logLineAndroidWrapper(int androidLogLevel,
                                                           const char androidLogTag[],
                                                           const char filePathStr[],
                                                           uif32 filePathStrLength,
                                                           uif32 lineNumber,
                                                           const char functionName[],
                                                           const char format[],
                                                           const Args&... args) {
        while(inLogLineAndroidWrapper.load());
        inLogLineAndroidWrapper.store(true);
        
        uif32 maxDateTimeLength = 100;
        char dateTimeStr[maxDateTimeLength];
        QXLogger::constructDateTimeString(dateTimeStr, maxDateTimeLength);
        
        fmt::memory_buffer line;
        //Context
        format_to(line,
                  QXLogger::CONTEXT_FORMAT_STR,
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
        // *** There is certainly a faster way...
        format_to(line, "*");
        *(line.end() - 1) = (char)'\0'; //Add terminator
        
        __android_log_print(androidLogLevel, androidLogTag, line.data());
        
        inLogLineAndroidWrapper.store(false);
    }
    
} //namespace QXLogger



#define logLineAndroidWrapperMacro(androidLogLevel, androidLogTag, format, ...) \
    logLineAndroidWrapper( \
        androidLogLevel, \
        androidLogTag, \
        __FILE__, \
        sizeof(__FILE__) - 1,
        __LINE__, \
        __func__, \
        format, \
        ##__VA_ARGS__)




#ifndef NDEBUG  //If not not debugging, i.e. if debugging...

    #define iLog(...) logLineAndroidWrapperMacro(ANDROID_LOG_INFO, LOG_TAG, format, ##__VA_ARGS__)

#else

    #define iLog(...) (QX_LOGGER_NOOP) //Silenced

#endif


#define eLog(...) logLineAndroidWrapperMacro(ANDROID_LOG_ERROR, LOG_TAG, format, ##__VA_ARGS__)



#endif  //__LOGGER_ANDROID_HPP__
