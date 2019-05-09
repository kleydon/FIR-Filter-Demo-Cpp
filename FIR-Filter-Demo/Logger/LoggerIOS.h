//LoggerIOS.h
//
//See: Logger.hpp

#ifndef __LOGGER_IOS_H__
#define __LOGGER_IOS_H__


#include "LoggerSharedFunctions.hpp"


#ifndef QX_IOS
    #error LoggerIOS.h included in a non-ios device/simulator build
#endif



namespace QXLogger {

    void logLineWithCrashlyticsForIOS(const char* line);
    
    static std::atomic<bool> inLogLineWithCrashlyticsForIOSWrapper(false);
    template <typename... Args> void logLineWithCrashlyticsForIOSWrapper(const char level[],
                                                                         const char filePathStr[],
                                                                         uif32 filePathStrLength,
                                                                         uif32 lineNumber,
                                                                         const char functionName[],
                                                                         const char format[],
                                                                         const Args&... args) {
        while(inLogLineWithCrashlyticsForIOSWrapper.load());
        inLogLineWithCrashlyticsForIOSWrapper.store(true);
        
        uif32 maxDateTimeLength = 100;
        char dateTimeStr[maxDateTimeLength];
        constructDateTimeString(dateTimeStr, maxDateTimeLength);
        
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
        //(line.data() is a char array, not a string)
        // *** There is certainly a faster way...
        format_to(line, "*");
        *(line.end() - 1) = (char)'\0'; //Add terminator
        
        line.write("\n");
        
        logLineWithCrashlyticsForIOS(line.data());
        
        inLogLineWithCrashlyticsForIOSWrapper.store(false);
    }
    
} //namespace QXLogger



#define logLineWithCrashlyticsForIOSWrapperMacro(level, format, ...) \
    QXLogger::logLineWithCrashlyticsForIOSWrapper( \
        level, \
        __FILE__, \
        sizeof(__FILE__) - 1, \
        __LINE__, \
        __func__, \
        format, \
        ##__VA_ARGS__);



//NOTE: Apple's XCode does not, by default enable NDEBUG for release builds. Instead, it defines
//"DEBUG" for debug builds. Be ready for either below, but complain if NDEBUG is not set.
//(Since standard functionality such as assert() statements rely on it
#ifndef DEBUG
    #ifndef NDEBUG
        #error Error: NDEBUG preprocessor macro not defined for an Apple IOS (XCode?) release build.
    #endif
#endif
#if !defined(NDEBUG) || defined(DEBUG)

    #define iLog(format, ...) logLineWithDefaultCppWrapperMacro(stdout, QX_LOG_LEVEL_INFO, format, ##__VA_ARGS__)

    #define eLog(format, ...) logLineWithDefaultCppWrapperMacro(stderr, QX_LOG_LEVEL_ERROR, format, ##__VA_ARGS__)

#else //Release

    #define iLog(format, ...) (QX_LOGGER_NOOP) //Silenced

    #define eLog(format, ...) logLineWithCrashlyticsForIOSWrapperMacro(QX_LOG_LEVEL_ERROR, format, ##__VA_ARGS__)

#endif



#endif  //__LOGGER_APPLE_H__
