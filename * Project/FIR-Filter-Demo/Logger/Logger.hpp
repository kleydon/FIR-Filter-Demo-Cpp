//Logger.hpp
//
//Cross-platform type-safe, thread-safe logger that plays nice with Crashlytics.
//Optimized for portability. May not be the fastest.
//
//Usage:
//
//Include header: "Logger.hpp" in files you want to log from.
//Note: Can only include from C++ or Objective-C++ files.
//As a result, DON'T USE THIS LOGGER IN AppDelegate!
//
//
//Implementation Notes"
//
//See:
//  http://www.sureshjoshi.com/mobile/cross-platform-mobile-logging-macro/
//  http://stackoverflow.com/questions/10921108/objective-c-preprocessor-definition-dynamic-c-string-to-nsstring-declaration
//
//Note 1: The use of macro functions here is extensive, so that __FILE__,
//__LINE__, etc will reflect the details of the caller, and not the
//logging functions themselves.
//
//Note 2: ##__VA_ARGS__ used instead of __VA_ARGS__ to allow for format-only
//strings, e.g to avoid a trailing comma: printf("Hi",) for some compilers
//which don't automatically omit it.
// GCC and VisualStudio support this. Clang does too, at least for the Mac.

//Note 3: A trick to having multiple statements, local vars, etc, in a macro is
//to put them in the body of a do...while loop.


#ifndef __LOGGER_H__
#define __LOGGER_H__



#ifdef QX_ANDROID

    #include "LoggerAndroid.hpp"

#elif defined(QX_IOS)

    #include "LoggerIOS.h"

#else

    #include "LoggerDefaultCpp.hpp"

#endif

#endif //__LOGGER_H__

