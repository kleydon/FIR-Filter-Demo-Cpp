//LoggerIOS.mm
//
//See: Logger.hpp

#import "LoggerIOS.h"


#ifndef QX_IOS
    #error LoggerIOS.mm included in a non-ios build
#else
    #include <Crashlytics/Crashlytics.h>
#endif



namespace QXLogger {
    
    void logLineWithCrashlyticsForIOS(const char line[]) {
        //printf("Before logging with crashlytics: %s", line);
        CLS_LOG("%s", line);
        //printf("After logging with crashlytics: %s", line);
    }
    
} //namespace QXLogger
