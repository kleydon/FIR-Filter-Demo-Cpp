//UtilityApple.h

#ifndef __APPLE__
    #error UtilityApple.h included from non-Apple build.
#endif

//Convert NSString* references to standard UTF8Strings
#define nsStr2StdStr(nsStr) (std::string([nsStr UTF8String]))