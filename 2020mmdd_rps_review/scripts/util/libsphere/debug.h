/**
 * Configuration options:
 *
 * DEBUG_NOFUNC
 *     If defined, then log messages will not include the function name in
 *     the message output.
 *
 * DEBUG_NOCOLOR
 *     If defined, then terminal color control characters will not be emitted.
 */

#ifndef DEBUG_H
#define DEBUG_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
    #include <stdio.h>

    #define GCC_VERSION (__GNUC__ * 10000 \
                         + __GNUC_MINOR__* 100 \
                         + __GNUC_PATCHLEVEL__)

    // Get the correct function name to be used
    #if GCC_VERSION > 30200
        #define DBGLOG_FUNC __PRETTY_FUNCTION__
    #else
        #define DBGLOG_FUNC __func__
    #endif

    #ifdef DEBUG_NOCOLOR
        #define COLOR_BEGIN ""
        #define COLOR_END   ""
    #else /* DEBUG_NOCOLOR */
        #define COLOR_BEGIN "\x1b[31m"
        #define COLOR_END   "\x1b[0m"
    #endif

    // If function annotations are disabled,
    #ifdef DEBUG_NOFUNC

        #define xdbglog(cb, ce, msg, file, line, ...) \
            fprintf(stderr, cb "(%s:%d): " msg ce, \
                    file, line, ##__VA_ARGS__)
        #define dbglog(msg, ...) \
            xdbglog(COLOR_BEGIN, COLOR_END, msg, __FILE__, \
                    __LINE__, ##__VA_ARGS__)

    #else /* DEBUG_NOFUNC */

        #define xdbglog(cb, ce, msg, func, file, line, ...) \
            fprintf(stderr, cb "%s (%s:%d): " msg ce, \
                    func, file, line, ##__VA_ARGS__)
        #define dbglog(msg, ...) \
            xdbglog(COLOR_BEGIN, COLOR_END, msg, DBGLOG_FUNC, __FILE__, \
                    __LINE__, ##__VA_ARGS__)

    #endif

#else /* NDEBUG */

    /* Give an empty definition of all macros */
    #define dbglog(msg, ...)

#endif /* NDEBUG */

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* DEBUG_H */
