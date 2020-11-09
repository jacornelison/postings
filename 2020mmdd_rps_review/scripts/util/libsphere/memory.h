#ifndef MEMORY_H
#define MEMORY_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h> /* size_t */

/**
 * malloc() helper wrapper which adds enhancements.
 */
void* xmalloc(size_t nmemb, size_t size, char* msg);
/**
 * realloc() helper wrapper which adds enhancements.
 */
void* xrealloc(void* ptr, size_t nmemb, size_t size, char* msg);
/**
 * free() helper wrapper which adds enhancements.
 */
void  xfree(void** ptr);

#define MKALLOC(PREFIX, TYPE)           MKALLOC_(, PREFIX, TYPE)
#define MKALLOC_PRIVATE(PREFIX, TYPE)   MKALLOC_(static, PREFIX, TYPE)

#define MKALLOC_(STATIC, PREFIX, TYPE) \
\
/*****************************************************************************\
 * ?_malloc() helper function                                                *\
 ****************************************************************************/\
STATIC TYPE* PREFIX ## _malloc(size_t n, char* msg) \
{\
    return (TYPE*)xmalloc(n, sizeof(TYPE), msg); \
}\
\
/*****************************************************************************\
 * ?_realloc() helper function                                               *\
 ****************************************************************************/\
STATIC TYPE* PREFIX ## _realloc(TYPE* ptr, size_t n, char* msg) \
{\
    return (TYPE*)xrealloc((void*)ptr, n, sizeof(TYPE), msg); \
}\
\
/*****************************************************************************\
 * ?_free() helper function                                                  *\
 ****************************************************************************/\
STATIC void PREFIX ## _free(TYPE** ptr) \
{\
    xfree((void**)ptr); \
}\
\
STATIC void PREFIX ## _set(TYPE* ptr, size_t n, TYPE copy) \
{\
    for (size_t ii=0; ii<n; ++ii)\
        ptr[ii] = copy;\
}\


#define MKALLOC_PROTO(PREFIX, TYPE) \
TYPE* PREFIX ## _malloc(size_t n, char* msg); \
TYPE* PREFIX ## _realloc(TYPE* ptr, size_t n, char* msg); \
void  PREFIX ## _free(TYPE** ptr); \
void  PREFIX ## _set(TYPE* ptr, size_t n, TYPE copy); \


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* MEMORY_H */
