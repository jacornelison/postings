#include <signal.h>     /* for raise() */
#include <stdio.h>      /* for fprintf() */
#include <stdlib.h>     /* for malloc() */
#include <string.h>     /* for memset() */
#include <stdint.h>     /* for SIZE_MAX */

#pragma GCC visibility push(hidden)

/* Want s1*s2 <= SIZE_MAX. We can avoid a more costly division checks if we
 * first can note that both {s1,s2} < sqrt(SIZE_MAX) guarantee s1*s2<=SIZE_MAX.
 * Let N = sizeof(size_t).
 *
 * sqrt(SIZE_MAX) = sqrt(2^8N) = (2^8N)^0.5 = 2^4N = (1 << 4N)
 */
#define XMALLOC_SQRT_SIZEMAX (1UL << (sizeof(size_t)*4))

/**
 * malloc() helper wrapper which adds enhancements.
 *
 * 1) Protects against integer overflow when size_t(TYPE)*N > SIZE_MAX
 *    -- Inspired by BSD's reallocarray() implementation.
 * 2) Makes number of elements and size of type explicit
 * 4) Can emit a signal if a memory allocation fails
 * 5) Automatically zeros the requested memory
 */
void* xmalloc(size_t nmemb, size_t size, char* msg)
{
    void*  ptr = NULL;
    size_t nbytes = 0;

    // Check that the byte allocation can't overflow. First check in a
    // fast-path that can guarantee no overflow when both nmemb and size are
    // less than sqrt(SIZE_MAX). Then if either is greater, do a more costly
    // division check.
    if ((nmemb >= XMALLOC_SQRT_SIZEMAX || size >= XMALLOC_SQRT_SIZEMAX)
        && (nmemb > 0 && SIZE_MAX/nmemb < size))
    {
        fprintf(stderr, "****Memory allocation request size overflows: %s\n",
            msg);
        raise(SIGABRT);
    }
    nbytes = nmemb * size;

    ptr = malloc(nbytes);

    // If there's an error, print an error message, and then raise an abort
    if (ptr == NULL)
    {
        fprintf(stderr, "****Memory allocation failed for %zu bytes: %s\n",
            nbytes, msg);
        raise(SIGABRT);
    }

    // Clear the memory
    memset(ptr, 0, nbytes);

    return ptr;
}

/**
 * realloc() helper wrapper which adds enhancements.
 *
 * 1) Protects against integer overflow when size_t(TYPE)*N > SIZE_MAX
 *    -- Inspired by BSD's reallocarray() implementation.
 * 2) Can emit a signal if a memory allocation fails
 */
void* xrealloc(void* ptr, size_t nmemb, size_t size, char* msg)
{
    void*  newptr = NULL;
    size_t nbytes = 0;

    // Check that the byte allocation can't overflow. First check in a
    // fast-path that can guarantee no overflow when both nmemb and size are
    // less than sqrt(SIZE_MAX). Then if either is greater, do a more costly
    // division check.
    if ((nmemb >= XMALLOC_SQRT_SIZEMAX || size >= XMALLOC_SQRT_SIZEMAX)
        && (nmemb > 0 && SIZE_MAX/nmemb < size))
    {
        fprintf(stderr, "****Memory reallocation request size overflows\n");
        raise(SIGABRT);
    }
    nbytes = nmemb * size;

    newptr = realloc(ptr, nbytes);

    // If there's an error, print an error message, and then raise an abort
    if (newptr == NULL)
    {
        fprintf(stderr, "****Memory reallocation failed for %zu bytes: %s",
            nbytes, msg);
        raise(SIGABRT);
    }

    return newptr;
}

/**
 * free() helper wrapper which adds enhancements.
 *
 * 1) Nulls the pointer to zero to avoid double free
 */
void xfree(void** ptr)
{
    if (*ptr != NULL)
    {
        free(*ptr);
        *ptr = NULL;
    }
}

#pragma GCC visibility pop

