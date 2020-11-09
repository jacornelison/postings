#ifndef __DATA_MODE_H
#define __DATA_MODE_H

#include <stdlib.h>
#include <stdint.h>

struct mode_key {
    struct {
        double g;
        int s1, s2;
    } fb, err;
    struct {
        int s1, s2;
    } num_fj;
};

#define DEFAULT_FILTER_GAIN 2044

#define UNPACK_INT32(w,s1,s2) (((int32_t)w << s2) >> s1)
#define PACK_INT32(r,s1,s2) ((((uint32_t)((int32_t)r) << s1)) >> s2)

int get_mode_key (int data_mode, double filter_gain, struct mode_key * key);
#define UNPACK_FB(w,k) k.fb.g*(UNPACK_INT32(w,k.fb.s1,k.fb.s2))
#define UNPACK_ERR(w,k) k.err.g*(UNPACK_INT32(w,k.err.s1,k.err.s2))
#define UNPACK_NUM_FJ(w,k) UNPACK_INT32(w,k.num_fj.s1,k.num_fj.s2)

void pack_all (double fb, double err, int num_fj, uint32_t * w, struct mode_key * key);

#endif
