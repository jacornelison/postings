#include <stdlib.h>
#include <stdio.h>
#include "data_mode.h"

static void set_fb (struct mode_key * key, double g, int s1, int s2)
{
    key->fb.g = g;
    key->fb.s1 = s1;
    key->fb.s2 = s2;
}

static void set_err (struct mode_key * key, double g, int s1, int s2)
{
    key->err.g = g;
    key->err.s1 = s1;
    key->err.s2 = s2;
}

static void set_num_fj (struct mode_key * key, int s1, int s2)
{
    key->num_fj.s1 = s1;
    key->num_fj.s2 = s2;
}

int get_mode_key (int data_mode, double filter_gain, struct mode_key * key)
{
    set_fb (key, 0, 0, 0);
    set_err (key, 0, 0, 0);
    set_num_fj (key, 0, 0);
    switch (data_mode)
    {
      case 0 : set_err (key, 1, 0, 32); break;
      case 1 : set_fb (key, 1.0/(double)(1<<12), 0, 32); break;
      case 2 : set_fb (key, 1.0/(double)filter_gain, 0, 32); break;
      case 4 : set_fb (key, 1, 14, 18); set_err (key, 1, 0, 14); break;
      case 9 : set_fb (key, 2.0/(double)filter_gain, 8, 24); set_num_fj (key, 0, 8); break;
      case 10 : set_fb (key, 8.0/(double)filter_gain, 7, 25); set_num_fj (key, 0, 7); break;
      default : return -1;
    }
    key->fb.s2 = 32 - key->fb.s2 - key->fb.s1;
    key->fb.s1 = key->fb.s1 + key->fb.s2;
    key->err.s2 = 32 - key->err.s2 - key->err.s1;
    key->err.s1 = key->err.s1 + key->err.s2;
    key->num_fj.s2 = 32 - key->num_fj.s2 - key->num_fj.s1;
    key->num_fj.s1 = key->num_fj.s1 + key->num_fj.s2;

    return 0;
}

void pack_all (double fb, double err, int num_fj, uint32_t * w, struct mode_key * key)
{
    uint32_t tmp;

    tmp = (uint32_t)(fb / key->fb.g);
    *w = PACK_INT32(tmp,key->fb.s1,key->fb.s2);
    tmp = (uint32_t)(err / key->err.g);
    *w |= PACK_INT32(tmp,key->err.s1,key->err.s2);
    *w |= PACK_INT32(num_fj,key->num_fj.s1,key->num_fj.s2);
}

