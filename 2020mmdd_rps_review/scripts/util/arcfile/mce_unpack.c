/*
 * RWO 090318 - use mex to speed up unpacking feedback, error,
 *              num flux jumps from MCE data files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#ifndef HAVE_OCTAVE
#  include "matrix.h"
#endif
#include "data_mode.h"

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    int data_mode;
    double filter_gain;
    struct mode_key key;
    void * raw_words;
    mxArray * fb;
    mxArray * aux;
    double * fb_dat;
    double * aux_dat;
    int num_raw_dim;
    int * raw_dim;
    int num_raw_elements;
    int r;
    unsigned long ii;

    if ((nrhs < 2) || (nrhs > 3))
        mexErrMsgTxt ("Usage: [fb aux] = mce_unpack (raw_words, data_mode, [filter_gain])");
    if (nlhs > 2)
        mexErrMsgTxt ("Usage: [fb aux] = mce_unpack (raw_words, data_mode, [filter_gain])");

    if (!mxIsNumeric (prhs[0]))
        mexErrMsgTxt ("raw words must be numeric.");
    if (!mxIsNumeric (prhs[1]) || mxIsEmpty (prhs[1]))
        mexErrMsgTxt ("data mode must be numeric.");
    if ((nrhs == 3) && (!mxIsNumeric(prhs[2]) || mxIsEmpty(prhs[2])))
    {
        mexErrMsgTxt ("filter gain must be numeric.");
        filter_gain = (double) mxGetScalar (prhs[2]);
    }
    else
      filter_gain = DEFAULT_FILTER_GAIN;

    data_mode = (int) mxGetScalar (prhs[1]);
    raw_words = mxGetData (prhs[0]);
    if (NULL == raw_words)
        mexErrMsgTxt ("No raw data words given.");

    if ((nlhs < 2) & (data_mode >= 4))
        mexErrMsgTxt ("Data modes 4 and up are mixed, and require two outputs.");

    r = get_mode_key (data_mode, filter_gain, &key);
    if (r != 0)
    {
        printf ("Unknown data mode %d.\n", data_mode);
        mexErrMsgTxt ("mce_unpack failed.");
    }

    num_raw_dim = mxGetNumberOfDimensions (prhs[0]);
    raw_dim = (int *)mxGetDimensions (prhs[0]);
    num_raw_elements = mxGetNumberOfElements (prhs[0]);
    mxClassID inpt_class;
    inpt_class = mxGetClassID (prhs[0]);

    fb = mxCreateNumericArray (num_raw_dim, raw_dim, mxDOUBLE_CLASS, mxREAL);
    fb_dat = mxGetPr (fb);
    aux = mxCreateNumericArray (num_raw_dim, raw_dim, mxDOUBLE_CLASS, mxREAL);
    aux_dat = mxGetPr (aux);

    /* printf ("Allocated everything.\n"); */

    switch (inpt_class)
    {
        case mxDOUBLE_CLASS :
            switch (data_mode)
            {
                case 5 :
                case 8 :
                case 9 :
                case 10 :
                for (ii=0; ii<num_raw_elements; ii++)
                {
                    fb_dat[ii] = UNPACK_FB ((uint32_t)*(((double *)raw_words) + ii), key);
                    aux_dat[ii] = UNPACK_NUM_FJ ((uint32_t)*(((double *)raw_words) + ii), key);
                }
                break;

                default :
                for (ii=0; ii<num_raw_elements; ii++)
                {
                    fb_dat[ii] = UNPACK_FB ((uint32_t)*(((double *)raw_words) + ii), key);
                    aux_dat[ii] = UNPACK_ERR ((uint32_t)*(((double *)raw_words) + ii), key);
                }
            }
            break;

        case mxUINT32_CLASS :
        case mxINT32_CLASS :
            switch (data_mode)
            {
                case 5 :
                case 8 :
                case 9 :
                case 10 :
                for (ii=0; ii<num_raw_elements; ii++)
                {
                    fb_dat[ii] = UNPACK_FB (*(((uint32_t *)raw_words) + ii), key);
                    aux_dat[ii] = UNPACK_NUM_FJ (*(((uint32_t *)raw_words) + ii), key);
                }
                break;

                default :
                for (ii=0; ii<num_raw_elements; ii++)
                {
                    fb_dat[ii] = UNPACK_FB (*(((uint32_t *)raw_words) + ii), key);
                    aux_dat[ii] = UNPACK_ERR (*(((uint32_t *)raw_words) + ii), key);
                }
            }
            break;

        default :
            mxDestroyArray (fb);
            mxDestroyArray (aux);
            mexErrMsgTxt ("mex_unpack takes data as class double or uint32.");
    }

    if (nlhs < 2)
    {
        if (data_mode == 0)
        {
            plhs[0] = aux;
            mxDestroyArray (fb);
        }
        else
        {
            plhs[0] = fb;
            mxDestroyArray (aux);
        }
    }
    else
    {
        plhs[0] = fb;
        plhs[1] = aux;
    }

    return;
}

