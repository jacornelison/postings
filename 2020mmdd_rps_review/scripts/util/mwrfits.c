/*
mwrfits.c 
Edward Wu March 22, 2006, QUaD collaboration
A mex function for matlab that more or less duplicates the functionality
of MWRFITTS in IDL. No support for multiple HDUs or strings yet.

Input: mwrfits(<structure>, <filename>)

The structure should be arrays of floats or integers of various types. This
program should automatically handle their typing provided it is made explicit
in matlab and allocate the FITS file efficiently.

The filename will overwrite an existing file only if you append "!" in
front of the file name. Also, I recommend ending it in ".fits"

Compilation:
Requires the cfitsio library installed on the machine in question. 
Compile with: mex mwrfits.c libcfitsio.a


NOTES:
Internally, the data is saved in a way that follows the ridiculous IDL 
MWRFITS convention of creating a 1 row binary table and then using
the TDIM keywords to specify the dimensionality of the data.  This is
a pain for cfitsio users, but parses automatically for IDL users. This is
how the Stanford TOD data is arranged currently, so I'm assuming 
anyone in the collaboration who needs to handle it can already do so.


*/


#include "fitsio.h"
#include "mex.h"
#define FILENAME_LENGTH 500
#define COLNAME_LENGTH 100

fitsfile * fptr;
int status;    /* FITS errors status reporting var */
int hdutype = BINARY_TBL;

void errorCheck (int i) {
    if (status != 0) {
    char error[1000];
    fits_get_errstatus(status, error);
    mexPrintf("%d: %s\n", i, error);
    mexErrMsgTxt("Something went wrong in CFITSIO."); 
  }
}

/* get an mxArray which presumably has descended from a structure. 
Get the relevant characteristics needed to convert it into a FITS column, like
data type, dimensionality, total size, etc. i is the column index to use */
void parseArray (mxArray *mxa, int i, 
		 char ** ttype,
		 char ** tform, 
		 char ** tunit,
		 char ** tdim, 
		 int * ctype,
		 long * size,
		 void ** data) {

  int dims = mxGetNumberOfDimensions(mxa);
  const int * dim_array = mxGetDimensions(mxa);

  /* get the actual data pointer and transfer it over.  Hope to hell that 
     the column major/row major stuff works, might have to be fixed */
  data[i] = mxGetData(mxa);


  if (dims > 1) {
    /* if we're 2D or more, we have to update tdim properly. First, seed
       tdim and the size variable. */
    sprintf (tdim[i], "(%d", dim_array[0]);

    size[i] = dim_array[0];
    

    int j;
    char buf[100];
    for (j = 1; j < dims; j++) {
      size[i] = size[i] * dim_array[j];
      sprintf (buf, ", %d",dim_array[j]);
      strcat (tdim[i], buf);
    }
    strcat (tdim[i], ")");

  }
  else {
    /* else, toss tdim for this element and set the size to 
       the first dimension's length for just a vector */

    mxFree(tdim[i]);
    tdim[i] = NULL;
    size[i] = dim_array[0];
  }

  /* let's update the rest of the column creation fields now that we know the sizes */


  if (mxIsChar(mxa)) {
    sprintf(tform[i], "%dB", size[i]);
    ctype[i] = TBYTE;
  }

  else if (mxIsSingle(mxa)) {
    sprintf(tform[i], "%dE", size[i]);
    ctype[i] = TFLOAT;
  }
  else if (mxIsDouble(mxa)) {
    sprintf(tform[i], "%dD", size[i]);
    ctype[i] = TDOUBLE;
  } 
  else if (mxIsInt8(mxa)) {
    mexErrMsgTxt("Signed 8-bit integers unsupported.");
    sprintf(tform[i], "%dS", size[i]);
    ctype[i] = TBYTE;
  } 
  else if (mxIsUint8(mxa)) {
    sprintf(tform[i], "%dB", size[i]);
    ctype[i] = TBYTE;
  } 
  else if (mxIsInt16(mxa)) {
    sprintf(tform[i], "%dI", size[i]);
    ctype[i] = TSHORT;
  } 
  else if (mxIsUint16(mxa)) {
    sprintf(tform[i], "%dU", size[i]);
    ctype[i] = TUSHORT;
  } 
  else if (mxIsInt32(mxa)) {
    sprintf(tform[i], "%dJ", size[i]);
    ctype[i] = TINT;
  } 
  else if (mxIsUint32(mxa)) {
    sprintf(tform[i], "%dJ", size[i]);
    ctype[i] = TUINT;
  } 
  else if (mxIsInt64(mxa)) {
    sprintf(tform[i], "%dK", size[i]);
    ctype[i] = TLONG;
  } 
  else if (mxIsUint64(mxa)) {
    sprintf(tform[i], "%dV", size[i]);
    ctype[i] = TULONG;
  } 
  else 
    mexErrMsgTxt("Unsupported type passed - logicals?.");




}

parseStructure (const mxArray * mxa, const char * tablename) {
  char ** ttype;  /* column name - this will correspond to the name in the str */
  char ** tform;  /* FITS column datatype description */
  char ** tunit;  /* no way to probe units yet, leave blank */
  char ** tdim;   /* dimensionality keywords to insert */
  int * ctype;    /* c-type to pass to FITS for any given column */
  long * size;   /* total number of elements in a column */
  void ** data;   /* points to the data for any given FITS column (or incoming structure) */
  
  int ncols;   /* number of columns in the data */
  int i; /* counter variable */

  int hdu;

  char extname[COLNAME_LENGTH]; /* name of the table */
  strcpy (extname, tablename);

  if (!mxIsStruct(mxa))
    mexErrMsgTxt("Non-structure passed to parseStructure - shouldn't happen.");

  if (mxGetNumberOfElements(mxa) > 1)
    mexErrMsgTxt("Only pass this thing one structure, not an array.");
  
  int nfields = mxGetNumberOfFields(mxa);
  ncols = nfields;
  
    /* tons of string memory allocation for table descriptors */

  ttype = (char **)mxCalloc(ncols, sizeof(char *));
  tform = (char **)mxCalloc(ncols, sizeof(char *));
  tunit = (char **)mxCalloc(ncols, sizeof(char *));
  tdim  = (char **)mxCalloc(ncols, sizeof(char *));
  ctype = (int *)mxCalloc(ncols, sizeof(int));
  data  = (void **)mxCalloc(ncols, sizeof(void *));
  size = (long *)mxCalloc(ncols, sizeof(long *));
  
  for (i = 0; i < ncols; i++) {
    /* initialize the column descriptor field pointers */
    ttype[i] = (char *)mxCalloc(COLNAME_LENGTH, sizeof(char *));
    tform[i] = (char *)mxCalloc(COLNAME_LENGTH, sizeof(char *));
    tunit[i] = (char *)mxCalloc(COLNAME_LENGTH, sizeof(char *));
    tdim[i]  = (char *)mxCalloc(COLNAME_LENGTH, sizeof(char *));
    strcpy(tunit[i], "(n/a)");  /* units are still unknown */
  }

  /* all right, let's put in our arrays */

  int j;  /* index over original structure fields */
  int k;  /* index over new table columns */
  k =0;
  mexPrintf("nfields: %d\n", nfields);

  for (j = 0; j < nfields; j++) {
    mxArray * mxaa; /* the current field we're looking at */
    mxaa = (mxArray *) mxGetFieldByNumber(mxa, 0, j);
    if (mxaa == NULL) 
      mexErrMsgTxt("Empty field passed in the structure");

    const char * fieldname;

    fieldname =  mxGetFieldNameByNumber(mxa, j);

    mexPrintf("field %d: %s\n", j, fieldname);
    if (fieldname != NULL) 
      strcpy(ttype[k], fieldname); 
    else
      strcpy(ttype[k], "");
    /* the FITS column name is the field name */
    
    if (mxIsStruct (mxaa)) {
      fits_get_hdu_num(fptr, &hdu);
      parseStructure(mxaa, fieldname);
      fits_movabs_hdu(fptr, hdu, &hdutype, &status);
      ncols--;
    }
    else {
      parseArray(mxaa, k, ttype, tform, tunit, tdim, ctype, size, data);
      k++;
    }
  }

  /* write the actual table */
  fits_create_tbl(fptr, BINARY_TBL, 1, ncols, ttype, tform, tunit, extname, &status);

  for (i = 0; i < ncols; i++) {
    /* insert the dimensionality for IDL if its 2 or higher */
    char tdimkey[100];
    sprintf(tdimkey, "TDIM%d", i + 1);
    
    if (tdim[i] != NULL) 
      fits_write_key(fptr, TSTRING, tdimkey, tdim[i], NULL, &status);
    errorCheck(1);
  }


  /* write out the data to the appropriate places */

  errorCheck(2);

  for (i = 0; i < ncols; i++)
    fits_write_col(fptr, ctype[i], i + 1, 1, 1, size[i], data[i], &status);

}



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) {

  status = 0;

  if (nrhs != 2)
    mexErrMsgTxt("Usage: mwrfits(<structure>, <filename>)");

  if (!mxIsStruct(prhs[0]))
    mexErrMsgTxt("You need to pass a structure in argument 1");
  if (!mxIsChar(prhs[1]) || mxGetM(prhs[1]) != 1)
    mexErrMsgTxt("You need to pass a filename in argument 2");

  char filename[FILENAME_LENGTH];
  
  mxGetString(prhs[1], filename, FILENAME_LENGTH);

  fits_create_file(&fptr, filename, &status);  
  /* overwite by using ! in front of the filename */

  if (status != 0) 
    mexErrMsgTxt("Error opening FITS file for write.");

  parseStructure(prhs[0], "data");
  fits_close_file(fptr, &status);
  

  errorCheck (0);
}
    
    


