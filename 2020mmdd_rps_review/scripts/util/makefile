FC = gfortran
# When compatibility with R2014b and older is no longer necessary,
# the -DMEX_DOUBLE_HANDLE argument can be dropped.
#   https://www.mathworks.com/help/matlab/matlab_external/upgrade-mex-files-to-use-graphics-objects.html
MEX ?= mex -DMEX_DOUBLE_HANDLE
export FC
export MEX

ifeq ($(FC),g77)
FLIB = -lg2c
else
FLIB = -lgfortran
endif

.PHONY: clean arcfile libsphere

# The default thing to make: everything
all: \
	arcfile \
	hfillc.mexa64 hfill2c.mexa64 \
	hprofc.mexa64 \
	matminc.mexa64 \
	azel2radecc.mexa64 radec2azelc.mexa64 \
	sunposc.mexa64 moonposc.mexa64 \
	healpix_neighborhood_c.mexa64 \
	calc_legendre_c.mexa64 \
	libsphere

hfillc.mexa64: hfillc.c
	$(MEX) hfillc.c

hfill2c.mexa64: hfill2c.c
	$(MEX) hfill2c.c

hprofc.mexa64: hprofc.c
	$(MEX) hprofc.c

matminc.mexa64: makefile
	make -C minuit

azel2radecc.mexa64: makefile point.h
	$(MEX) CFLAGS='-std=gnu99' CFLAGS='-fPIC' azel2radecc.c

radec2azelc.mexa64: makefile point.h
	$(MEX) CFLAGS='-std=gnu99' CFLAGS='-fPIC' radec2azelc.c

sunposc.mexa64: makefile point.h
	$(MEX) CFLAGS='-std=gnu99' CFLAGS='-fPIC' sunposc.c

moonposc.mexa64: makefile point.h
	$(MEX) CFLAGS='-std=gnu99' CFLAGS='-fPIC' moonposc.c

healpix_neighborhood_c.mexa64: healpix_neighborhood_c.c
	$(MEX) healpix_neighborhood_c.c

calc_legendre_c.mexa64: makefile calc_legendre_c.c
	$(MEX) CFLAGS='$$CFLAGS -Wall' -DLEGENDRE_MAJOR_ORDER=1  calc_legendre_c.c

libsphere:
	make -C libsphere

arcfile:
	make -C arcfile

tab_by_bin.mexa64: makefile tab_by_bin.c
	$(MEX) tab_by_bin.c

clean:
	make -C arcfile clean
	make -C minuit clean
	make -C libsphere clean
	rm *.mex*

