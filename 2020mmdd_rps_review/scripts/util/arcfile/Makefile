MEX ?= mex

CFILES = databuf.c dataset.c arcfile.c endian.c fileset.c namelist.c \
		 readarc.c reglist.c utcrange.c handlesig.c
HFILES = $(patsubst %.c,%.h,$(CFILES))
OBJS = $(patsubst %.c,%.o,$(CFILES))

%.o: %.c %.h
	gcc -c -fPIC -std=gnu99 -o $@ $<

all: ../readarc.mexa64 ../listarc.mexa64 ../mce_unpack.mexa64

../readarc.mexa64: $(OBJS)
	$(MEX) LDFLAGS='$$LDFLAGS -lz -lbz2' -outdir .. -output readarc \
		mex_readarc.c $(OBJS)

../listarc.mexa64: $(OBJS)
	$(MEX) LDFLAGS='$$LDFLAGS -lz -lbz2' -outdir .. -output listarc \
		mex_listarc.c $(OBJS)

../mce_unpack.mexa64: mce_unpack.c data_mode.c data_mode.h
	$(MEX) -outdir .. -output mce_unpack mce_unpack.c data_mode.c
