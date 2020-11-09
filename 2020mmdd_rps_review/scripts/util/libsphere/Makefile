CC ?= gcc
CFLAGS ?= -g -O3 -ffast-math -Wall
DEBUG ?= 0

override CFLAGS += -std=gnu11 -fPIC -fno-exceptions
ifeq ($(DEBUG),0)
	override CFLAGS += -DNDEBUG
endif

.PHONY: clean
targets: ../libsphere.so
%.o: %.c
	$(CC) -c $(CFLAGS) -o $@ $^
%.so:
	$(CC) -shared $(CFLAGS) $(LDFLAGS) -o $@ $^

../libsphere.so: override LDFLAGS += -lm
../libsphere.so: libsphere.o legendre.o memory.o

clean:
	rm -f *.o *.so ../libsphere.so
	rm -f *_proto.m
