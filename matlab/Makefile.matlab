MEX=mex
CFLAGS=-I../src -g
LDFLAGS=../src/denoising.o ../src/averaging.o ../src/mathadd.o ../src/helper.o #../src/.libs/libeegtools #-leegtools
MEXFLAGS=-v
SRCS := $(wildcard ml_padt*.c)

all: 
	cd ../src; make; cd ../matlab;
	for i in ${SRCS}; do \
		$(MEX) $(CFLAGS) $(MEXFLAGS) $$i $(LDFLAGS); \
	done

$(SRCS):
	$(MEX) $(CFLAGS) $(MEXFLAGS) $< $(LDFLAGS)

.c:
	$(mex) $(CFLAGS) $@.c $(LDFLAGS) -o $@