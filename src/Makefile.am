# set the include path found by configure
AM_CPPFLAGS =  -Isrc $(all_includes) -D_GNU_SOURCE -fPIC --std=gnu99 -Wall @MATIO_CFLAGS@ @DEBUG_FLAGS@ `@GSL_CONFIG@ --cflags`

# the library search path.
lib_LTLIBRARIES = libeegtools.la
noinst_HEADERS = 


if EXPERIMENTAL_SUPPORT
experimental_sources=hmm.c nnsearch.c isomap.c gapstat.c som.c
experimental_headers=hmm.h som.h nnsearch.h isomap.h gapstat.h
else
experimental_sources=
experimental_headers=
endif

stable_sources=mathadd.c averaging.c averaging.h filter.c wavelet.c \
	helper.c io.c io_wav.c io_matlab.c io_eegtools.c clustering.c tools.c time_frequency.c warping.c \
	nonlinear.c recurrence_plot.c regularization.c distances.c fidlib/fidlib.c \
	eeg.c optarg.c mpitools.c pqueue.c array.c linalg.c \
	complex.c imageproc.c  slist.c

stable_headers= averaging.h filter.h wavelet.h helper.h mathadd.h \
	definitions.h clustering.h escape_codes.h tools.h time_frequency.h chaninfo.h warping.h \
	nonlinear.h recurrence_plot.h regularization.h distances.h \
	fidlib/fidlib.h fidlib/fidmkf.h fidlib/fidrf_combined.h fidlib/fidrf_cmdlist.h \
	eeg.h optarg.h mpitools.h \
	pqueue.h array.h linalg.h complex.h imageproc.h io.h \
	io_matlab.h io_wav.h io_eegtools.h slist.h

libeegtools_la_SOURCES = $(stable_sources) $(experimental_sources)

pkginclude_HEADERS = $(stable_headers) $(experimental_headers)

libeegtools_la_LDFLAGS =  `@GSL_CONFIG@ --libs` @MATIO_LIBS@

## own rules are ok
tag:
	$(shell CURDIR=`pwd`)
	cd $(top_srcdir)/src; find . -name "*.[chCH][chCH]" -print | etags -; cd $(CURDIR);
