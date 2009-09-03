# set the include path found by configure
AM_CPPFLAGS =  -DT_LINUX -Isrc $(all_includes) -fPIC --std=gnu99 -Wall @MATIO_CFLAGS@ @DEBUG_FLAGS@ `@GSL_CONFIG@ --cflags`


# the library search path.
lib_LTLIBRARIES = libeegtools.la
noinst_HEADERS = 

libeegtools_la_SOURCES = mathadd.c averaging.c averaging.h denoising.c \
	helper.c reader.c clustering.c writer.c tools.c time_frequency.c warping.c \
	nonlinear.c recurrence_plot.c regularization.c distances.c hmm.c fidlib/fidlib.c \
	som.c eeg.c
pkginclude_HEADERS = averaging.h denoising.h helper.h mathadd.h \
	definitions.h reader.h clustering.h writer.h\
	escape_codes.h tools.h time_frequency.h chaninfo.h warping.h \
	nonlinear.h recurrence_plot.h regularization.h distances.h hmm.h \
	fidlib/fidlib.h fidlib/fidmkf.h fidlib/fidrf_combined.h fidlib/fidrf_cmdlist.h \
	fidlib/fidrf_jit.h som.h eeg.h

libeegtools_la_LDFLAGS =  `@GSL_CONFIG@ --libs` @MATIO_LIBS@

## own rules are ok
tag:
	$(shell CURDIR=`pwd`)
	cd $(top_srcdir)/src; find . -name "*.[chCH][chCH]" -print | etags -; cd $(CURDIR);
