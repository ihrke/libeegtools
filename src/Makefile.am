# set the include path found by configure
AM_CPPFLAGS =  -I$(build_dir)/src $(all_includes) --std=gnu99 -Wall @DEBUG_FLAGS@ `@GSL_CONFIG@ --cflags`

# the library search path.
lib_LTLIBRARIES = libeegtools.la
noinst_HEADERS = 

libeegtools_la_SOURCES = mathadd.c averaging.c averaging.h denoising.c \
	helper.c reader.c clustering.c writer.c tools.c time_frequency.c warping.c
pkginclude_HEADERS = averaging.h denoising.h helper.h mathadd.h definitions.h reader.h clustering.h writer.h\
	escape_codes.h tools.h time_frequency.h chaninfo.h warping.h
libeegtools_la_LDFLAGS =  `@GSL_CONFIG@ --libs`
