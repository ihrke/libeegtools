# set the include path found by configure
AM_CPPFLAGS =  $(all_includes) `@GSL_CONFIG@ --cflags` #-I/home/ihrke/local//include

# the library search path.
lib_LTLIBRARIES = libeegtools.la
noinst_HEADERS = 

libeegtools_la_SOURCES = mathadd.c averaging.c averaging.h denoising.c \
	helper.c reader.c clustering.c writer.c tools.c time_frequency.c
include_HEADERS = averaging.h denoising.h helper.h mathadd.h definitions.h reader.h clustering.h writer.h\
	escape_codes.h tools.h time_frequency.h
libeegtools_la_LDFLAGS =  `@GSL_CONFIG@ --libs`
