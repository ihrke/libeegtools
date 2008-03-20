
# set the include path found by configure
INCLUDES = -I/home/ihrke/local//include $(all_includes)

# the library search path.
lib_LTLIBRARIES = libeegtools.la
noinst_HEADERS = 

libeegtools_la_SOURCES = mathadd.c averaging.c averaging.h denoising.c \
	denoising.h helper.c helper.h mathadd.h definitions.h
include_HEADERS = averaging.h denoising.h helper.h mathadd.h definitions.h
libeegtools_la_LIBADD = -lm -lgsl -lgslcblas
