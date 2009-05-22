METASOURCES = AUTO
AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags`
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

if BUILD_CHECKS
noinst_PROGRAMS = checkall
noinst_HEADERS = check_all.h
checkall_SOURCES = checkfcts.c check_all.c check_checkfcts.c check_denoising.c \
	 check_distances.c check_hmm.c
checkall_LDADD = ../src/libeegtools.la
endif