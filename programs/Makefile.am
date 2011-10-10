AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags`
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

noinst_PROGRAMS = devtest #distmatrix

devtest_SOURCES = devtest.c
devtest_LDADD = ../src/libeegtools.la

#distmatrix_SOURCES = distmatrix.c
#distmatrix_LDADD = ../src/libeegtools.la
