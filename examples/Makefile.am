AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags`
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

noinst_PROGRAMS = nnsearch array

nnsearch_SOURCES = nnsearch.c
nnsearch_LDADD = $(top_builddir)/src/libeegtools.la

array_SOURCES = array.c
array_LDADD = $(top_builddir)/src/libeegtools.la
