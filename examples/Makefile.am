AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags`
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO


if EXPERIMENTAL_SUPPORT
	experimental_progs=nnsearch
	nnsearch_SOURCES = nnsearch.c
	nnsearch_LDADD = $(top_builddir)/src/libeegtools.la
endif

noinst_PROGRAMS =  array $(experimental_progs)


array_SOURCES = array.c
array_LDADD = $(top_builddir)/src/libeegtools.la
