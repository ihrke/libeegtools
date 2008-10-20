AM_CPPFLAGS = -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags` # -I/home/ihrke/local/include
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

bin_PROGRAMS = rawplot
noinst_HEADERS = clusterplot.R paramplot.R paramscan.py

rawplot_SOURCES = rawplot.c
rawplot_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
rawplot_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
rawplot_LDADD = $(top_builddir)/src/libeegtools.la
