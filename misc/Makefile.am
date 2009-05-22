AM_CPPFLAGS = -I../src `@GSL_CONFIG@ --cflags` # -I/home/ihrke/local/include
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

bin_PROGRAMS = rawplot
noinst_HEADERS = clusterplot.R paramplot.R paramscan.py channels.R clusterplot.R functions.R rawls.py

rawplot_SOURCES = rawplot.c
rawplot_CFLAGS = @PLOTTER_CFLAGS@
rawplot_LDFLAGS = @PLOTTER_LDFLAGS@
rawplot_LDADD = ../src/libeegtools.la
