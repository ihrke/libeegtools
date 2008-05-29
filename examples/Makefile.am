AM_CPPFLAGS = -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags` # -I/home/ihrke/local/include
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

bin_PROGRAMS = t_warpavg t_denoise t_denoisestat t_timewarpstat t_loocv \
	t_diffmatrix t_readeeglab t_testfilter t_rawreadertest t_clustering t_gapstat


t_warpavg_SOURCES = argspec.c argspec.h t_warpavg.c
t_warpavg_LDADD = $(top_builddir)/src/libeegtools.la 
t_denoise_SOURCES = argspec.c argspec.h t_denoise.c
t_denoise_LDADD = $(top_builddir)/src/libeegtools.la 
t_denoisestat_SOURCES = argspec.c argspec.h t_denoisestat.c
t_denoisestat_LDADD = $(top_builddir)/src/libeegtools.la 
t_timewarpstat_SOURCES = argspec.c argspec.h t_timewarpstat.c
t_timewarpstat_LDADD = $(top_builddir)/src/libeegtools.la
t_loocv_SOURCES = argspec.c argspec.h t_loocv.c
t_loocv_LDADD = $(top_builddir)/src/libeegtools.la
t_diffmatrix_SOURCES = t_diffmatrix.c argspec.c argspec.h
t_diffmatrix_LDADD = $(top_builddir)/src/libeegtools.la
t_readeeglab_SOURCES = t_readeeglab.c
t_readeeglab_LDADD = $(top_builddir)/src/libeegtools.la

t_rawreadertest_SOURCES = t_rawreadertest.c
t_rawreadertest_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_rawreadertest_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_rawreadertest_LDADD = $(top_builddir)/src/libeegtools.la

t_clustering_SOURCES = t_clustering.c
t_clustering_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_clustering_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_clustering_LDADD = $(top_builddir)/src/libeegtools.la

t_gapstat_SOURCES = t_gapstat.c
t_gapstat_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_gapstat_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_gapstat_LDADD = $(top_builddir)/src/libeegtools.la


t_testfilter_SOURCES = t_testfilter.c 
t_testfilter_CFLAGS =  `@PLOTTER_CONFIG@ --with-capi --cflags`
t_testfilter_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_testfilter_LDADD = $(top_builddir)/src/libeegtools.la 