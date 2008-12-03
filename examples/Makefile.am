AM_CPPFLAGS = -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags` # -I/home/ihrke/local/include
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

METASOURCES = AUTO

bin_PROGRAMS = t_denoise t_denoisestat \
	 t_readeeglab t_testfilter t_rawreadertest t_clustering t_gapstat \
	t_dtw t_padtw t_adtw_markers t_timefreq


t_denoise_SOURCES = argspec.c argspec.h t_denoise.c
t_denoise_LDADD = $(top_builddir)/src/libeegtools.la 
t_denoisestat_SOURCES = argspec.c argspec.h t_denoisestat.c
t_denoisestat_LDADD = $(top_builddir)/src/libeegtools.la 

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


t_dtw_SOURCES = t_dtw.c
t_dtw_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_dtw_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_dtw_LDADD = $(top_builddir)/src/libeegtools.la


t_padtw_SOURCES = t_padtw.c
t_padtw_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_padtw_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_padtw_LDADD = $(top_builddir)/src/libeegtools.la

t_adtw_markers_SOURCES = t_adtw_markers.c
t_adtw_markers_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_adtw_markers_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_adtw_markers_LDADD = $(top_builddir)/src/libeegtools.la

t_timefreq_SOURCES = t_timefreq.c
t_timefreq_CFLAGS = `@PLOTTER_CONFIG@ --with-capi --cflags`
t_timefreq_LDFLAGS = `@PLOTTER_CONFIG@ --with-capi --ldflags`
t_timefreq_LDADD = $(top_builddir)/src/libeegtools.la
