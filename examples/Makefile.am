INCLUDES = -I$(top_srcdir)/src -I/home/ihrke/local/include
METASOURCES = AUTO
bin_PROGRAMS = t_warpavg t_denoise t_denoisestat t_timewarpstat t_loocv \
	t_diffmatrix
t_warpavg_SOURCES = argspec.c argspec.h t_warpavg.c
t_warpavg_LDADD = $(top_builddir)/src/libeegtools.la -lm -lgsl -lgslcblas
t_denoise_SOURCES = argspec.c argspec.h t_denoise.c
t_denoise_LDADD = $(top_builddir)/src/libeegtools.la -lm -lgsl -lgslcblas
t_denoisestat_SOURCES = argspec.c argspec.h t_denoisestat.c
t_denoisestat_LDADD = $(top_builddir)/src/libeegtools.la -lm -lgsl -lgslcblas
t_timewarpstat_SOURCES = argspec.c argspec.h t_timewarpstat.c
t_timewarpstat_LDADD = $(top_builddir)/src/libeegtools.la -lm -lgsl -lgslcblas
t_loocv_SOURCES = argspec.c argspec.h t_loocv.c
t_loocv_LDADD = $(top_builddir)/src/libeegtools.la
t_diffmatrix_SOURCES = t_diffmatrix.c argspec.c argspec.h
t_diffmatrix_LDADD = $(top_builddir)/src/libeegtools.la
