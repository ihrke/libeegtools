METASOURCES = AUTO
AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags` -DCHECKDATADIR=\"$(datadir)/tests/\"
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

checkextradir = $(datadir)/tests/
checkextra_DATA = \
	./data/bindoubles1.dat \
	./data/eeg061206_4_resampled500hz_filtered_DT.set \
	./data/vp19_tt.set

if BUILD_CHECKS
noinst_PROGRAMS = checkall
noinst_HEADERS = check_all.h

if EXPERIMENTAL_SUPPORT
experimental_tests=check_hmm.c 
endif

checkall_SOURCES = checkfcts.c check_all.c check_checkfcts.c check_denoising.c \
	 check_distances.c check_list.c check_other.c check_nonlinear.c \
	check_array.c check_linalg.c check_warping.c check_cluster.c check_average.c \
	check_io.c check_optarg.c check_eeg.c $(experimental_tests)

checkall_LDADD = ../src/libeegtools.la
endif
