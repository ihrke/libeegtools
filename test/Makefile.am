METASOURCES = AUTO
AM_CPPFLAGS = -I../src -I$(top_srcdir)/src `@GSL_CONFIG@ --cflags` -DCHECKDATADIR=\"$(datadir)\" 
AM_LDFLAGS = `@GSL_CONFIG@ --libs`

checkextradir = $(datadir)
checkextra_DATA = \
	./data/bindoubles1.dat \
	./data/eeg061206_4_resampled500hz_filtered_DT.set \
	./data/vp19_tt.set

if BUILD_CHECKS
noinst_PROGRAMS = checkall
noinst_HEADERS = check_all.h
checkall_SOURCES = checkfcts.c check_all.c check_checkfcts.c check_denoising.c \
	 check_distances.c check_hmm.c check_list.c check_other.c check_nonlinear.c \
	check_array.c
checkall_LDADD = ../src/libeegtools.la
endif