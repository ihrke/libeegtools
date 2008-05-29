#COMPILE = $(MEX) -v -g $@

#METASOURCES = AUTO
noinst_DATA = artificial_data.m closest.m fBm.m \
			get_monotonous_function.m ml_make.m	plot_spectrum.m rmse.m	snr.m testwa.m warpavgtest.m writesig.m #Makefile.matlab
#	ml_denoise.c ml_dwt.c ml_timewarp.c ml_warpavg.c \
#		 ml_testwarppath.c ml_warptwo.c

#bin_PROGRAMS = script
#script_SOURCES = ml_testwarppath.c
#script_COMPILE = mex -v -g $@
#script_LINK =
SOURCES = Makefile.matlab	ml_denoise.c ml_dwt.c ml_timewarp.c ml_warpavg.c \
		 ml_testwarppath.c ml_warptwo.c
all:
	make -f Makefile.matlab