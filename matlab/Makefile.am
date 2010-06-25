METASOURCES = AUTO

if MATLAB_SUPPORT
  MEX_DIR = mex
endif

SUBDIRS = $(MEX_DIR)

m_SCRIPTS = artificial_data.m artificial_data2.m closest.m eegset_to_raw.m get_epoch_indices.m \
	get_monotonous_function.m ldist.m  raw_to_eegset.m read_raw.m rmse.m select_trials_from_eegset.m \
	snr.m write_art_data_raw.m write_raw.m ml_template.c
mdir = $(MATLAB_DIR)

EXTRA_DIST = $(m_SCRIPTS)