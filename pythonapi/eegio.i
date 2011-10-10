%{
#include <libeegtools/reader.h>
#include <libeegtools/writer.h>
#include <libeegtools/helper.h>
%}

EEGdata_trials* read_eegtrials_from_raw(const char *file);
void write_eegtrials_to_raw_file( const EEGdata_trials *eeg, const char *file );
