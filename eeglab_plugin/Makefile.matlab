EEGLAB_PATH=/home/ihrke/install/eeglab6.01b/
PROJNAME=warpavg

all:
	rm -rf $(EEGLAB_PATH)/plugins/$(PROJNAME)
	mkdir $(EEGLAB_PATH)/plugins/$(PROJNAME)
	cp *.m $(EEGLAB_PATH)/plugins/$(PROJNAME)
	cp ../matlab/ml* $(EEGLAB_PATH)/plugins/$(PROJNAME)
