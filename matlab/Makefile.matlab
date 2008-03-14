all:
	cd ../src; make; cd ../matlab;
#	mex -v -I../src ml_denoise.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o
#	mex -v -I../src ml_timewarp.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o
#	mex -v -I../src ml_testwarppath.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o
#	mex -v -I../src ml_warptwo.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o
#	mex -v -I../src ml_warpavg.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o
	mex -v -I../src ml_running_median.c ../src/helper.o ../src/mathadd.o ../src/denoising.o ../src/averaging.o