matlabmex_SCRIPTS= signaldistmatrix.$(MEXEXT) distancetransform.$(MEXEXT) \
	dtwpath.c
matlabmexdir=$(MATLAB_DIR)
mex_cfiles = $(mexglx_SCRIPTS:.$(MEXEXT)=.c) mex_utils.c mex_utils.h

EXTRA_DIST = $(mex_cfiles)

mexcflags=-I../../src -I$(top_srcdir)/src -I${srcdir}/src -D_MATLAB_VER -g
mexldflags=-L../../src/.libs -leegtools

%.$(MEXEXT): %.c
	$(MEX) $(mexcflags) $^ $(mexldflags)


signaldistmatrix.$(MEXEXT): signaldistmatrix.c mex_utils.c
distancetransform.$(MEXEXT): distancetransform.c mex_utils.c
dtwpath.$(MEXEXT): dtwpath.c mex_utils.c


clean-local:
	rm -f *.$(MEXEXT)

.PHONY: clean
