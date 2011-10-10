matlabmex_SCRIPTS=signaldistmatrix.$(MEXEXT) distancetransform.$(MEXEXT) dtwpath.$(MEXEXT) crossrecplot.$(MEXEXT) dtwadd.$(MEXEXT) loscalc.$(MEXEXT) distmatrix.$(MEXEXT) tdelayreconstruction.$(MEXEXT) bresenham.$(MEXEXT)

matlabmexdir=$(MATLAB_DIR)
mex_cfiles = $(matlabmex_SCRIPTS:.$(MEXEXT)=.c) mex_utils.c mex_utils.h



mexcflags=-I../../src -I$(top_srcdir)/src -I${srcdir}/src -D_MATLAB_VER -g
mexldflags=-L../../src/.libs -leegtools

%.$(MEXEXT): %.c
	$(MEX) $(mexcflags) $^ $(mexldflags)


signaldistmatrix.$(MEXEXT): signaldistmatrix.c mex_utils.c
distmatrix.$(MEXEXT): distmatrix.c mex_utils.c
distancetransform.$(MEXEXT): distancetransform.c mex_utils.c
dtwpath.$(MEXEXT): dtwpath.c mex_utils.c
dtwadd.$(MEXEXT): dtwadd.c mex_utils.c
crossrecplot.$(MEXEXT): crossrecplot.c mex_utils.c
loscalc.$(MEXEXT): loscalc.c mex_utils.c
tdelayreconstruction.$(MEXEXT): tdelayreconstruction.c mex_utils.c
bresenham.$(MEXEXT): bresenham.c mex_utils.c

m_SCRIPTS=$(matlabmex_SCRIPTS:.$(MEXEXT)=.m)
mdir = $(MATLAB_DIR)

EXTRA_DIST = $(mex_cfiles) $(m_SCRIPTS)

clean-local:
	rm -f *.$(MEXEXT)

.PHONY: clean
