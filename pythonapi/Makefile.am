# set the include path found by configure
AM_CPPFLAGS = -Isrc $(all_includes) -fPIC --std=gnu99 -Wall @MATIO_CFLAGS@ @DEBUG_FLAGS@ `@GSL_CONFIG@ --cflags`



if BUILD_PYAPI

# BUILT_SOURCES contains source files that must be built prior
# to any compilation. Here it is used to call swig on the 
# interface files.
BUILT_SOURCES = pyeegtools_wrap.c
INTERFACES = pyeegtools.i i numpy.i pythonfunctions.i imageproc.i
EXTRA_DIST=$(INTERFACES)

pyeegtools_wrap.c: $(INTERFACES)
	echo "Swigging interface files"
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -I$(top_srcdir)/src -outdir ./ -o $@ $<

all-local:
	$(PYTHON) setup.py build_ext --inplace

clean-local:
	-rm -rf build
	-rm pyeegtools_wrap.c _pyeegtools.so

endif
