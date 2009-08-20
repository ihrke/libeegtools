#!/usr/bin/env python

"""
setup.py file for SWIG wrapper
"""

from distutils.core import setup, Extension
import os

gslcflags  = os.popen( "gsl-config --cflags" ).read().strip();
gslldflags = os.popen( "gsl-config --libs" ).read().strip();
gsllibdirs = []
for l in gslldflags.split(" "):
    if l.strip().startswith("-L"):
        gsllibdirs.append(l.strip()[2:]);

libeegtools_module = Extension('_libeegtools',
                               sources=['libeegtools_wrap.c'],
                               include_dirs=['../../src'],
                               library_dirs=['../../src/.libs'],
                               extra_compile_args=[gslcflags],
                               extra_link_args=[gslldflags, "-leegtools"],
                               runtime_library_dirs=gsllibdirs,
                               )

setup (name = 'libeegtools',
       version = '0.1',
       author      = "Matthias Ihrke",
       description = """Wrapper for LibEEGTools""",
       ext_modules = [libeegtools_module],
       py_modules = ["libeegtools"],
       )