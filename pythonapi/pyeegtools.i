%define DOCSTRING
"
This is a python-wrapper for some functionality provided by the libeegtools library.
"
%enddef

%module(docstring=DOCSTRING) pyeegtools
%{
#define SWIG_FILE_WITH_INIT

%}

%include "typemaps.i"
%include "numpy.i"


%header %{
#include "array.h"     // my own array
#include "numpy/arrayobject.h" // numpy arrays
int dtype_npy_map[DT_END]={ NPY_INT16, NPY_UINT, NPY_INT, NPY_LONG, NPY_ULONG, NPY_FLOAT, NPY_DOUBLE };

int npy_to_dtype( int nptype ){
    int r;
    switch(nptype){
     case NPY_DOUBLE:
        r=DOUBLE; break;
     case NPY_INT:
        r=INT; break;
     default:
        r=DOUBLE;
    }
    return r;
}

 PyArrayObject* myarray_to_pyarray( Array *a ){
        npy_intp *dims=(npy_intp*)malloc( (a->ndim)*sizeof(npy_intp));
        int i;
        for(i=0; i<a->ndim; i++ ){
            dims[i]=a->size[i];
        }

        PyArrayObject* out=PyArray_SimpleNewFromData( a->ndim, dims, dtype_npy_map[a->dtype], (void*)(a->data));
        return out;
 }

 Array* pyarray_to_myarray( PyObject *obj ){
         if(!PyArray_Check(obj)){
            return NULL;
         }
         PyArrayObject *a=(PyArrayObject*)obj;

        int nd=PyArray_NDIM(a);
        DType dt=npy_to_dtype(PyArray_TYPE(a));
        Array *out=array_fromptr( dt, nd, PyArray_DATA(a), (const unsigned int *)PyArray_DIMS(a) );
        if( PyArray_ITEMSIZE(a)!=out->dtype_size || PyArray_NBYTES(a)!=out->nbytes ){
            array_free(out);
            return NULL;
         }

        return out;
 }
%}


%init %{ // this is executed when the module is imported
  	 /* need to call this to enable numpy */
	 import_array(); 
%}



%pythoncode %{
import numpy as np

%}

%feature("autodoc", "1");

/* the modules */
%include imageproc.i
%include pythonfunctions.i

//%include eeg.i
//%include denoising.i
//%include distances.i
//%include som.i

