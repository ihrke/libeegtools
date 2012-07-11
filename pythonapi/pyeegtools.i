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
#include "imageproc.h"
#include "helper.h"
#include "array.h"             // my own array
#include "numpy/arrayobject.h" // numpy arrays

int dtype_npy_map[DT_END]={ NPY_INT16, NPY_UINT, NPY_INT, NPY_LONG, NPY_ULONG, NPY_FLOAT, NPY_DOUBLE };
DType npy_to_dtype( int nptype ){
    int r;
    switch(nptype){
     case NPY_DOUBLE:
        r=DOUBLE; break;
     case NPY_INT:
        r=INT; break;
     case NPY_LONG:
        r=LONG; break;
     default:
        warnprintf("no valid nptype encountered (=%i), using DOUBLE\n", nptype);
        warnprintf("NPY_INT=%i, NPY_LONG=%i\n", NPY_INT, NPY_LONG);
        r=DOUBLE; break;
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
        uint *dims=(uint*)malloc(nd*sizeof(uint));
        npy_intp *npdims=PyArray_DIMS(a);
        int i;
        for(i=0; i<nd; i++ )
            dims[i]=(uint)npdims[i];
        Array *out=array_fromptr( dt, nd, PyArray_DATA(a), dims );

        if( PyArray_ITEMSIZE(a)!=out->dtype_size || PyArray_NBYTES(a)!=out->nbytes ){
            errprintf("nd=%i, dt=%i, np_dt=%c\n", nd, dt, PyArray_DESCR(a)->kind);
            for( i=0; i<nd; i++ )
               errprintf("dims[%i]=%i\n", i, dims[i]);
            errprintf("np_itemsize=%i, ar_dtypesize=%i, np_nbytes=%li, arr_nbytes=%li\n",
                PyArray_ITEMSIZE(a), out->dtype_size, PyArray_NBYTES(a), out->nbytes);
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


%feature("autodoc", "1");

/* the modules */
%include imageproc.i
%include pythonfunctions.i

//%include eeg.i
//%include denoising.i
//%include distances.i
//%include som.i

