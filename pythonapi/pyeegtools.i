%module pyeegtools
%{
#define SWIG_FILE_WITH_INIT

%}

%include "typemaps.i"
%include "carrays.i"
%include "numpy.i"


%header %{
#define true 1
#define false 0
double* python_list_to_doubleptr( PyObject *p ){
  if (PyList_Check(p)) {
	 double *out;
	 int l = PyList_Size( p );
	 int i;
	 out = (double*)malloc( l*sizeof( double ) );
	 for( i=0; i<l; i++ ){
		out[i] = PyFloat_AsDouble( PyList_GetItem(p,i) );
	 }
	 return out;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

 char* python_string_to_charptr( PyObject *O ){
	if (PyString_Check(O)) {
	  return PyString_AsString(O);
	} else {
    PyErr_SetString(PyExc_TypeError,"not a string");
    return NULL;
	}
 } 
 PyObject* doubleptr_to_pythonlist( double *p, int n ){
	if( !p ){
	  PyErr_SetString(PyExc_TypeError,"cannot convert");
	  return NULL;
	}
	PyObject *l = PyList_New( n );
	int i;
	for( i=0; i<n; i++ ){
		 PyList_SetItem( l, i, PyFloat_FromDouble( p[i] ));
	}
   return l;
 }

 PyArrayObject* doubleptrptr_to_pyarray( double **X, int d1, int d2 ){
	double *d = (double*)malloc( d1*d2*sizeof(double));
	int i ;
	for( i=0; i<d1; i++ ){
	  memcpy( d+i*d2, X[i], d2*sizeof(double));
	}
	npy_intp dims[2];
	dims[0]=d1; dims[1]=d2;
	PyArrayObject* a=PyArray_SimpleNewFromData( 2, dims, NPY_DOUBLE, (void*)d);
	return a;
 }

 PyArrayObject* doubleptrptrptr_to_pyarray( double ***X, int d1, int d2, int d3 ){
	double *d = (double*)malloc( d1*d2*d3*sizeof(double));
	int i,j ;
	for( i=0; i<d1; i++ ){
	  for( j=0; j<d2; j++ ){
		 memcpy( d+(i*d2*d3)+(j*d3), X[i][j], d3*sizeof(double));
	  }
	}
	npy_intp dims[3];
	dims[0]=d1; dims[1]=d2; dims[2]=d3;
	PyArrayObject* a=PyArray_SimpleNewFromData( 3, dims, NPY_DOUBLE, (void*)d);
	return a;
 }

 
 PyArrayObject* myarray_to_pyarray( Array *a ){
	int i ;
	npy_intp dims[2];
	dims[0]=d1; dims[1]=d2;
	PyArrayObject* a=PyArray_SimpleNewFromData( 2, dims, NPY_DOUBLE, (void*)d);
	return a;
 }


%}

 /* converting python lists to double arrays */
%typemap(in) double* {
  $1 = python_list_to_doubleptr( $input );
}
%typemap(in) char* {
  $1 = python_string_to_charptr( $input );
}
%typemap(freearg) double* {
   free($1);
}
%typemap(freearg) char* {
   free($1);
}


%typemap(in) Array* {
  $1 = myarray_to_pyarray( $input );
}

/* low-level access to int* and double* */
%array_class(int, intArray);
%array_class(double, doubleArray);
%array_class(doubleArray, doublePtrArray);

%pythoncode %{
def buildMatrix( d1, d2 ):
  r = doublePtrArray(d1);
  for i in range(0,d1):
     r[i]=doubleArray(d2);
  return r;

%}

%init %{ // this is executed when the module is imported
  	 /* need to call this to enable numpy */
	 import_array(); 

%}

/* the modules */
//%include eeg.i
//%include denoising.i
//%include distances.i

 //%include pythonfunctions.i
