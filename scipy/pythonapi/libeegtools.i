%module libeegtools
%{
#define SWIG_FILE_WITH_INIT
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

%}

%include "typemaps.i"
%include "carrays.i"


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


/* low-level access to int* and double* */
%array_class(int, intArray);
%array_class(double, doubleArray);

/* the modules */
%include mathadd.i
%include eeg.i
%include eegio.i
%include denoising.i
 /*%include som.i*/


%include pythonfunctions.i
