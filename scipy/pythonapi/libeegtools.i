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
		printf("  Recv: %f\n", out[i]);
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


%array_class(int, intArray);
%array_class(double, doubleArray);

%include mathadd.i
%include eeg.i
%include eegio.i
