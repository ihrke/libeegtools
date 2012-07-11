/** wrapping imageproc.h */

%typemap(in) (const Array *in, Array *dt) {
    if(!PyArray_Check($input)){
        PyErr_SetString(PyExc_ValueError, "Expecting a Numpy-Array");
        return NULL;
    }
    $1 = pyarray_to_myarray($input);
    if($1==NULL){
        PyErr_SetString(PyExc_ValueError, "conversion from np to myarray failed");
        return NULL;
    }
    $2 = NULL;
}
%typemap(out) Array* {
    if( !$1 ){
        PyErr_SetString(PyExc_ValueError, "Something wrong in function");
        return NULL;
    }
    $result = myarray_to_pyarray($1);
}
Array*  disttransform_deadreckoning( const Array *in, Array *dt );
%clear Array*;
