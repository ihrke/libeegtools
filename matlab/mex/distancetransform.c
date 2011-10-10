/**\file distancetransform.c
	\brief Matlab wrapper for \ref disttransform_deadreckoning().

	 This distance-transform (DT) of a binary image includes at each
	 pixel the shortest distance to the closest nonzero point in
	 the input image.

	 This implementation uses the dead-reckoning algorithm from
	 \verbatim
	 @article{grevera2004dead,
  	   title={{The" Dead reckoning" signed distance transform}},
	   author={Grevera, G.J.},
	   journal={Computer Vision and Image Understanding},
	   volume={95},
	   number={3},
	   pages={317--333},
	   year={2004}
	 }
	 \endverbatim
	The MATLAB-interface is documented in the corresponding M-file.
*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "distances.h"
#include "imageproc.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
  char msg[2000];
  const char *docstring = 
	 "[ dt ] = distancetransform( binarymatrix );\n"
	 " The input array is treated as binary, i.e. if a value is >0, it is\n"
	 " assumed to be 1, else it is 0.\n";
  
  /* check proper input and output */
  if(nrhs<1){
    sprintf(msg, ">> Need 1 input.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!is_mex_matrix(prhs[0])){
    sprintf(msg, ">> First Input must be a double-precision matrix.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } 

  /* getting the appropriate array representation */
  Array *binimg = mex_mxarray_to_array( prhs[0] );

  /* computation */
  array_typecast( binimg, INT );
  Array *dt = disttransform_deadreckoning( binimg, NULL );
  plhs[0]=mex_array_to_mxarray( dt ); 

  /* cleaning up */
  array_free( binimg );
  array_free( dt );

  return;
}
        
