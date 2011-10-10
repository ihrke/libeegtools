[ dt ] = distancetransform( binarymatrix );
% DISTANCETRANSFORM
%	 This distance-transform (DT) of a binary image includes at each
%	 pixel the shortest distance to the closest nonzero point in
%	 the input image.
%
% [ dt ] = distancetransform( binarymatrix );
%  The input array is treated as binary, i.e. if a value is >0, it is
%  assumed to be 1, else it is 0. Input must be double.
%
%	 This implementation uses the dead-reckoning algorithm from
%	 @article{grevera2004dead,
% 	   title={{The" Dead reckoning" signed distance transform}},
%	   author={Grevera, G.J.},
%	   journal={Computer Vision and Image Understanding},
%	   volume={95},
%	   number={3},
%	   pages={317--333},
%	   year={2004}
%	 }
%
