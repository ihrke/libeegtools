#include "writer.h"

void write_double_matrix_ascii(const char *fname, const double **d, int xdim, int ydim){
  FILE *f;
  int x, y;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fprintf(f, "%f ", d[y][x]);
	 }
	 fprintf(f, "\n");
  }
  fclose(f);
}
