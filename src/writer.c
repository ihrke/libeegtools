#include "writer.h"

void write_double_matrix_ascii_file(const char *fname, const double **d, int xdim, int ydim){
  FILE *f;
  int x, y;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_double_matrix_ascii(f, d, xdim, ydim);

  fclose(f);
}

void write_double_matrix_ascii(FILE *out, const double **d, int xdim, int ydim){
  int x, y;

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fprintf(out, "%f ", d[y][x]);
	 }
	 fprintf(out, "\n");
  }
}
