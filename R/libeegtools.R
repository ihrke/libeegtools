#
# These are some convenience functions to enable
# communication between R and libeegtools.
#

DT_CHAR=0;
DT_UINT=1;
DT_INT=2;
DT_LONG=3;
DT_ULONG=4;
DT_FLOAT=5;
DT_DOUBLE=6;

array.sizeof <- function( dtype ){
   if( dtype==DT_CHAR ){
	  return(1);
   } else if( dtype==DT_UINT ){
	  return(4);
   } else if( dtype==DT_INT ){
	 return(4);
   } else if( dtype==DT_LONG ){
	 return(8);
   } else if( dtype==DT_ULONG ){
	 return(8);
   } else if( dtype==DT_FLOAT ){
	 return(4);
   } else if( dtype==DT_DOUBLE ){
	 return(8);
   }
}

read.array.from.zz <- function( zz ){
   dtype <- readBin( zz, integer(), 1, 2 );
   if(length(dtype)==0 ) return( NULL );

   sizeof_dt <- readBin( zz, integer(), 1, 2 );
   ndim <- readBin( zz, integer(), 1, 2 );
   size <- rep(0,ndim);
   for( i in 1:ndim ){
	  size[i] <- readBin( zz, integer(), 1, 2 );
   }
   nbytes <- readBin( zz, integer(), 1, 4 );
   nel <- prod(size);
   data <- readBin( zz, double(), nel, sizeof_dt );

   array <- list(
	  size=size,
	  dtype=dtype,
	  sizeof_dt=sizeof_dt,
	  data=data
	  );
   class(array) <- list("list", "let_array");
   return( array );
}

read.array <- function( fname ){
   zz <- file( fname, "rb" );
   arrays = list();
   i=1;
   repeat{
	  a= read.array.from.zz( zz );
	  if( is.null(a) ) break;
	  arrays[[i]] = a;
	  i = i+1;
   }
   close( zz );
   return( arrays );
}

write.array.to.zz <- function( array, zz ){
   writeBin( as.integer(array$dtype), zz, size=2 );
   writeBin( as.integer(array$sizeof_dt), zz, size=2 );
   writeBin( as.integer(length(array$size)), zz, size=2 );
   writeBin( as.integer(array$size), zz, size=2 );
   writeBin( as.integer( prod(array$size)*array$sizeof_dt), zz, size=4 );
   writeBin( array$data, zz, size=array$sizeof_dt );
}

write.array <- function( array, fname, append=F ){
   if( append==T ){
	  zz <- file( fname, "ab");
   } else {
	  zz <- file( fname, "wb");
   }
   write.array.to.zz( array, zz );
   close( zz );
}

write.array.matrix <- function( mat, fname, append=F ){
   x = list();
   x$dtype=DT_DOUBLE;
   x$sizeof_dt=array.sizeof(DT_DOUBLE);
   x$size=rev(dim(mat)); # column-major stuff
   x$data=as.double(mat);
   write.array( x, fname, append=append );
}

as.matrix.let_array <- function( array ){
   if( length( array$size)!=2 ){
	  warning("The LET-array does not have 2 dimensions; using the first two dims");
   }
   x=matrix( nrow=array$size[2], ncol=array$size[1], data=array$data, byrow=F );
   return(x);
}
