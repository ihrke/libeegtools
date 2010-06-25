#!/bin/sh

# this script returns the correct revision for each file
# it is called by doxygen to put the tag into the documentation


EXTENSIONS="h c doc"

CURDIR=`pwd`
#echo $CURDIR
#echo $1
F=${1/$CURDIR/"."}
EXT=${F##*.}

FF=`basename $F`
#if [[ "$EXTENSIONS" =~ "${EXT}"  ]];
#then
cat filestatus | grep  -1 "$FF" |sed -n 's/^[ \]*Working revision:[ \t]*\([0-9][0-9\.]*\).*/\1/p'
#	cvs status $F | sed -n 's/^[ \]*Working revision:[ \t]*\([0-9][0-9\.]*\).*/\1/p'
#fi
