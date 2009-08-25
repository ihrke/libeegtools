#!/bin/sh



CURDIR=`pwd`
echo $CURDIR
cd ./src
find . -name "*.[chCH]" -print | etags -
cd $CURDIR
