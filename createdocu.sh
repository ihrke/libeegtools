#!/bin/sh

DOXYGEN="$HOME/install/doxygen-svn/bin/doxygen"

echo ">>> Removing current documentation in ./doc/html <ENTER>"
read test
rm -rf doc/html

if ! test -n $1 ; then
	echo ">>> Getting CVS status for correct version numbers"
	cvs status > filestatus
fi

echo ">>> Generating Documentation"
$DOXYGEN &> doxygen.log

echo ">>> Deleting filestatus"
rm filestatus
