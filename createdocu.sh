#!/bin/sh

DOXYGEN="$HOME/install/doxygen-svn/bin/doxygen"

echo ">>> Removing current documentation in ./doc/html <ENTER>"
read test
rm -rf doc/html

echo ">>> Getting CVS status for correct version numbers"
cvs status > filestatus

echo ">>> Generating Documentation"
$DOXYGEN

echo ">>> Deleting filestatus"
rm filestatus
