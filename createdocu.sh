#!/bin/sh

DOXYGEN="$HOME/install/doxygen-svn/bin/doxygen"

echo ">>> Removing current documentation in ./doc/html <ENTER>"
read test
rm -rf doc/html

if ! test "${1}" = nocvs ; then
	echo ">>> Getting CVS status for correct version numbers"
	cvs status > filestatus
fi

remove_experimental=yes
if test "${remove_experimental}" = yes ; then

echo ">>> Exclude Experimental Files from Documentation"

experimental_files=""

cfiles=`ls src/*.[hc]`
for cfile in $cfiles ; do
	 z=`cat $cfile | grep experimental_do_not_document`
	 if test ${#z} -gt 0 ; then
		  experimental_files="$experimental_files $cfile"
	 fi
done

echo ">>> Excluding files: $experimental_files"

#experimental_files=${experimental_files//\//\\/}
#experimental_files=${experimental_files//\ /\\ /}

echo ">>> Excluding files: $experimental_files"
cat Doxyfile | sed s/EXCLUDE\ *=/EXCLUDE="$experimental_files"/ | $DOXYGEN - &> doxygen.log

else

echo ">>> Generating Documentation"
$DOXYGEN &> doxygen.log

fi

echo ">>> Deleting filestatus"
rm filestatus
