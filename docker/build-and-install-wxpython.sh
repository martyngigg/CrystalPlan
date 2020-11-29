#!/bin/bash
# Build and install wxpython based on the
# CentOS 7 version of wxGTK and conda python.
# It is assumed that wxGTK==2.8 is already
# installed along with devtoolset-7, scl & python.
set -e
SRC=https://sourceforge.net/projects/wxpython/files/wxPython/2.8.12.1/wxPython-src-2.8.12.1.tar.bz2/download
DEST=wxPython-src-2.8.12.1.tar.bz2
WXDIR=$(basename -s .tar.bz2 $DEST)
PYTHON=/opt/conda/bin/python

function cleanup() {
  rm -fr /tmp/$WXDIR
  rm -fr /tmp/$DEST
}
trap cleanup EXIT

# download and extract
cd /tmp
curl --progress-bar --location $SRC -o $DEST
tar --extract --bzip2 -f $TARGET
cd $WXDIR/wxPython
scl enable devtoolset-7 "$PYTHON setup.py install"

