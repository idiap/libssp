#!/bin/sh
#
# Copyright 2013 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

# CMake finders
gitraw=https://raw.githubusercontent.com
if [ ! -e FindLibUBE.cmake ]
then
    wget $gitraw/pgarner/libube/master/cmake/FindLibUBE.cmake
fi

# Download a test file
if [ ! -e arctic_a0001.wav ]
then
    arctic=http://www.speech.cs.cmu.edu/cmu_arctic
    wget $arctic/cmu_us_bdl_arctic/wav/arctic_a0001.wav
fi
