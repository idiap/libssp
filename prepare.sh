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
if [ ! -e cmake/FindLibUBE.cmake ]
then
    curl -LO $gitraw/pgarner/libube/master/cmake/FindLibUBE.cmake
    mkdir -p cmake
    mv FindLibUBE.cmake cmake
fi

# Download a test file
if [ ! -e test/test.wav ]
then
    arctic=http://www.speech.cs.cmu.edu/cmu_arctic
    curl -LO $arctic/cmu_us_bdl_arctic/wav/arctic_a0001.wav
    mv arctic_a0001.wav test/test.wav
fi
