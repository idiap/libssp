#!/bin/sh
#
# Copyright 2013 by Idiap Research Institute, http://www.idiap.ch
#
# See the file COPYING for the licence associated with this software.
#
# Author(s):
#   Phil Garner, December 2013
#

rm -rf CMakeCache.txt CMakeFiles cmake_install.cmake

# Download a test file
if [ ! -e arctic_a0001.wav ]
then
    arctic=http://www.speech.cs.cmu.edu/cmu_arctic
    wget $arctic/cmu_us_bdl_arctic/wav/arctic_a0001.wav
fi

export CC=clang
export CXX=clang++

export CPATH=~/local/include

cmake \
    -D CMAKE_BUILD_TYPE=debug \
    -D CMAKE_INSTALL_PREFIX=~/local \
    .
