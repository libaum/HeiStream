#!/bin/bash

NCORES=12
ENABLE_TIME_MEASUREMENTS=${1:-OFF}  # Default: OFF, kann über Parameter überschrieben werden
unamestr=`uname`

rm -rf build
mkdir build
cd build 
cmake ../
cmake ../ -DENABLE_TIME_MEASUREMENTS=$ENABLE_TIME_MEASUREMENTS
make -j $NCORES
cd ..

mkdir deploy
cp ./build/heistream deploy/
cp ./build/heistream_edge deploy/
rm -r build