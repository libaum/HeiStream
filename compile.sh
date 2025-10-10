#!/bin/bash

NCORES=12
ENABLE_TIME_MEASUREMENTS=${1:-OFF}  # Default: OFF, kann über Parameter überschrieben werden
unamestr=`uname`

# rm -rf build
if [ ! -d "build" ]; then
  mkdir build
fi

cd build
cmake ../ -DENABLE_TIME_MEASUREMENTS=$ENABLE_TIME_MEASUREMENTS
make -j $NCORES
cd ..

if [ ! -d "deploy" ]; then
  mkdir deploy
fi

cp ./build/buffcut deploy/
cp ./build/par_buffcut deploy/
# cp ./build/heistream_edge deploy/
# rm -r build

echo "Compiled with DENABLE_TIME_MEASUREMENTS=$ENABLE_TIME_MEASUREMENTS"
