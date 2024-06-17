#!/bin/bash

NCORES=12
unamestr=`uname`

if [ ! -d "build" ]; then
  mkdir build
fi

cd build
cmake ../
make -j $NCORES
cd ..

if [ ! -d "deploy" ]; then
  mkdir deploy
fi

cp ./build/heistream deploy/