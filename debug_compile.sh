#!/bin/bash

NCORES=12
unamestr=`uname`

if [ ! -d "debug_build" ]; then
  mkdir debug_build
fi

cd debug_build
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-mno-avx -mno-avx2 -mno-fma"  ../ # -fsanitize=address,leak -fno-omit-frame-pointer"
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-mno-avx -mno-avx2 -mno-fma" ../
make -j $NCORES
cd ..

if [ ! -d "deploy" ]; then
  mkdir deploy
fi

cp ./debug_build/heistream deploy/debug_heistream
cp ./debug_build/par_heistream deploy/debug_par_heistream
