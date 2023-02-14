#!/bin/bash
# Usage: docker build -t packaging-alpine3.17 -f packaging-dockerfile-alpine3.17 . 
#        docker run -v `pwd`:/app -w /app packaging-alpine3.17 ./package-linux.sh
set -euo pipefail

src_dir=`pwd`
build_dir=${src_dir}/pkg-build-alpine

rm -rf ${build_dir}
mkdir ${build_dir}
cd ${build_dir}

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"

cmake -P ../dependencies.cmake deps/ -DSHRINKWRAP_PREFER_STATIC=ON
rm deps/lib/libz.so
unset CFLAGS
unset CXXFLAGS

arc=`uname -m`

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=$(pwd)/deps \
  -DCMAKE_CXX_FLAGS="-I$(pwd)/deps/include" \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DCMAKE_EXE_LINKER_FLAGS="-static" \
  -DCPACK_GENERATOR="STGZ" \
  -DCPACK_SYSTEM_NAME="Linux-${arc}" \
  -DCPACK_RESOURCE_FILE_LICENSE=${src_dir}/LICENSE \
  -DCPACK_PACKAGE_CONTACT="csg-devel@umich.edu" \
  ${src_dir}

#-DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" 

make
make manuals
make package

mkdir -p ${src_dir}/pkg/
cp minimac4-*.sh ${src_dir}/pkg/
cd ${src_dir}
rm -r ${build_dir}
