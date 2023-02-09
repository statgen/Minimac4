#!/bin/bash
# Usage: docker build -t packaging-ubuntu20 -f packaging-dockerfile-ubuntu20 . 
#        docker run -v `pwd`:/app -w /app packaging-ubuntu20 ./package-linux.sh
set -euo pipefail

src_dir=`pwd`
build_dir=${src_dir}/pkg-build

rm -rf ${build_dir}
mkdir ${build_dir}
cd ${build_dir}

export CFLAGS="-fPIC"
export CXXFLAGS="-fPIC"
#cget ignore xz
cget install zlib,http://zlib.net/fossils/zlib-1.2.11.tar.gz
rm cget/lib/libz.so
cget install -f ${src_dir}/requirements.txt
unset CFLAGS
unset CXXFLAGS


cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=cget/cget/cget.cmake \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=ON \
  -DCMAKE_EXE_LINKER_FLAGS="-static -static-libgcc -static-libstdc++" \
  -DCPACK_GENERATOR="STGZ;DEB;RPM" \
  -DCPACK_PACKAGE_CONTACT="csg-devel@umich.edu" \
  ${src_dir}

make
#make manuals
make package

mkdir -p ${src_dir}/pkg/
cp minimac4-*.{sh,deb,rpm} ${src_dir}/pkg/
cd ${src_dir}
rm -r ${build_dir}
