#!/bin/bash
set -euo pipefail

src_dir=$1

#linux_version=$(lsb_release -si | tr '[:upper:]' '[:lower:]')-$(lsb_release -sr)

cget install -DCMAKE_C_FLAGS="-fPIC" -DCMAKE_CXX_FLAGS="-fPIC" -f ${src_dir}/requirements.txt --prefix /cget


# sav cli statically linked
mkdir -p ${src_dir}/linux-build
cd ${src_dir}/linux-build && rm -rf ./*
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=/cget/cget/cget.cmake \
  -DCMAKE_CXX_FLAGS="-static -static-libgcc -static-libstdc++" \
  -DCMAKE_BUILD_WITH_INSTALL_RPATH=1 \
  -DCPACK_GENERATOR="STGZ;DEB;RPM" \
  -DCPACK_PACKAGE_CONTACT="csg-devel@umich.edu" \
  ${src_dir}

make all manuals
make package
#cp minimac4-*.{sh,deb,rpm} ${out_dir}/
cd ..
