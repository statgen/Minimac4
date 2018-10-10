#!/bin/sh

rm -rf cget/ release-build/ install.log
echo -e "Installing Dependencies - Libstatgen ..."
cget install -f requirements.txt &> install.log
mkdir release-build
cd release-build/
echo -e "Generating MakeFiles ..."
cmake -DCMAKE_CXX_FLAGS=-pg -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
make
echo "Binary created at /release-build/minimac4"

