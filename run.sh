mkdir -p build && cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
./src/numsim "../$1"
cd ..
