mkdir -p build && cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
time ./src/numsim "../$1"
cd ..
