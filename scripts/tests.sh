mkdir -p build && cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
./src/runTests
cd ..
