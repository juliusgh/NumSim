mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
./src/runTests
cd ..
