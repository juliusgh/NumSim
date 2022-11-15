if command -v module &> /dev/null
then
    source ./load.sh
fi
mkdir -p build && cd build || exit
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
time ./src/numsim "../$1"
cd ..
