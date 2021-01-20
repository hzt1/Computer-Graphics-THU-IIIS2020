#!/usr/bin/env bash

# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    mkdir -p build
    cd build
    cmake ..
    cd ..
fi

# Build project.
cd build
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output
# bin/PA1 testcases/cube_glass.txt output/cube_glass.bmp
# time bin/PA1 testcases/bunny.txt output/bunny.bmp
#time bin/PA1 testcases/cornell_js.txt output/cornell_js.bmp
#time bin/PA1 testcases/cornell_js2.txt output/focus.bmp
# time  bin/PA1 testcases/cornell.txt output/cornell.bmp
 #time  bin/PA1 testcases/cornell_tt.txt output/cornell_tt_s.bmp
# time bin/PA1 testcases/cornell_s.txt output/cornell_s.bmp
# bin/PA1 testcases/scene01_basic.txt output/scene01.bmp
# time bin/PA1 testcases/cornell_g.txt output/cornell_g.bmp
#  time bin/PA1 testcases/cornell_g_p.txt output/cornell_g_p.bmp
# time bin/PA1 testcases/glass.txt output/glass.bmp
time bin/PA1 testcases/arsenal.txt output/arsenal.bmp