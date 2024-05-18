#!/bin/bash

make clean
make

echo > output.txt
./tema3_blas ../input/my_input >> output.txt
./tema3_neopt ../input/my_input >> output.txt
./tema3_opt_m ../input/my_input >> output.txt


