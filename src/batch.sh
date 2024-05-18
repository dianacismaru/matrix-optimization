#!/bin/bash

rm out*
make clean
make

echo > output.txt

echo "Running BLAS..."
./tema3_blas ../input/input >> output.txt
echo "Comparing results..."
./compare out1 ../references/out1 0.001
./compare out2 ../references/out2 0.001
./compare out3 ../references/out3 0.001
echo "Running BLAS memcheck"
valgrind --tool=memcheck --leak-check=full ./tema3_blas ../input/input_valgrind > ../memory/blas.memory 2>&1
echo "Running BLAS cachegrind"
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema3_blas ../input/input_valgrind > ../cache/blas.cache 2>&1

rm out*
echo "Running NEOPT..."
./tema3_neopt ../input/input >> output.txt
echo "Comparing results..."
./compare out1 ../references/out1 0.001
./compare out2 ../references/out2 0.001
./compare out3 ../references/out3 0.001
echo "Running NEOPT memcheck"
valgrind --tool=memcheck --leak-check=full ./tema3_neopt ../input/input_valgrind > ../memory/neopt.memory 2>&1
echo "Running NEOPT cachegrind"
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema3_neopt ../input/input_valgrind > ../cache/neopt.cache 2>&1

rm out*
echo "Running OPT..."
./tema3_opt_m ../input/input >> output.txt
echo "Comparing results..."
./compare out1 ../references/out1 0.001
./compare out2 ../references/out2 0.001
./compare out3 ../references/out3 0.001
echo "Running OPT memcheck"
valgrind --tool=memcheck --leak-check=full ./tema3_opt_m ../input/input_valgrind > ../memory/opt_m.memory 2>&1
echo "Running OPT cachegrind"
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema3_opt_m ../input/input_valgrind > ../cache/opt_m.cache 2>&1

