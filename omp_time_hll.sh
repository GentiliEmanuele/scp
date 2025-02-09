if [ "$#" -ne 3 ]; then
    echo "usage: $0 <file> <max number of threads> <output_dir>"
    exit 1
fi
for i in 1 12 24;
do
    echo "iteration number: $i"
    ./build/omp_time_hll.exe "$1" 50 $i 32 "$3/hll_$i"
done
