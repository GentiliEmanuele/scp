if [ "$#" -ne 3 ]; then
    echo "usage: $0 <file> <max number of threads> <output_dir>"
    exit 1
fi
for i in $(seq 1 $2);
do
    echo "iteration number: $i"
    ./build/omp_time_hll.exe "$1" 32 $i 100 "$3/hll_$i"
done
