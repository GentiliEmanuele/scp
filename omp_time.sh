if [ $# -ne 4 ] && [ $# -ne 5 ]; then
    echo "usage: $0 file (csr|hll) {hack_size} max_num_threads output_dir"
    exit 1
fi
if [ "$2" != "csr" ] && [ "$2" != "hll" ]; then
    echo "expected either csr or hll but got $2"
    exit 1
fi
if [ "$2" == "csr" ]; then
    max_num_threads=$3
else
    max_num_threads=$4
fi
for i in $(seq 1 $max_num_threads);
do
    export OMP_NUM_THREADS=$i
    echo "iteration number $i with #threads=$OMP_NUM_THREADS"
    if [ "$2" == "csr" ]; then
        ./build/omp_time_"$2".exe "$1" $i 100 "$4/$2_$i"
    else
        ./build/omp_time_"$2".exe "$1" $3 $i 100 "$4/$2_$i"
    fi
done
