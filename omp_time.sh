if [ "$#" != 4 ]; then
    echo "usage: $0 <file> (csr|hll) <max number of threads> <output dir>"
    exit 1
fi
if ["$2" != "csr" ] && ["$2" != "hll"]; then
    echo "usage: $0 <file> (csr|hll)"
    exit 1
fi
for i in 1 12 24
do
    export OMP_NUM_THREADS=$i
    echo "iteration number $i with #threads=$OMP_NUM_THREADS"
    ./build/omp_time_"$2".exe "$1" $i 50 "$4/$2_$i"
done
