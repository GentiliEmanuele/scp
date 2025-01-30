if [ "$#" -ne 4 ]; then
    echo "usage: $0 <file> (csr|hll) <output> <max number of threads>"
    exit 1
fi
if ["$2" -ne "csr" ] && ["$2" -ne "hll"]; then
    echo "usage: $0 <file> (csr|hll)"
    exit 1
fi
for i in $(seq 1 $4);
do
    echo "iteration number: $i"
    ./build/omp_time_"$2".exe "$1" $i 10 "$3"
done
