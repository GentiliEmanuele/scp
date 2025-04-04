if [ "$#" -ne 3 ]; then
    echo "usage: $0 <name_file> <max number of threads> <output_dir>"
    exit 1
fi
for i in $(seq 1 $2);
do
    echo "iteration number: $i"
    ./build/omp_time_csr.exe "$1" $i 1000 "$3/csr_$i.csv"
done
