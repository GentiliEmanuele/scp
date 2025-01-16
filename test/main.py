import numpy as np
from io import StringIO
from scipy.io import mmread
from scipy.sparse import coo_matrix
import os
import sys

def mm_name(mm_path):
    return os.path.basename(mm_path)

def test(mm_path):
    with open(mm_path) as fp:
        try:
            mm = mmread(StringIO(fp.read()))
            v = np.random.rand(mm.shape[1], 1)
            result = mm.dot(v)
            vpath = os.path.join(OUTPUT_DIR, f'{mm_name(mm_path)}.vector')
            np.savetxt(vpath, v, delimiter='\n', fmt='%f')
            rpath = os.path.join(OUTPUT_DIR, f'{mm_name(mm_path)}.result')
            np.savetxt(rpath, result, delimiter='\n', fmt='%f')
        except ValueError as e:
            print(f'(matrix={mm_name(mm_path)}) {e}')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("see usage: program matrices_directory output_directory")
    MATRICES = sys.argv[1]
    OUTPUT_DIR = sys.argv[2]
    for mm_path in list(filter(lambda x: x.endswith('.mtx'), os.listdir(MATRICES))):
        test(os.path.join(MATRICES, mm_path))