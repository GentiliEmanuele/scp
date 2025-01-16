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
        mm = mmread(StringIO(fp.read()))
        v = np.random.rand(mm.shape[0])
        path = os.path.join(OUTPUT_DIR, f'{mm_name(mm_path)}.vector')
        np.savetxt(path, v, delimiter='\n', fmt='%f')
        result = mm.dot(v)
        path = os.path.join(OUTPUT_DIR, f'{mm_name(mm_path)}.result')
        np.savetxt(path, result, delimiter='\n', fmt='%f')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("see usage: program matrices_directory output_directory")
    MATRICES = sys.argv[1]
    OUTPUT_DIR = sys.argv[2]
    for mm_path in list(filter(lambda x: x.endswith('.mtx'), os.listdir(MATRICES))):
        test(os.path.join(MATRICES, mm_path))