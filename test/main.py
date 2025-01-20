import numpy as np
from io import StringIO
from scipy.io import mminfo, mmread
from scipy.sparse import coo_matrix
import os
import sys

names = '''1138_bus.mtx
3x3.mtx
4x4.mtx
af23560.mtx
cage4.mtx
mbeacxc.mtx
sym3.mtx
sym4.mtx'''

def mm_name(mm_path):
    return os.path.basename(mm_path)

def test(mm_path, output_dir):
    with open(mm_path) as fp:
        try:
            text = fp.read()
            (_, _, _, _, field, _) = mminfo(StringIO(text))
            if field == 'complex':
                print(f'({mm_name(mm_path)}) complex matrix not supported')
                return
            mm = mmread(StringIO(text))
            v = np.random.rand(mm.shape[1], 1)
            result = mm.dot(v)
            vpath = os.path.join(output_dir, f'{mm_name(mm_path)}.vector')
            np.savetxt(vpath, v, delimiter='\n', fmt='%f')
            rpath = os.path.join(output_dir, f'{mm_name(mm_path)}.result')
            np.savetxt(rpath, result, delimiter='\n', fmt='%f')
        except ValueError as e:
            print(f'(matrix={mm_name(mm_path)}) {e}')

if __name__ == '__main__':
    matrices_dir = sys.argv[1]
    output_dir = sys.argv[2]

    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for mm_path in list(filter(lambda x: x.endswith('.mtx'), os.listdir(matrices_dir))):
            test(os.path.join(matrices_dir, mm_path), output_dir)
    except OSError as e:
        print(f'cannot create directory {output_dir}, got error {e}')
