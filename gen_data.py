#!/usr/bin/python
import numpy as np
import sys

def gen_data(argv):
    n = int(argv[1])
    p = 100*np.random.ranf((n, 3))-50
    v = 10*np.random.ranf((n, 3))-5
    m = 100*np.random.ranf((n, 1))
    pvm = np.hstack((p, v, m))
    np.savetxt(argv[2], pvm, delimiter=',')

if __name__ == "__main__":
    gen_data(sys.argv)
