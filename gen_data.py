#!/usr/bin/python
import numpy as np
import sys

def gen_data(argv):
    n = int(argv[1])
    maxlen = float(argv[2])
    p = maxlen*np.random.ranf((n, 3))-0.5*maxlen
    v = 0.01*maxlen*np.random.ranf((n, 3))-(0.005*maxlen)
    m = maxlen*np.random.ranf((n, 1))
    pvm = np.hstack((p, v, m))
    np.savetxt(argv[3], pvm, delimiter=',')

if __name__ == "__main__":
    gen_data(sys.argv)
