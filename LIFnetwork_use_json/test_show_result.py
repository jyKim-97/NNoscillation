import numpy as np
import matplotlib.pyplot as plt


def load_result(fname):
    ts = []
    vs = []
    with open(fname, "r") as fid:
        line = fid.readline()
        while line:
            vals = line.split(",")
            ts.append(float(vals[0]))
            vs.append(np.array([float(v) for v in vals[1:-1]]))
            line = fid.readline()
    ts = np.array(ts)
    vs = np.array(vs)

    return ts, vs


if __name__ == "__main__":

    """ load data """
    ts, vs = load_result("test_out.csv")
    ncells = len(vs[0])
    print(vs.shape)

    """ Show result """
    plt.figure(dpi=150, figsize=(6, 6))

    plt.subplot(211)
    for i in range(ncells):
        plt.plot(ts, vs[:, i], lw=0.5)
    plt.xlabel("times (ms)", fontsize=12)
    plt.ylabel("voltage (mV)", fontsize=12)
    
    plt.subplot(212)
    v_sum = np.average(vs, axis=1)
    plt.plot(ts, v_sum, "k", lw=1)
    plt.xlabel("times (ms)", fontsize=12)
    plt.ylabel("voltage (mV)", fontsize=12)

    plt.tight_layout()
    plt.show()
    