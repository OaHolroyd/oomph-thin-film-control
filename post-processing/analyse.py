import numpy as np
import matplotlib.pyplot as plt


def shift_hmax_to_zero(x, h):
    # find the index of the largest value of h
    i = np.argmax(h)
    xmax = x[i]
    x = ((x - xmax) + 30) % 30
    return x, h


def get_oomphlib_data(output_dir):
    # timeseries of dh
    t = []
    dh = []
    k = 0
    for i in range(0, 100_000):
        filename = f"{output_dir}/1d_{i}.dat"
        try:
            with open(filename) as fp:
                # dummy data (to be overwritten)
                dh.append(0.0)
                t.append(0.0)

                # read the info from the file
                for line in fp:
                    if line[0] == '#':
                        t[k] = float(line.rstrip().split(' ')[2])
                    else:
                        # get h from the node
                        numbers = [float(s) for s in line.rstrip().split(' ')]
                        dh[k] = max(dh[k], numbers[1] - 1.0)
        except FileNotFoundError:
            break
        k = k + 1
    t = np.array(t)
    dh = np.array(dh)

    return t, dh


def main():
    fig, ax = plt.subplots()

    for nx, nz in [(100, 6)]:
        t, dh = get_oomphlib_data(f"output-{nx}-{nz}")
        ax.semilogy(t - 200, dh, label=f'{nx} x {nz}')

    ax.legend()
    ax.set_xlabel('t')
    ax.set_ylabel('max(h-1)')
    plt.show()


if __name__ == '__main__':
    main()
