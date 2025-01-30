import numpy as np
import matplotlib.pyplot as plt


def shift_hmax_to_zero(x, h):
    # find the index of the largest value of h
    i = np.argmax(h)
    xmax = x[i]
    x = ((x - xmax) + 30) % 30
    return x, h


def get_oomphlib_data(output_dir, dt, step):
    # timeseries of dh
    t = []
    dh = []
    k = 0
    for i in range(0, 100_000, step):
        filename = f"{output_dir}/spine_interface_{i}.dat"
        try:
            with open(filename) as fp:
                # dummy data (to be overwritten)
                dh.append(0.0)
                t.append(i * dt)

                # read the info from the file
                for line in fp:
                    if len(line) > 20:
                        # get h from the node
                        numbers = [float(s) for s in line.rstrip().split(' ')]
                        dh[k] = max(dh[k], numbers[1] - 1.0)
        except FileNotFoundError:
            final_i = i - step
            break
        k = k + 1
    t = np.array(t)
    dh = np.array(dh)


    # final interface
    x = []
    h = []
    filename = f"{output_dir}/spine_interface_{final_i}.dat"
    with open(filename) as fp:
        # read the info from the file
        for line in fp:
            if len(line) > 20:
                # get h from the node
                row = [float(s) for s in line.rstrip().split(' ')]
                x.append(row[0])
                h.append(row[1])
    x = np.array(x)
    h = np.array(h)
    x, h = shift_hmax_to_zero(x, h)

    return t, dh, x, h



def main():
    ########### GET OOMPH-LIB DATA ###########
    t, dh, x, h = get_oomphlib_data("output", 0.1, 10)
    # t250, dh250, x250, h250 = get_oomphlib_data("output-250", 0.1, 20)


    ########### GET BASILISK DATA ###########
    basilisk_out_dir = "../2d-film-no-control/out"

    # timeseries of dh
    t_basilisk = []
    dh_basilisk = []
    with open(basilisk_out_dir + "/ns-0.dat") as fp:
        for line in fp:
            # the row contains t, dh (2-norm), dh (inf-norm)
            row = [float(s) for s in line.rstrip().split(' ')]
            t_basilisk.append(row[0])
            dh_basilisk.append(row[2])
    t_basilisk = np.array(t_basilisk)
    dh_basilisk = np.array(dh_basilisk)

    # final interface
    x_basilisk = []
    h_basilisk = []
    with open(basilisk_out_dir + f"/ns-1-{470:010d}.dat") as fp:
        for line in fp:
            # skip first line with timestamp
            if line[0] == '#':
                continue

            # the row contains x, h, q
            row = [float(s) for s in line.rstrip().split(' ')]
            x_basilisk.append(row[0])
            h_basilisk.append(row[1])
    x_basilisk = np.array(x_basilisk)
    h_basilisk = np.array(h_basilisk)
    x_basilisk, h_basilisk = shift_hmax_to_zero(x_basilisk, h_basilisk)


    fig, ax = plt.subplots(1, 2)

    # plot growth
    ax[0].semilogy(t_basilisk, dh_basilisk, label='basilisk')
    ax[0].semilogy(t, dh, label='oomph-lib')
    # ax[0].semilogy(t250, dh250, linestyle='--', label='oomph-lib (n=250)')
    ax[0].legend()
    ax[0].set_xlabel('t')
    ax[0].set_ylabel('max(h-1)')

    # plot interface
    ax[1].scatter(x_basilisk, h_basilisk, label='basilisk')
    ax[1].scatter(x, h, label='oomph-lib')
    # ax[1].scatter(x250, h250, label='oomph-lib (n=250)')
    ax[1].legend()
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('h')

    plt.show()

if __name__ == '__main__':
    main()
