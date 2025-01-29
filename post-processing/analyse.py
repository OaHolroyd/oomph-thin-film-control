import numpy as np
import matplotlib.pyplot as plt


def shift_hmax_to_zero(x, h):
    # find the index of the largest value of h
    i = np.argmax(h)
    xmax = x[i]
    x = ((x - xmax) + 20) % 20
    return x, h


def main():
    ########### GET OOMPH-LIB DATA ###########
    dt = 0.1
    step = 20
    t = []
    dh = []

    # timeseries of dh
    k = 0
    for i in range(0, 100_000, step):
        filename = f"output/spine_interface_{i}.dat"
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
    filename = f"output/spine_interface_{final_i}.dat"
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
    ax[0].legend()
    ax[0].set_xlabel('t')
    ax[0].set_ylabel('max(h-1)')

    # plot interface
    ax[1].scatter(x_basilisk, h_basilisk, label='basilisk')
    ax[1].scatter(x, h, label='oomph-lib')
    ax[1].legend()
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('h')

    plt.show()



    return

    k = 0;
    for i in range(0, 100_000, 10):
        filename = f"output/spine_interface_{i}.dat"

        hmax.append(0.0)
        t.append(i * dt)

        try:
            with open(filename) as fp:
                for line in fp:
                    if len(line) > 20:
                        # get h from the node
                        numbers = [float(s) for s in line.rstrip().split(' ')]
                        hmax[k] = max(hmax[k], abs(numbers[1] - 1.0))
        except FileNotFoundError:
            print(f"file {filename} not found")
            final_i = i
            hmax.pop()
            t.pop()
            break

        k = k + 1

    plt.semilogy(t, hmax, label='oomphlib')
    plt.semilogy(9999+hmax_basilisk[:, 0], hmax_basilisk[:, 1], label='basilisk')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('interfacial perturbation')
    plt.show()

    # also plot the travelling waves

if __name__ == '__main__':
    main()
