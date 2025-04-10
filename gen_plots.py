from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    nx = 40
    ny = 20
    n = nx * ny

    x = np.arange(nx) + 0.5
    y = np.arange(ny) + 0.5
    x, y = np.meshgrid(x, y)

    A = np.loadtxt("output/A_wr.dat")
    K = np.loadtxt("output/K.dat")
    B = np.loadtxt("output/B_wr.dat")
    F = np.loadtxt("output/F.dat")

    print(A.shape)
    print(K.shape)
    print(B.shape)
    print(F.shape)

    i = 10

    # plot a row of K
    k = K[i, :n]
    k = k.reshape(ny, nx)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, k)
    fig.savefig("output/Kh.png")
    plt.close(fig)

    k = K[i, n:]
    k = k.reshape(ny, nx)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, k)
    fig.savefig("output/Kq.png")
    plt.close(fig)

    # plot a column of B
    b = B[:n, i]
    b = b.reshape(ny, nx)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, b)
    fig.savefig("output/Bh.png")
    plt.close(fig)

    b = B[n:, i]
    b = b.reshape(ny, nx)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, b)
    fig.savefig("output/Bq.png")
    plt.close(fig)

    # plot a column of F
    f = F[:n, i]
    f = f.reshape(ny, nx)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.plot_surface(x, y, f)
    fig.savefig("output/F.png")
    plt.close(fig)

    # sparsity patterns of A
    fig, ax = plt.subplots()
    ax.spy(A)
    ax.plot([-0.5, 2*n-0.5], [n-0.5, n-0.5])
    ax.plot([n-0.5, n-0.5], [-0.5, 2*n-0.5])
    fig.savefig("output/A.png")
    plt.close(fig)

    # get the data at a given timestep
    tmax = []
    hmax = []
    i = 0
    while Path(f"output/surface_{i}.dat").exists():
        t = None
        file = f"output/surface_{i}.dat"
        h = np.zeros((n,))
        q = np.zeros((n,))
        f = np.zeros((n,))
        j = 0
        with open(file, "r") as fp:
            while line := fp.readline():
                line = line.strip()

                # empty line
                if len(line) < 1:
                    continue

                # header comment
                if line[0] == "#":
                    t = float(line.split(' ')[2])
                else:
                    data = [float(x) for x in line.split(' ')]
                    h[j] = data[2]
                    q[j] = data[3]
                    f[j] = data[5]
                    j += 1

        ff = -F @ (K @ np.concatenate([h - 1.0, q - 2.0/3.0]))

        h = h.reshape(ny, nx)
        q = q.reshape(ny, nx)
        f = f.reshape(ny, nx)
        ff = ff.reshape(ny, nx)

        tmax.append(t)
        hmax.append(np.max(np.abs(h-1)))

        if i == 2001:
            fig, ax = plt.subplots(1, 3, figsize=(12, 7), subplot_kw={"projection": "3d"})
            ax[0].plot_surface(x, y, h)
            ax[1].plot_surface(x, y, f)
            ax[2].plot_surface(x, y, ff)
            fig.suptitle(f"t = {t}")
            fig.savefig(f"output/plot-{i}.png")
            plt.close(fig)

        i += 1
        print(i)

    fig, ax = plt.subplots()
    ax.semilogy(tmax, hmax)
    fig.savefig("output/hmax.png")
    plt.close(fig)


    # try:
    #     # plot columns of L
    #     L = np.loadtxt(f'out/L.dat')
    #     n = len(L[:, 0]) // 2
    #     x = np.arange(2*n)
    #     fig, ax = plt.subplots()
    #     plt.plot(x, L)
    #     fig.savefig("plots/L.png")
    #     plt.close(fig)
    # except:
    #     pass

    # try:
    #     # plot rows of K
    #     K = np.loadtxt(f'out/K.dat')
    #     x = np.arange(max(*K.shape))
    #     fig, ax = plt.subplots()
    #     plt.plot(x, K.T)
    #     fig.savefig("plots/K.png")
    #     plt.close(fig)
    # except:
    #     pass

    # try:
    #     # plot columns of B
    #     B = np.loadtxt(f'out/Bcf.dat')
    #     x = np.arange(2*n)
    #     fig, ax = plt.subplots()
    #     plt.plot(x, B)
    #     fig.savefig("plots/B.png")
    #     plt.close(fig)
    # except:
    #     pass

    # try:
    #     # sparsity patterns of A
    #     A = np.loadtxt(f'out/Acf.dat')
    #     Awr = np.loadtxt(f'out/A_wr.dat')
    #     fig, ax = plt.subplots(1, 2)
    #     ax[0].spy(A)
    #     ax[0].plot([-0.5, 2*n-0.5], [n-0.5, n-0.5])
    #     ax[0].plot([n-0.5, n-0.5], [-0.5, 2*n-0.5])
    #     ax[1].spy(Awr)
    #     ax[1].plot([-0.5, 2*n-0.5], [n-0.5, n-0.5])
    #     ax[1].plot([n-0.5, n-0.5], [-0.5, 2*n-0.5])
    #     fig.savefig("plots/A.png")
    #     plt.close(fig)
    # except:
    #     pass


    # pde = 'ns'

    # # Get 1D data
    # data = np.loadtxt(f'out/{pde}-0.dat')
    # t = data[:, 0]
    # dh = data[:, 1]
    # de = data[:, 2]
    # dc = data[:, 3]
    # c = data[:, 4]


    # # Plot 1D data
    # fig, ax = plt.subplots()
    # ax.semilogy(t, dh)
    # ax.semilogy(t, de)
    # fig.savefig("plots/lines.png")
    # plt.close(fig)


    # # Plot 2D frames
    # # plot dummy data
    # fig, ax = plt.subplots()
    # data = np.loadtxt(f'out/{pde}-1-{0:010d}.dat')
    # x = data[:, 0]
    # dx = x[1] - x[0]

    # hplot, = plt.plot(x, x)
    # fplot, = plt.plot(x, x)
    # zplot, = plt.plot(x, x)

    # plt.axis([0, 30, 0.3, 1.7])


    # for i in range(len(t)):
    #     # Get 2D data for the ith step
    #     data = np.loadtxt(f'out/{pde}-1-{i:010d}.dat')
    #     x = data[:, 0]
    #     h = data[:, 1]
    #     f = data[:, 2]
    #     z = data[:, 3]
    #     # q = data[:, 4]

    #     # TODO: could be faster if we just change the ydata
    #     hplot.set_ydata(h)
    #     fplot.set_ydata(1+f)
    #     zplot.set_ydata(z)
    #     plt.title(f'time {t[i]} [step {i}]')

    #     fig.savefig(f"plots/{i}.png")

    # Turn plots into a gif
    # TODO


if __name__ == '__main__':
    main()
