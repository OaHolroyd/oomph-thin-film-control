import numpy as np
import matplotlib.pyplot as plt


def main():
    Nx = 5
    Ny = 4
    Np = 3
    Xperiodic = False
    Yperiodic = True

    Xmax = Nx * (Np-1)
    Ymax = Ny * (Np-1)

    spines = []
    count = 0

    # loop over the rows of elements in the grid
    for i in range(Ny):
        # loop over the elements in the row
        for j in range(Nx):
            e = i * Nx + j

            # loop over the rows of nodes in the element
            for ly in range(Np):
                # loop over the nodes in the row
                for lx in range(Np):
                    # construct the node
                    node = [(Np-1) * j + lx, (Np-1) * i + ly, count, e]

                    # account for periodicity if requried
                    if Xperiodic:
                        node[0] = node[0] % Xmax
                    if Yperiodic:
                        node[1] = node[1] % Ymax

                    # add to the list if it's not there already
                    to_add = True
                    for spine in spines:
                        if (spine[0] == node[0]) and (spine[1] == node[1]):
                            to_add = False
                            break

                    if to_add:
                        spines.append(node)
                        count = count + 1

    spines = np.array(spines)

    fig, ax = plt.subplots()

    # plot grid lines
    for i in range(Ny+1):
        ax.plot([0, Xmax], [i*(Np-1), i*(Np-1)], 'k')
    for j in range(Nx+1):
        ax.plot([j*(Np-1), j*(Np-1)], [0, Ymax], 'k')

    # plot points
    ax.scatter(spines[:, 0], spines[:, 1], c=spines[:, 3], cmap='prism')

    # label points
    for k in range(count):
        ax.annotate(f"{spines[k, 2]}", spines[k, 0:2])

    plt.show()
    plt.close()

if __name__ == '__main__':
  main()
