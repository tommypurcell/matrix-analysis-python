import numpy as np
import matplotlib.pyplot as plt

def plotstructure(coord, ends):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Convert ends to zero-based indexing in place
    ends -= 1

    # Plot nodes
    ax.scatter(coord[:, 0], coord[:, 1], coord[:, 2], color='red')

    # Plot elements
    for i in range(ends.shape[0]):
        start_node = ends[i, 0]  # Use direct indexing if `ends` is zero-based
        end_node = ends[i, 1]    # Use direct indexing if `ends` is zero-based
        x = [coord[start_node, 0], coord[end_node, 0]]
        y = [coord[start_node, 1], coord[end_node, 1]]
        z = [coord[start_node, 2], coord[end_node, 2]]
        ax.plot(x, y, z, color='blue')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
