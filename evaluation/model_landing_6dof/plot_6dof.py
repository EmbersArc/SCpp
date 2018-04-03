import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from mpl_toolkits import mplot3d


figures_i = 0
figures_N = 40


def key_press_event(event):
    global figures_i
    fig = event.canvas.figure

    if event.key == 'q' or event.key == 'escape':
        plt.close(event.canvas.figure)
        return

    if event.key == 'right':
        figures_i = (figures_i + 1) % figures_N
    elif event.key == 'left':
        figures_i = (figures_i - 1) % figures_N

    fig.clear()
    my_plot(fig, figures_i)
    plt.draw()


def my_plot(fig, figures_i):
    ax = fig.add_subplot(111, projection='3d')

    iteration = str(figures_i).zfill(3)

    X = np.loadtxt("../../output/model_landing_6dof/iteration" + iteration + "_X.txt")
    U = np.loadtxt("../../output/model_landing_6dof/iteration" + iteration + "_U.txt")

    ax.set_xlabel('Z')
    ax.set_ylabel('Y')
    ax.set_zlabel('X')

    K = X.shape[1]

    for k in range(K):
        rz, ry, rx = X[1:4, k]
        vz, vy, vx = X[4:7, k]

        # speed vector
        ax.quiver(rx, ry, rz, vx, vy, vz, length=0.25)

        # # attitude vector
        # TODO

        # thrust vector
        # TODO

    ax.axis('equal')
    ax.set_title("iter " + str(figures_i))
    ax.plot(X[3, :], X[2, :], X[1, :], color='r')

def main():
    fig = plt.figure(figsize=(10, 12))
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
