import os
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

    ax.set_zlabel('X, up')
    ax.set_xlabel('Y, east')
    ax.set_ylabel('Z, north')

    K = X.shape[1]

    for k in range(K):
        rx, ry, rz = X[1:4, k]
        vx, vy, vz = X[4:7, k]
        qw, qx, qy, qz = X[7:11, k]

        CBI = np.array([
            [1 - 2 * (qy ** 2 + qz ** 2), 2 * (qx * qy + qw * qz), 2 * (qx * qz - qw * qy)],
            [2 * (qx * qy - qw * qz), 1 - 2 * (qx ** 2 + qz ** 2), 2 * (qy * qz + qw * qx)],
            [2 * (qx * qz + qw * qy), 2 * (qy * qz - qw * qx), 1 - 2 * (qx ** 2 + qy ** 2)]
        ])

        Fx, Fy, Fz = np.dot(np.transpose(CBI), U[:, k])
        dx, dy, dz = np.dot(np.transpose(CBI), np.array([1., 0., 0.]))

        # speed vector
        ax.quiver(ry, rz, rx, vy, vz, vx, length=0.25, color='green')

        # attitude vector
        ax.quiver(ry, rz, rx, dy, dz, dx, length=0.25, arrow_length_ratio=0.0, color='blue')

        # thrust vector
        ax.quiver(ry, rz, rx, -Fy, -Fz, -Fx, length=0.25, arrow_length_ratio=0.0, color='red')

    ax.axis('equal')
    ax.set_title("iter " + str(figures_i))
    ax.plot(X[2, :], X[3, :], X[1, :], color='black')


def main():
    # find the number of saved iterations
    global figures_i, figures_N
    figures_N = max([int(string.lstrip("iteration").rstrip("_X.txt"))
                     for string in os.listdir("../../output/model_landing_6dof/")
                     if string.endswith("X.txt")]) + 1
    figures_i = figures_N - 1

    fig = plt.figure(figsize=(10, 12))
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
