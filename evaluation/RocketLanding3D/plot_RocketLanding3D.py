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
    iteration = str(figures_i).zfill(3)

    X = np.loadtxt(f"output/RocketLanding3D/iteration{iteration}_X.txt")
    U = np.loadtxt(f"output/RocketLanding3D/iteration{iteration}_U.txt")

    K = X.shape[1]

    # 3D
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_xlabel('X, east')
    ax.set_ylabel('Y, north')
    ax.set_zlabel('Z, up')

    for k in range(K):
        rx, ry, rz = X[1:4, k]
        # vx, vy, vz = X[4:7, k]
        qw, qx, qy, qz = X[7:11, k]

        CBI = np.array([
            [1 - 2 * (qy ** 2 + qz ** 2),
             2 * (qx * qy + qw * qz),
             2 * (qx * qz - qw * qy)],
            [2 * (qx * qy - qw * qz),
             1 - 2 * (qx ** 2 + qz ** 2),
             2 * (qy * qz + qw * qx)],
            [2 * (qx * qz + qw * qy),
             2 * (qy * qz - qw * qx),
             1 - 2 * (qx ** 2 + qy ** 2)]
        ])

        Fx, Fy, Fz = np.dot(np.transpose(CBI), U[0:3, k])
        dx, dy, dz = np.dot(np.transpose(CBI), np.array([0., 0., 1.]))
        tx, ty, tz = np.dot(np.transpose(CBI), np.array([1., 0., 0.]))

        # # speed vector
        # ax.quiver(rx, ry, rz, vx, vy, vz, length=0.1, color='green')

        # attitude vector
        ax.quiver(rx, ry, rz, dx, dy, dz, length=0.04,
                  arrow_length_ratio=0.0, color='blue')

        # up vector
        ax.quiver(rx, ry, rz, tx, ty, tz, length=0.01,
                  arrow_length_ratio=0.0, color='green')

        # thrust vector
        ax.quiver(rx, ry, rz, -Fx, -Fy, -Fz, length=0.1,
                  arrow_length_ratio=0.0, color='red')

    scale = np.max(X[1:4, :])
    ax.set_xlim3d(-scale/2, scale/2)
    ax.set_ylim3d(-scale/2, scale/2)
    ax.set_zlim3d(0, scale)
    ax.set_aspect('equal')
    ax.plot(X[1, :], X[2, :], X[3, :], color='black')

    fig.suptitle("iter " + str(figures_i), fontsize=14)


def main():
    global figures_i, figures_N
    figures_N = sum(f.endswith("X.txt")
                    for f in os.listdir("output/RocketLanding3D/"))

    fig = plt.figure(figsize=(10, 10))
    figures_i = figures_N - 1
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
