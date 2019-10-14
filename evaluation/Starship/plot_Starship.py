import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from mpl_toolkits import mplot3d


figures_i = 0
figures_N = 100
FOLDER = ""


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
    my_plot(fig)
    plt.draw()


def my_plot(fig):
    global figures_i
    X = np.loadtxt(f"{FOLDER}/{figures_i}/X.txt", delimiter=",")
    U = np.loadtxt(f"{FOLDER}/{figures_i}/U.txt", delimiter=",")

    K = X.shape[0]

    # 3D
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_xlabel('X, east')
    ax.set_ylabel('Y, north')
    ax.set_zlabel('Z, up')

    for k in range(K):
        rx, ry, rz = X[k, 1:4]
        # vx, vy, vz = X[k, 4:7]
        qw, qx, qy, qz = X[k, 7:11]

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
        CIB = np.transpose(CBI)

        Fx, Fy, Fz = np.dot(CIB, U[k, :])
        dx, dy, dz = np.dot(CIB, np.array([0., 0., 1.]))
        # tx, ty, tz = np.dot(CIB, np.array([1., 0., 0.]))

        # # speed vector
        # ax.quiver(rx, ry, rz, vx, vy, vz, length=0.1, color='green')

        # attitude vector
        ax.quiver(rx, ry, rz, dx, dy, dz, length=30,
                  arrow_length_ratio=0.0, color='blue')

        # # up vector
        # ax.quiver(rx, ry, rz, tx, ty, tz, length=0.01,
        #           arrow_length_ratio=0.0, color='green')

        # thrust vector
        ax.quiver(rx, ry, rz, -Fx, -Fy, -Fz, length=1/30000,
                  arrow_length_ratio=0.0, color='red')

    scale = np.abs(np.max(X[:, 1:4]))
    ax.set_xlim3d(-scale/2, scale/2)
    ax.set_ylim3d(-scale/2, scale/2)
    ax.set_zlim3d(0, scale)
    # ax.set_aspect('equal')
    ax.plot(X[:, 1], X[:, 2], X[:, 3], color='gray')

    fig.suptitle("iter " + str(figures_i), fontsize=14)


def main():
    global figures_i, figures_N, FOLDER
    model_folder = "output/Starship/SC"
    folder_num = sorted(map(int, os.listdir(model_folder)))[-1]
    print(folder_num)
    FOLDER = f"{model_folder}/{folder_num}"

    figures_N = len(os.listdir(FOLDER))

    fig = plt.figure(figsize=(15, 15))
    figures_i = figures_N - 1
    my_plot(fig)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
