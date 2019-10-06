import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc
from mpl_toolkits import mplot3d
from sympy import rot_axis1, rot_axis2


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

    T_max = np.max(np.linalg.norm(U, axis=1))

    for k in range(K):
        rx, ry, rz = X[k, 0:3]
        # vx, vy, vz = X[k, 3:6]
        phi, theta = X[k, 6:8]

        Rx = rot_axis1(phi).T
        Ry = rot_axis2(theta).T
        R = Rx * Ry

        Fx, Fy, Fz = np.dot(R, U[k, :] / T_max)
        dx, dy, dz = np.dot(R, np.array([0., 0., 1.]))
        # tx, ty, tz = np.dot(R, np.array([1., 0., 0.]))

        # # speed vector
        # ax.quiver(rx, ry, rz, vx, vy, vz, length=0.1, color='green')

        # attitude vector
        ax.quiver(rx, ry, rz, dx, dy, dz, length=0.5,
                  arrow_length_ratio=0.0, color='blue')

        # # up vector
        # ax.quiver(rx, ry, rz, tx, ty, tz, length=0.01,
        #           arrow_length_ratio=0.0, color='green')

        # thrust vector
        ax.quiver(rx, ry, rz, -Fx, -Fy, -Fz, length=0.5,
                  arrow_length_ratio=0.0, color='red')

    scale = 1.2 * np.abs(np.max(X[:, 0:3]))
    ax.set_xlim3d(-scale, scale)
    ax.set_ylim3d(-scale, scale)
    ax.set_zlim3d(0, scale)
    ax.plot(X[:, 0], X[:, 1], X[:, 2], color='gray')

    # fig.suptitle(f"Iteration {figures_i}", fontsize=14)
    # plt.savefig(f"{FOLDER}/{figures_i}.png")

def main():
    global figures_i, figures_N, FOLDER
    model_folder = "output/RocketHover/MPC"
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
