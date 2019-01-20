import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc

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
    ax = fig.add_subplot(111)

    iteration = str(figures_i).zfill(3)

    X = np.loadtxt(f"output/RocketLanding2D/iteration{iteration}_X.txt")
    U = np.loadtxt(f"output/RocketLanding2D/iteration{iteration}_U.txt")

    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    lines = []
    line_colors = []

    K = X.shape[1]

    for k in range(K):
        rx = X[0, k]
        ry = X[1, k]
        vx = X[2, k]
        vy = X[3, k]
        theta = X[4, k]
        throttle = U[0, k]
        gimbal = U[1, k]

        # speed vector
        speed_scale = 0.8
        lines.append(
            [(rx, ry), (rx + speed_scale * vx, ry + speed_scale * vy)])
        line_colors.append((0, 1, 0, 1))

        # attitude vector
        heading_scale = 0.1
        c_theta = heading_scale * np.cos(theta)
        s_theta = heading_scale * np.sin(theta)
        lines.append([(rx, ry), (rx + s_theta, ry + c_theta)])
        line_colors.append((0, 0, 1, 1))

        # thrust vector
        throttle_scale = 0.1
        Tx = throttle_scale * throttle * np.sin(theta + gimbal)
        Ty = throttle_scale * throttle * np.cos(theta + gimbal)
        lines.append([(rx, ry), (rx - Tx, ry - Ty)])
        line_colors.append((1, 0, 0, 1))

    lc = mc.LineCollection(lines, colors=line_colors, linewidths=1.5)

    ax.add_collection(lc)
    ax.axis('equal')
    ax.set_title("iter " + str(figures_i))


def main():
    global figures_i, figures_N
    figures_N = sum(f.endswith("X.txt")
                    for f in os.listdir("output/RocketLanding3D/"))
    fig = plt.figure(figsize=(10, 10))
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
