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


def axis_map_x(subplot_pos, val_up, val_east, val_north):
    if subplot_pos == 1: # top view (east, north)
        return val_east
    if subplot_pos == 3: # front view (east, up)
        return val_east
    if subplot_pos == 4: # side view (north, up)
        return val_north


def axis_map_y(subplot_pos, val_up, val_east, val_north):
    if subplot_pos == 1: # top view (east, north)
        return val_north
    if subplot_pos == 3: # front view (east, up)
        return val_up
    if subplot_pos == 4: # side view (north, up)
        return val_up


def my_plot(fig, figures_i):

    iteration = str(figures_i).zfill(3)

    X = np.loadtxt("../../output/model_landing_6dof/iteration" + iteration + "_X.txt")
    U = np.loadtxt("../../output/model_landing_6dof/iteration" + iteration + "_U.txt")


    K = X.shape[1]


    for subplot_pos in [1,3,4]:

        ax = fig.add_subplot(2,2,subplot_pos)
        ax.set_xlabel(axis_map_x(subplot_pos, 'x, up', 'y, east', 'z, north'))
        ax.set_ylabel(axis_map_y(subplot_pos, 'x, up', 'y, east', 'z, north'))

        rx = axis_map_x(subplot_pos, X[1,:], X[2,:], X[3,:])
        ry = axis_map_y(subplot_pos, X[1,:], X[2,:], X[3,:])

        ax.plot(rx, ry, color='black')
        ax.axis('equal')

        ax.set_title(['top view','','front view','side view'][subplot_pos-1])


        for k in range(K):

            rx = axis_map_x(subplot_pos, X[1,k], X[2,k], X[3,k])
            ry = axis_map_y(subplot_pos, X[1,k], X[2,k], X[3,k])

            vx = axis_map_x(subplot_pos, X[4,k], X[5,k], X[6,k])
            vy = axis_map_y(subplot_pos, X[4,k], X[5,k], X[6,k])

            ux = axis_map_x(subplot_pos, U[0,k], U[1,k], U[2,k])
            uy = axis_map_y(subplot_pos, U[0,k], U[1,k], U[2,k])


            # speed vector
            ax.quiver(rx, ry, vx, vy, scale=5, color='green', width=0.002)

            # # attitude vector
            w, x, y, z = X[7:11, k]
            dx_ = 2 * (x * z - w * y)
            dy_ = 2 * (y * z + w * x)
            dz_ = 1 - 2 * (x * x + y * y)

            dx = axis_map_x(subplot_pos, dx_, dy_, dz_)
            dy = axis_map_y(subplot_pos, dx_, dy_, dz_)

            ax.quiver(rx, ry, dx, dy, scale=5, color='blue', width=0.002)

            # thrust vector
            ax.quiver(rx, ry, -ux, -uy, scale=5, color='red', width=0.002)

    ax.axis('equal')
    fig.suptitle("iter " + str(figures_i), fontsize=14)


def main():
    fig = plt.figure(figsize=(16, 12))
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
