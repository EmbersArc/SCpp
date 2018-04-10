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

        lines = []
        line_colors = []

        ax = fig.add_subplot(2,2,subplot_pos)
        ax.set_xlabel(axis_map_x(subplot_pos, 'x, up', 'y, east', 'z, north'))
        ax.set_ylabel(axis_map_y(subplot_pos, 'x, up', 'y, east', 'z, north'))

        Rx = axis_map_x(subplot_pos, X[1,:], X[2,:], X[3,:])
        Ry = axis_map_y(subplot_pos, X[1,:], X[2,:], X[3,:])

        ax.plot(Rx, Ry, color='black')
        ax.axis('equal')

        ax.set_title(['top view','','front view','side view'][subplot_pos-1])


        for k in range(K):


            rx, ry, rz = X[1:4, k]
            vx, vy, vz = X[4:7, k]
            qw, qx, qy, qz = X[7:11, k]

            CBI = np.array([
                [1-2*(qy**2+qz**2), 2*(qx*qy+qw*qz), 2*(qx*qz-qw*qy)],
                [2*(qx*qy-qw*qz), 1-2*(qx**2+qz**2), 2*(qy*qz+qw*qx)],
                [2*(qx*qz+qw*qy), 2*(qy*qz-qw*qx), 1-2*(qx**2+qy**2)]
            ])

            Fx, Fy, Fz = np.dot(np.transpose(CBI), U[:, k])
            dx, dy, dz = np.dot(np.transpose(CBI), np.array([1.,0.,0.]))


            rx_ = axis_map_x(subplot_pos, rx, ry, rz)
            ry_ = axis_map_y(subplot_pos, rx, ry, rz)

            vx_ = axis_map_x(subplot_pos, vx, vy, vz)
            vy_ = axis_map_y(subplot_pos, vx, vy, vz)

            dx_ = axis_map_x(subplot_pos, dx, dy, dz)
            dy_ = axis_map_y(subplot_pos, dx, dy, dz)

            Fx_ = axis_map_x(subplot_pos, Fx, Fy, Fz)
            Fy_ = axis_map_y(subplot_pos, Fx, Fy, Fz)

            scale = 0.5

            # speed vector
            lines.append([(rx_, ry_), (rx_ + scale * vx_, ry_ + scale * vy_)])
            line_colors.append((0, 1, 0, 1))


            # attitude vector
            lines.append([(rx_, ry_), (rx_ + scale * dx_, ry_ + scale * dy_)])
            line_colors.append((0, 0, 1, 1))

            # thrust vector
            lines.append([(rx_, ry_), (rx_ - scale * Fx_, ry_ - scale * Fy_)])
            line_colors.append((1, 0, 0, 1))

        lc = mc.LineCollection(lines, colors=line_colors, linewidths=1.5)
        ax.add_collection(lc)

    ax.axis('equal')
    fig.suptitle("iter " + str(figures_i), fontsize=14)


def main():
    fig = plt.figure(figsize=(16, 12))
    my_plot(fig, figures_i)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.show()


if __name__ == '__main__':
    main()
