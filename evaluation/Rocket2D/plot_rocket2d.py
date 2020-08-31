import numpy as np
import matplotlib.pyplot as plt
import os


figures_i = 0
figures_N = 100
FOLDER = ""

rocket_length = 20
thrust_length = 20

def my_plot(fig):
    global figures_i

    X = np.loadtxt(f"{FOLDER}/{figures_i}/X.txt", delimiter=",")
    U = np.loadtxt(f"{FOLDER}/{figures_i}/U.txt", delimiter=",")

    ax = fig.add_subplot(111)

    ax.set_xlabel('X, east')
    ax.set_ylabel('Y, up')

    rx = X[:, 0]
    ry = X[:, 1]

    dx = -np.sin(X[:, 4]) * rocket_length
    dy = np.cos(X[:, 4]) * rocket_length

    Fx = -np.sin(X[:, 4] + U[:, 0]) * U[:, 1]
    Fy = np.cos(X[:, 4] + U[:, 0]) * U[:, 1]

    max_thrust = np.max(np.sqrt(Fx**2+Fy**2))
    Fx /= max_thrust
    Fy /= max_thrust

    Fx *= thrust_length
    Fy *= thrust_length

    # position vector
    ax.plot(X[:, 0], X[:, 1], color='lightgrey', zorder=0)

    # attitude vector
    ax.quiver(rx, ry, dx, dy, color='blue', width=0.004, headwidth=1,
              headlength=0, pivot="mid", scale_units="xy", scale=1)

    # force vector
    ax.quiver(rx-dx/2, ry-dy/2, Fx, Fy, color='red', width=0.004,
              pivot="tip", scale_units="xy", scale=1, headwidth=1, headlength=0)

    ax.set_title("Iteration " + str(figures_i))
    ax.set_aspect("equal")


def key_press_event(event):
    global figures_i, figures_N

    fig = event.canvas.figure

    if event.key == 'q' or event.key == 'escape':
        plt.close(event.canvas.figure)
        return
    if event.key == 'right':
        figures_i += 1
        figures_i %= figures_N
    elif event.key == 'left':
        figures_i -= 1
        figures_i %= figures_N
    fig.clear()
    my_plot(fig)
    plt.draw()


def main():
    global figures_i, figures_N, FOLDER
    model_folder = "output/Rocket2D/SC"
    folder_num = sorted(os.listdir(model_folder))[-1]
    print(folder_num)
    FOLDER = f"{model_folder}/{folder_num}"

    figures_N = len(os.listdir(FOLDER))

    fig = plt.figure()
    figures_i = figures_N - 1
    my_plot(fig)
    cid = fig.canvas.mpl_connect('key_press_event', key_press_event)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
