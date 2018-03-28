import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
from math import cos, sin

def main():
    X = np.loadtxt("../../output/model_landing_3dof/iteration040_X.txt")
    U = np.loadtxt("../../output/model_landing_3dof/iteration040_U.txt")


    plt.figure(figsize=(10,12))
    plt.plot(X[0,:], X[1,:])

    plt.xlabel('X')
    plt.ylabel('Y')


    lines = []
    line_colors = []


    K = X.shape[1]
    
    for k in range(K):
        rx = X[0,k]
        ry = X[1,k]
        vx = X[2,k]
        vy = X[3,k]
        theta = X[4,k]
        throttle = U[0,k]
        gimbal = U[1,k]

        # speed vector
        speed_scale = 0.8
        lines.append([(rx,ry),(rx+speed_scale*vx,ry+speed_scale*vy)])
        line_colors.append((0,1,0,1))
        
        # attitude vector
        heading_scale = 6
        c_theta = heading_scale * cos(theta)
        s_theta = heading_scale * sin(theta)
        lines.append([(rx,ry),(rx+s_theta,ry+c_theta)])
        line_colors.append((0,0,1,1))

        # thrust vector
        throttle_scale = 6
        Tx = throttle_scale * throttle * sin(theta+gimbal)
        Ty = throttle_scale * throttle * cos(theta+gimbal)
        lines.append([(rx,ry),(rx-Tx,ry-Ty)])
        line_colors.append((1,0,0,1))



    lc = mc.LineCollection(lines, colors=line_colors, linewidths=1.5)

    ax = plt.gca()
    ax.add_collection(lc)
    plt.axis('equal')



    def quit_figure(event):
        if event.key == 'q' or event.key == 'escape':
            plt.close(event.canvas.figure)

    cid = plt.gcf().canvas.mpl_connect('key_press_event', quit_figure)

    plt.show()


if __name__ == '__main__': main()
    
