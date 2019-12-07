import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline, BSpline

model_folder = "output/Starship/SC"
folder_num = sorted(os.listdir(model_folder))[-1]
folder = f"{model_folder}/{folder_num}"
acc_b = np.loadtxt(f"{folder}/x/acc_passenger_b.txt", delimiter=",")

acc_b = acc_b[1:-1, :]
K = acc_b.shape[0]

g = 9.81
g_b = acc_b / g

xnew = np.linspace(0, K, 300)
spl = make_interp_spline(np.arange(0, K), g_b[:, 0], k=3)
smooth_x = spl(xnew)
spl = make_interp_spline(np.arange(0, K), g_b[:, 1], k=3)
smooth_y = spl(xnew)
spl = make_interp_spline(np.arange(0, K), g_b[:, 2], k=3)
smooth_z = spl(xnew)
spl = make_interp_spline(
    np.arange(0, K), np.linalg.norm(g_b, axis=1), k=3)
smooth_g = spl(xnew)


fig = plt.figure(figsize=(15, 5))
ax = fig.add_subplot(1, 1, 1)
ax.plot(smooth_x, color="r", label="x")
ax.plot(smooth_y, color="g", label="y")
ax.plot(smooth_z, color="b", label="z")
ax.plot(smooth_g, color="black", label="norm")
ax.hlines(0, xmin=0, xmax=xnew.shape[0], linestyles="--", color="gray")
# ax.hlines(1, xmin=0, xmax=xnew.shape[0], linestyles="--", color="gray")
ax.set_xlim([0, xnew.shape[0]])
ax.legend()
plt.tight_layout()
plt.show()
