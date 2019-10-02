import numpy as np
import os

figures_N = sum(f.endswith("X.txt")
                for f in os.listdir("output/RocketLanding3D/"))

iteration0 = str(figures_N-1).zfill(3)
iteration1 = str(figures_N-2).zfill(3)

X0 = np.loadtxt(f"output/RocketLanding3D/iteration{iteration0}_X.txt")
U0 = np.loadtxt(f"output/RocketLanding3D/iteration{iteration0}_U.txt")
X1 = np.loadtxt(f"output/RocketLanding3D/iteration{iteration1}_X.txt")
U1 = np.loadtxt(f"output/RocketLanding3D/iteration{iteration1}_U.txt")
K = X0.shape[1]


print(np.linalg.norm(X0 - X1, 1))
print(np.linalg.norm(U0 - U1, 1))
