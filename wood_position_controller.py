%reset

# %%
from sympy.physics.vector import dynamicsymbols
from sympy import *
import numpy as np

from utils import matrix_c_code

# %%

# %%

t = symbols('t')

# Trajectory
# omega_hat = dynamicsymbols('omega_hat')
# T_f_ = symbols('T_f_', real=True)
# r_x_, r_y_, r_z_ = symbols('r_x_, r_y_, r_z_', real=True)
# omega_hat = pi * (1. - cos(pi * t/T_f_))
# lambda_d = Matrix([r_x_ * cos(omega_hat),
#                    r_y_ * sin(omega_hat),
#                    r_z_ * sin(omega_hat)])
# dot_lambda_d = simplify(diff(lambda_d, t, 1))
# dot_dot_lambda_d = simplify(diff(lambda_d, t, 2))
# dot_dot_dot_lambda_d = simplify(diff(lambda_d, t, 3))


# Orientation, rotation matrix and its Jacobian
phi_sym, theta_sym, psi_sym = dynamicsymbols('phi, theta, psi')
eta_sym = Matrix([phi_sym, theta_sym, psi_sym])

lamda_sym = Matrix(dynamicsymbols('lamda:3'))
sigma_sym = Matrix(dynamicsymbols('sigma:3'))
omega_sym = Matrix(dynamicsymbols('omega:3'))
k_1_, k_2_, k_3_, k_4_, kg_, km_, kI1_, kl_ = symbols(
    'k_1_, k_2_, k_3_, k_4_, kg_, km_, kI1_, kl_', real=True)


def R(v):
    return Matrix([[cos(v[1]) * cos(v[2]),
                    -cos(v[1]) * sin(v[2]),
                    sin(v[1])],
                   [sin(v[0]) * sin(v[1]) * cos(v[2]) + cos(v[0]) * sin(v[2]),
                    cos(v[0]) * cos(v[2]) - sin(v[0]) * sin(v[1]) * sin(v[2]),
                    -sin(v[0]) * cos(v[1])],
                   [-cos(v[0]) * sin(v[1]) * cos(v[2]) + sin(v[0]) * sin(v[2]),
                    sin(v[0]) * cos(v[2]) + cos(v[0]) * sin(v[1]) * sin(v[2]),
                    cos(v[0]) * cos(v[1])]])


def T(v):
    return 1. / cos(v[1]) * Matrix([[cos(v[2]), -sin(v[2]), 0.],
                                    [cos(v[1]) * sin(v[2]),
                                     cos(v[1]) * cos(v[2]), 0],
                                    [-sin(v[1]) * cos(v[2]),
                                     sin(v[1]) * sin(v[2]), cos(v[1])]])


E_x = Matrix([1, 0, 0])
E_y = Matrix([0, 1, 0])
E_z = Matrix([0, 0, 1])

# e_x = R(eta_sym) * E_x
# e_y = R(eta_sym) * E_y
# e_z = R(eta_sym) * E_z

I3x2 = Matrix([[1, 0], [0, 1], [0, 0]])
I2x3 = Matrix([[1, 0, 0], [0, 1, 0]])

# symbols to be replaced later
a_d_sym = Matrix(dynamicsymbols('a_d:3'))

# desired orientation and its derivatives
phi_d = atan2(-a_d_sym[1], (a_d_sym[2] + kg_))
theta_d = atan2(a_d_sym[0], sqrt(a_d_sym[1]**2 + (a_d_sym[2] + kg_)**2))
eta_d = Matrix([phi_d, theta_d, 0])
dot_eta_d = simplify(eta_d.diff(t))
dot_dot_eta_d = simplify(dot_eta_d.diff(t))

# a_d, dot_a_d and dot_dot_a_d to replace symbols later
a_d = -(k_1_ + k_2_) * sigma_sym - (k_1_ * k_2_ + 1) * lamda_sym

# singularity
# a_d[2] += 0.01 * Piecewise((1., (a_d + kg_ * E_z).norm() >= 0.98 * kg_),
#                            (0., True))
# a_d[2] = 0.9 * kg_ * tanh(a_d[2])

# Normal control law:
bar_T_z = km_ * sqrt(a_d_sym[0]**2 + a_d_sym[1]**2 + (a_d_sym[2] + kg_)**2)
dot_bar_T_z = bar_T_z.diff(t)

# Improved control law:
# condition = (R(eta_sym) * E_z).dot(R(eta_d) * E_z) >= 0  # switching condition
# bar_T_z = Piecewise(((R(eta_sym) * E_z).dot(a_d + kg_ * E_z), condition),
#                     (0., True))
a = 1. / km_ * bar_T_z * R(eta_sym) * E_z - kg_ * E_z
dot_a_d = -(k_1_ + k_2_) * sigma_sym - (k_1_ * k_2_ + 1) * a
# dot_bar_T_z = Piecewise(((R(eta_sym) * E_z).dot(a_d + kg_ * E_z) +
#                          (R(eta_sym) * E_z).dot(dot_a_d), condition),
#                         (0., True))

dot_a = 1. / km_ * dot_bar_T_z * R(eta_sym) * E_z \
    + 1. / km_ * bar_T_z * R(eta_sym) * omega_sym.cross(E_z)
dot_dot_a_d = -(k_1_ + k_2_) * dot_a - (k_1_ * k_2_ + 1) * a

T_inv = simplify(T(eta_sym).inv())
TyTx = kI1_ / kl_ * I2x3 * T_inv * \
    (
        dot_dot_eta_d
        - T(eta_sym).diff(t) * omega_sym
        - (k_3_ + k_4_) * (T(eta_sym) * omega_sym - dot_eta_d)
        - (k_3_ * k_4_ + 1) * (eta_sym - eta_d)
    )

# to replace dot_eta_sym
dot_eta = T(eta_sym) * omega_sym

subst = [
    *zip(a_d_sym.diff(t, 2), dot_dot_a_d),
    *zip(a_d_sym.diff(t, 1), dot_a_d),
    *zip(a_d_sym, a_d),
    *zip(eta_sym.diff(t), dot_eta),
    *zip(lamda_sym, Matrix(symbols('lambda:3', real=True))),
    *zip(sigma_sym, Matrix(symbols('sigma:3', real=True))),
    *zip(eta_sym, Matrix(symbols('eta:3', real=True))),
    *zip(omega_sym, Matrix(symbols('omega:3', real=True))),
]
control_vector = Matrix([-TyTx[1], TyTx[0], bar_T_z]).subs(subst)

# init_printing()

# %%

inputs = {
    'lambda': Matrix(symbols('lambda:3', real=True)),
    'sigma': Matrix(symbols('sigma:3', real=True)),
    'eta': Matrix(symbols('eta:3', real=True)),
    'omega': Matrix(symbols('omega:3', real=True)),
}

print(matrix_c_code(control_vector, 'Eigen::Vector3d', 'control_vector', inputs))

# %%
