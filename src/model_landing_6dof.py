from sympy import *


def matrix_c_code(M, type, name, inputs):
    result = ""

    for input_name, input_vector in inputs.items():
        for i, input_sym in enumerate(input_vector):
            if sum([len(e.find(input_sym)) for e in M]) > 0:
                result += "    const double " + str(input_sym) + " = " + input_name + "(" + str(i) + ", 0);\n"

    replacements, M = cse(M)
    M = Matrix(M)

    for lhs, rhs in replacements:
        result += "    const double " + ccode(lhs) + " = " + ccode(rhs) + ";\n"

    result += "\n    " + type + " " + name + ";\n    " + name + ".setZero();\n"

    for i in range(M.rows):
        for j in range(M.cols):
            if M[i, j] != 0:
                result += "    " + name + "(" + str(i) + ", " + str(j) + ") = "
                result += ccode(M[i, j]) + ";\n"
    return result


def c_code_postprocessing(code):

    return code


def skew(v):
    return Matrix([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])


def dir_cosine(q):
    return Matrix([
        [1 - 2 * (q[2] ** 2 + q[3] ** 2), 2 * (q[1] * q[2] + q[0] * q[3]), 2 * (q[1] * q[3] - q[0] * q[2])],
        [2 * (q[1] * q[2] - q[0] * q[3]), 1 - 2 * (q[1] ** 2 + q[3] ** 2), 2 * (q[2] * q[3] + q[0] * q[1])],
        [2 * (q[1] * q[3] + q[0] * q[2]), 2 * (q[2] * q[3] - q[0] * q[1]), 1 - 2 * (q[1] ** 2 + q[2] ** 2)]
    ])


def omega(w):
    return Matrix([
        [0, -w[0], -w[1], -w[2]],
        [w[0], 0, w[2], -w[1]],
        [w[1], -w[2], 0, w[0]],
        [w[2], w[1], -w[0], 0],
    ])

def main():
    f = zeros(14, 1)

    x = Matrix(symbols('m rx ry rz vx vy vz q0 q1 q2 q3 wx wy wz', real=True))
    u = Matrix(symbols('ux uy uz',  real=True))

    g_I = Matrix(symbols('g_I:3'))
    r_T_B = Matrix(symbols('r_T_B:3'))
    J_B = symbols('J_B:3')

    alpha_m = symbols('alpha_m')

    C_B_I = dir_cosine(x[7:11, 0])
    C_I_B = C_B_I.transpose()

    f[0, 0] = - alpha_m * u.norm()
    f[1:4, 0] = x[4:7, 0]
    f[4:7, 0] = 1 / x[0, 0] * C_I_B * u + g_I
    f[7:11, 0] = 1 / 2 * omega(x[11:14, 0]) * x[7: 11, 0]
    f[11:14, 0] = Matrix.diag(J_B) ** -1 * (skew(r_T_B) * u - skew(x[11:14, 0]) * Matrix.diag(J_B) * x[11:14, 0])

    inputs = {'x': x, 'u': u, 'g_I': g_I, 'r_T_B': r_T_B, 'J_B': J_B}

    A = f.jacobian(x)
    B = f.jacobian(u)

    print(f)

    print(B)

    print("\node:")
    print(c_code_postprocessing(matrix_c_code(f, "StateVector", "f", inputs)))

    print("\nstate_jacobian:")
    print(c_code_postprocessing(matrix_c_code(A, "StateMatrix", "A", inputs)))

    print("\ncontrol_jacobian:")
    print(c_code_postprocessing(matrix_c_code(B, "ControlMatrix", "B", inputs)))


if __name__ == "__main__": main()
