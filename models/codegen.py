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


def main():
    f = zeros(14, 1)

    x = Matrix(symbols('m rx ry rz vx vy vz q0 q1 q2 q3 wx wy wz', real=True))
    u = Matrix(symbols('ux uy uz',  real=True))

    f[0, 0] = 0
  
    inputs = {'x': x, 'u': u}

    A = f.jacobian(x)
    B = f.jacobian(u)


    print("\node:")
    print(c_code_postprocessing(matrix_c_code(f, "StateVector", "f", inputs)))

    print("\nstate_jacobian:")
    print(c_code_postprocessing(matrix_c_code(A, "StateMatrix", "A", inputs)))

    print("\ncontrol_jacobian:")
    print(c_code_postprocessing(matrix_c_code(B, "ControlMatrix", "B", inputs)))


if __name__ == "__main__": main()