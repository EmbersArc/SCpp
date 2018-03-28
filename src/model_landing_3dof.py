from sympy import *
import sympy


def matrix_c_code(M, type, name, inputs):
    result = ""

    for input_name, input_vector in inputs.items():
        for i, input_sym in enumerate(input_vector):
            if sum([len(e.find(input_sym)) for e in M]) > 0:
                result += "    const double " + str(input_sym) + " = " + input_name + "(" + str(i) + ", 0);\n"

    replacements, M = sympy.cse(M)
    M = Matrix(M)

    for lhs, rhs in replacements:
        result += "    const double " + sympy.ccode(lhs) + " = " + sympy.ccode(rhs) + ";\n"

    result += "\n    " + type + " " + name + ";\n    " + name + ".setZero();\n"

    for i in range(M.rows):
        for j in range(M.cols):
            if(M[i,j] != 0):
                result += "    " + name + "(" + str(i) + ", " + str(j) + ") = "
                result += sympy.ccode(M[i,j]) + ";\n"
    return result

def c_code_postprocessing(code):
    code = code.replace("pow(rG, 2)", "(rG*rG)")
    return code

def main():

    f = sympy.zeros(6,1)

    rx, ry, vx, vy, theta, dtheta = symbols('rx ry vx vy theta dtheta')
    x = Matrix([rx, ry, vx, vy, theta, dtheta])

    throttle, gimbalAngle = symbols('throttle gimbalAngle')
    u = [throttle, gimbalAngle]

    TWR, g, rG, rTB = symbols('TWR g rG rTB')

    f[0,0] = vx
    f[1,0] = vy
    f[2,0] = g * TWR * sympy.sin(theta + gimbalAngle)
    f[3,0] = g * (TWR * sympy.cos(theta + gimbalAngle) - 1)
    f[4,0] = dtheta
    f[5,0] = -g * TWR * sympy.sin(gimbalAngle) * rTB / sympy.Mul(rG, rG, evaluate=False)

    inputs = {'x':x, 'u':u}



    A = f.jacobian(x)
    B = f.jacobian(u)



    print("\node:")
    print(c_code_postprocessing(matrix_c_code(f, "StateVector", "f", inputs)))

    print("\nstate_jacobian:")
    print(c_code_postprocessing(matrix_c_code(A, "StateMatrix", "A", inputs)))

    print("\ncontrol_jacobian:")
    print(c_code_postprocessing(matrix_c_code(B, "ControlMatrix", "B", inputs)))

    


if __name__ == "__main__": main()