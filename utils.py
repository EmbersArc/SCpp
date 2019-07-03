from sympy import Matrix, ccode, cse


def matrix_c_code(M, type, name, inputs, return_value=True):
    result = "\n\n\n"

    # rename inputs
    for input_name, input_vector in inputs.items():
        for i, input_sym in enumerate(input_vector):
            if any([coefficient.find(input_sym) for coefficient in M]):
                result += f'const double {input_sym} = {input_name}({i}, 0);\n'

    # replace common subexpressions
    replacements, M = cse(M)
    M = Matrix(M)
    for lhs, rhs in replacements:
        result += f'const double {ccode(lhs)} = {ccode(rhs)};\n'

    # assemble final matrix
    result += '\n'

    if return_value:
        result += f'{type} {name};\n'

    result += f'{name}.setZero();\n\n'

    for i in range(M.rows):
        for j in range(M.cols):
            if M[i, j] != 0:
                result += f'{name}({i},{j}) = {ccode(M[i, j])};\n'

    if return_value:
        result += '\n'
        result += f'return {name};\n\n\n'

    return result
