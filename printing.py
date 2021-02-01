import sympy as sym

def print_matrix(M):
    strings = [[str(sym.simplify(c)) for c in M.row(i)] for i in range(M.rows)]
    for col in range(M.cols-1):
        max_len = 0
        for row in range(M.rows):
            if len(strings[row][col]) > max_len:
                max_len = len(strings[row][col])
        for row in range(M.rows):
            l = len(strings[row][col])
            strings[row][col] += ",   " + " "*max(0, max_len-l)
    max_len = max(len("".join(row)) for row in strings)
    print("-" * max_len)
    for row in strings:
        print("".join(row))
    print("-" * max_len)
        


def print_coeffs(poly):
    print(", ".join(str(c) for c in poly.all_coeffs()[::-1]))
    
