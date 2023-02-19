"assignment 2"

import numpy as np
np.set_printoptions(precision=7, suppress=True, linewidth=100)

q1_coordinates = np.array([[3.6, 1.675, -1.195], [3.8, 1.436, -1.188], [3.9, 1.318, -1.182]])
q1_x_val = 3.7
q2_coordinates = np.array([[7.2, 23.5492], [7.4, 25.3913], [7.5, 26.8224], [7.6, 27.4589]])
q3_x_val = 7.3
q5_coordinates = np.array([[2 ,3], [5, 5], [8, 7], [10, 9]])

def neville(raw, x_val):
    "performs neville's method"

    output = np.zeros((len(raw), len(raw) + 1))
    for i in range(len(raw)):
        output[(i, 0)] = raw[(i, 0)]
        output[(i, 1)] = raw[(i, 1)]

    for i in range(2, output.shape[1]):
        for j in range(i - 1, output.shape[0]):
            neville_recursion(output, i, j, x_val)

    print(output[output.shape[0] - 1, output.shape[1] - 1])

def neville_recursion(raw, i_val, j_val, x_val):
    "Performs the recursive method"
    x0_val = raw[(j_val - i_val + 1, 0)]
    x1_val = raw[(j_val, 0)]
    y0_val = raw[(j_val - 1, i_val - 1)]
    y1_val = raw[(j_val, i_val - 1)]
    p_x = "(1 / (x1_val - x0_val)) * (((x_val - x0_val) * y1_val) - ((x_val - x1_val) * y0_val))"
    raw[(j_val, i_val)] = eval(p_x)

def newtons_forward_diff(raw, x_val):
    "performs neville's method"

    output = np.zeros((len(raw), len(raw) + 1))
    for i in range(raw.shape[0]):
        output[(i, 0)] = raw[(i, 0)]
        output[(i, 1)] = raw[(i, 1)]

    for i in range(2, output.shape[1]):
        for j in range(i - 1, output.shape[0]):
            forward_diff(output, i, j)
            if j == i - 1:
                print(output[(j, i)])

    print()
    newtons_forward_diff_eval(output, x_val)

def forward_diff(raw, i_val, j_val):
    "Gets the next forward difference in our matrix"
    x0_val = raw[(j_val - i_val + 1, 0)]
    x1_val = raw[(j_val, 0)]
    y0_val = raw[(j_val - 1, i_val - 1)]
    y1_val = raw[(j_val, i_val - 1)]
    p_x = "(y1_val - y0_val) / (x1_val - x0_val)"
    raw[(j_val, i_val)] = eval(p_x)

def newtons_forward_diff_eval(raw, x_val):
    "getting the forward difference"
    aprrox = 0
    approx = raw[(0, 1)]
    for i in range(3):
        temp = raw[(i + 1, i + 2)]
        for j in range(i, -1, -1):
            temp *= (x_val - raw[(j, 0)])
        approx += temp
    print(approx)

def hermite_approx_matrix(raw):
    "hermite's method"
    output = np.zeros((raw.shape[0] * 2, raw.shape[0] * 2 + 1))
    for i in range(raw.shape[0] * 2):
        output[(i, 0)] = raw[(int(i / 2), 0)]
        output[(i, 1)] = raw[(int(i / 2), 1)]
        if i & 1:
            output[(i, 2)] = raw[(int(i / 2), 2)]
        elif i != 0:
            output[i, 2] = ((output[(i, 1)] - output[(i - 1, 1)]) /
            (output[(i, 0)] - output[(i - 1, 0)]))
    for i in range(3, output.shape[1]):
        for j in range(i - 1, output.shape[0]):
            forward_diff(output, i, j)
    print(output)

def cubic_spline(raw):
    "performs cubic spline"
    output = np.zeros((raw.shape[0], raw.shape[0]))
    output[(0,0)] = 1
    output[(output.shape[0] - 1, output.shape[1] - 1)] = 1
    for i in range(1, output.shape[0] - 1):
        output[(i, i - 1)] = raw[(i, 0)] - raw[(i - 1, 0)]
        output[(i, i + 1)] = raw[(i + 1, 0)] - raw[(i, 0)]
        output[(i, i)] = 2 * (output[(i, i - 1)] + output[(i, i + 1)])
    print(output)
    vector_b = np.zeros(output.shape[0])
    for i in range(1, vector_b.shape[0] - 1):
        vector_b[i] = ((3 / (raw[(i + 1, 0)] - raw[(i, 0)])) * (raw[(i + 1, 1)] - raw[(i, 1)]) -
        (3 / (raw[(i, 0)] - raw[(i - 1, 0)])) * (raw[(i, 1)] - raw[(i - 1, 1)]))
    print(vector_b)
    print(np.linalg.solve(output, vector_b))

if __name__ == "__main__":
    neville(q1_coordinates, q1_x_val)
    print()
    newtons_forward_diff(q2_coordinates, q3_x_val)
    print()
    hermite_approx_matrix(q1_coordinates)
    print()
    cubic_spline(q5_coordinates)
