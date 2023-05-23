from numpy import *
from matplotlib.pyplot import *

NUM = 20


def decomposition(a):
    """
    get a matrix and return L,U,D matrix
    :param a: the A matrix
    :return:
    l_matrix: a matrix whose elements are all down-diagonal elements of matrix A
    u_matrix: a matrix whose elements are all above-diagonal elements of matrix A
    d_matrix: a matrix whose elements are all diagonal elements of matrix A
    """
    k = len(a)
    l_matrix = zeros((k, k))
    d_matrix = zeros((k, k))
    u_matrix = zeros((k, k))
    for row in range(k):
        for column in range(k):
            if row == column:
                d_matrix[row, column] = a[row, column]
            if row > column:
                u_matrix[row, column] = a[row, column]
            if row < column:
                l_matrix[row, column] = a[row, column]
    return l_matrix, d_matrix, u_matrix


def jacobi(a, x, b, n=-1):
    """
    sends to jacobi_iteration all the data its needed and return the graph
    :param a: the A matrix
    :param x: estimated vector
    :param b: solution to the equation: A*x
    :param n: number of iterations requested
    :return: array that show the relative error of X in each Jacobi iteration
    """
    l_matrix, d_matrix, u_matrix = decomposition(a)
    k = len(a)
    d_1_matrix = zeros((k, k))
    graph = [(max(abs(sul - x)) / max(abs(sul)))[0, 0]]
    for i in range(k):
        d_1_matrix[i, i] = 1 / d_matrix[i, i]
    return jacobi_iteration(l_matrix, d_1_matrix, u_matrix, x, b, n, 1, graph)


def jacobi_iteration(l_matrix, d_1_matrix, u_matrix, x, b, n, num, graph):
    """
    does the jacobi iteration until it do n iteration or gets the right solution
    :param l_matrix: a matrix whose elements are all down-diagonal elements of matrix A
    :param d_1_matrix: the inverse matrix of D matrix
    :param u_matrix: a matrix whose elements are all above-diagonal elements of matrix A
    :param x: estimated vector that approaches the correct answer as the iterations pasted
    :param b: solution to the equation: A*x
    :param n: the number of iteration that remain if it ander 0 it means the user didn't
    give n, and the function will finish as the solution will be close enough
    :param num: number of the iteration yet
    :param graph: array of the relative error of X in each Jacobi iteration
    :return: when finished return the graph
    """
    if n == 0:
        return graph
    x_1 = d_1_matrix @ b - d_1_matrix @ (l_matrix + u_matrix) @ x
    graph.append((max(abs(sul - x_1)) / max(abs(sul)))[0, 0])
    if max(abs(x - x_1)) < 0.1 ** 10:
        return graph
    if n > 0:
        print("x", num, " = ", array(x_1.T)[0], ", absolute error: ", max(abs(sul - x_1))[0, 0],
              ", relative error: ", (max(abs(sul - x_1)) / max(abs(sul)))[0, 0], sep="")
    return jacobi_iteration(l_matrix, d_1_matrix, u_matrix, x_1, b, n - 1, num + 1, graph)


def gauss(a, x, b, n=-1):
    """
    sends to gauss_iteration all the data its needed and return the graph
    :param a: the A matrix
    :param x: estimated vector
    :param b: solution to the equation: A*x
    :param n: number of iterations requested
    :return: array that show the relative error of X in each gauss-seidel iteration
    """
    l_matrix, d_matrix, u_matrix = decomposition(a)
    k = len(a)
    d_1_matrix = zeros((k, k))
    graph = [(max(abs(sul - x)) / max(abs(sul)))[0, 0]]
    for i in range(k):
        d_1_matrix[i, i] = 1 / d_matrix[i, i]
    return gauss_iteration(l_matrix, d_1_matrix, u_matrix, x, b, n, 1, graph)


def gauss_iteration(l_matrix, d_1_matrix, u_matrix, x, b, n, num, graph):
    """
    does the gauss-seidel iteration until it do n iteration or gets the right solution
    :param l_matrix: a matrix whose elements are all down-diagonal elements of matrix A
    :param d_1_matrix: the inverse matrix of D matrix
    :param u_matrix: a matrix whose elements are all above-diagonal elements of matrix A
    :param x: estimated vector that approaches the correct answer as the iterations pasted
    :param b: solution to the equation: A*x
    :param n: the number of iteration that remain if it ander 0 it means the user didn't
    give n, and the function will finish as the solution will be close enough
    :param num: number of the iteration yet
    :param graph: array of the relative error of X in each gauss-seidel iteration
    :return: when finished return the graph
    """
    if n == 0:
        return graph
    x_1 = -linalg.inv(l_matrix + linalg.inv(d_1_matrix)) @ u_matrix @ x + linalg.inv(
        l_matrix + linalg.inv(d_1_matrix)) @ b
    graph.append((max(abs(sul - x_1)) / max(abs(sul)))[0, 0])
    if max(abs(x - x_1)) < 0.1 ** 10:
        return graph
    if n > 0:
        print("x", num, " = ", array(x_1.T)[0], ", absolute error: ", max(abs(sul - x_1))[0, 0],
              ", relative error: ", (max(abs(sul - x_1)) / max(abs(sul)))[0, 0], sep="")
    return gauss_iteration(l_matrix, d_1_matrix, u_matrix, x_1, b, n - 1, num + 1, graph)


# the matrix
m = matrix([[7, -2, 1], [14, -7, -3], [-7, 11, 18]])

# solution to the equation: A*x
vector = matrix([[12], [17], [5]])

# estimated vector
sul = matrix([[3], [4], [-1]])
print("for the matrix:\n", m, "\n and the vector b:\n", vector)
print("jacobi with (0,0,0) as the guessed vector:")
jacobi_first = jacobi(m, matrix([[0], [0], [0]]), vector, NUM)
print("\njacobi with (10,10,10) as the guessed vector:")
jacobi_second = jacobi(m, matrix([[10], [10], [10]]), vector, NUM)
print("\ngauss-seidel with (0,0,0) as the guessed vector:")
gauss_first = gauss(m, matrix([[0], [0], [0]]), vector, NUM)
print("\ngauss-seidel with (10,10,10) as the guessed vector:")
gauss_second = gauss(m, matrix([[10], [10], [10]]), vector, NUM)

# and now the graph
x_k = range(NUM + 1)
xlim(0, NUM)
xlabel("numer of iterations")
ylabel("relative error")
plot(x_k, jacobi_first, label="jacobi with the guess (0, 0, 0)")
plot(x_k, jacobi_second, label="jacobi with the guess (10, 10, 10)")
plot(x_k, gauss_first, label="gauss-seidel with the guess (0, 0, 0)")
plot(x_k, gauss_second, label="gauss-seidel with the guess (10, 10, 10)")
legend()
show()
