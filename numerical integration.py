"""
Eden Barda
208932202
"""

from matplotlib.pyplot import *
import numpy as np

# approximate solution of the integral from 0 to 2*pi of e^cos(x)dx
SOLUTION = 7.954926521012

# max numbers of samples
NUM = 20


def trapezoidal_rule(start, end, fun, n):
    """
    do integration by trapezoidal rule
    :param start: the start bound of the integral
    :param end: the end bound of the integral
    :param fun: the function
    :param n: numbers of samples
    :return: the approximate integral, his absolute error, and his relative error
    """
    delta_x = (end - start) / n
    x_n = start
    x_n_1 = start + delta_x
    integral = 0
    for j in range(n):
        integral += (fun(x_n) + fun(x_n_1)) / 2 * delta_x
        x_n = x_n_1
        x_n_1 += delta_x
    integral = float(integral)
    absolute_error = abs(SOLUTION - integral)
    relative_error = absolute_error / abs(SOLUTION)
    return integral, absolute_error, relative_error


def simpson_rule(start, end, fun, n):
    """
    do integration by simpson rule
    :param start: the start bound of the integral
    :param end: the end bound of the integral
    :param fun: the function
    :param n: numbers of samples
    :return: the approximate integral, his absolute error, and his relative error
    """
    delta_x = (end - start) / n
    x_n = start
    x_n_1 = start + delta_x
    integral = 0
    for j in range(n):
        integral += (fun(x_n) + fun(x_n_1) + 4 * fun(x_n_1 + x_n)) / 6 * delta_x
        x_n = x_n_1
        x_n_1 += delta_x
    integral = float(integral)
    absolute_error = abs(SOLUTION - integral)
    relative_error = absolute_error / abs(SOLUTION)
    return integral, absolute_error, relative_error


def gaussian_quadrature(start, end, fun, n):
    """
    do integration by gaussian quadrature rule
    :param start: the start bound of the integral
    :param end: the end bound of the integral
    :param fun: the function
    :param n: numbers of samples
    :return: the approximate integral, his absolute error, and his relative error
    """
    xi, ai = np.polynomial.legendre.leggauss(n)
    xi = (end - start) / 2 * xi + (end + start) / 2
    ai = (end - start) / 2 * ai
    integral = sum(ai * fun(xi))
    absolute_error = abs(SOLUTION - integral)
    relative_error = absolute_error / abs(SOLUTION)
    return integral, absolute_error, relative_error


def f(x):
    """
    do the function e^cos(x) on variable
    :param x: variable
    :return: the f(x)
    """
    return np.exp(np.cos(x))


a = 0
b = 2 * np.pi
trapezoidal = []
simpson = []
gauss = []
x_value = range(NUM)

print("for trapezoidal rule:")
for i in range(NUM):
    approx, absolute, relative = trapezoidal_rule(a, b, f, i + 1)
    trapezoidal.append(relative)
    print(i + 1, ": approx integral: ", approx, ", absolute error: ", absolute, ", relative error: ", relative, sep="")

print("\nfor simpson rule:")
for i in range(NUM):
    approx, absolute, relative = simpson_rule(a, b, f, i + 1)
    simpson.append(relative)
    print(i + 1, ": approx integral: ", approx, ", absolute error: ", absolute, ", relative error: ", relative, sep="")

print("\nfor gaussian quadrature rule:")
for i in range(NUM):
    approx, absolute, relative = gaussian_quadrature(a, b, f, i + 1)
    gauss.append(relative)
    print(i + 1, ": approx integral: ", approx, ", absolute error: ", absolute, ", relative error: ", relative, sep="")

# prints the graph of the relative error
xlabel("numer of samples")
ylabel("relative error")
plot(x_value, trapezoidal, label="trapezoidal")
plot(x_value, simpson, label="simpson")
plot(x_value, gauss, label="gauss")
legend()
show()
