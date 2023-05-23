import sympy as sy
from sympy import symbols
from matplotlib.pyplot import *
from numpy import linspace

SIZE = 1


def lagrange_poly(fun, values):
    """
    returns the function as a polynom
    :param fun: a function
    :param values: the x values of the dots we know
    :return: the polynomial of the function by lagrange
    """
    num = len(values)
    polynomial_fun = 0
    for j in range(num):
        polynomial_fun += (p_n(values, j) * fun.subs(x, values[j])) #p(x) = \sigma (
    return polynomial_fun


def p_n(values, num):
    """
    returns the p_n value
    :param values: the x values of the dots we know
    :param num: value of j
    :return: the p_n value
    """
    poly_n = 1
    for j in range(len(values)):
        if num != j:
            poly_n *= (x - values[j]) / (values[num] - values[j])
    return poly_n


x = symbols("x")
n = [10, 20, 30, 40, 50]
x_lists = []
y_lists = []
real_fun = sy.exp(-1 / (x ** 2))
fun_x = linspace(-SIZE, SIZE, 80)
fun_y = []

for i in n:
    x_lists.append(linspace(-SIZE, SIZE, i))
    y = lagrange_poly(real_fun, x_lists[-1])
    temp = []
    for k in fun_x:
        temp.append(y.subs(x, k))
    y_lists.append(temp)

for i in range(len(y_lists)):
    plot(fun_x, y_lists[i], label="poly for n = " + str(n[i]))

for i in fun_x:
    fun_y.append(real_fun.subs(x, i))

plot(fun_x, fun_y, label="exp(-1 / x ^ 2)")
xlabel("x")
ylabel("y")
grid()
legend()
show()
