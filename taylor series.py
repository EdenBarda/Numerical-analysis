from sympy import *
from numpy import linspace
from matplotlib.pyplot import *


def taylor(fun, x_0, n):
    """
    the function prints the taylor taylor_series and show it on the screen.
    @param fun: the function of the taylor_series.
    @param x_0: the x_0 parameter.
    @param n: the numbers of degrees of the taylor_series.
    :return: none.
    """
    taylor_series = 0
    for i in range(n + 1):
        taylor_series += fun.diff(x, i).subs(x, x_0) * (x - x_0) ** i / factorial(i)
    print(taylor_series)

    # making x and y points of the taylor_series.
    x_list = linspace(-20, 20, 100)
    y_list = []
    for i in x_list:
        y_list.append(taylor_series.subs(x, i))

    # create the taylor_series on the screen.
    xlim(-20, 20)
    ylim(-20, 20)
    plot(x_list, y_list, label="abc")
    xlabel("x")
    ylabel("y")
    grid()
    legend()
    show()


x = symbols("x")
taylor(sin(x), 0, 20)
