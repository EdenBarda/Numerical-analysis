"""
Eden Barda
208932202
"""


def bisection_method(n):
    """
    finds sqrt(n) with bisection method.
    :param n: a number that we need to find is root.
    :return: sqrt(n).
    """
    start = 1
    end = n
    i = 1
    while end - start > 10 ** -6:
        middle = (end + start) / 2
        print("step:", i, "Interval: [{}, {}]".format(start, end), "f_start:", start ** 2 - n, "f_end:",
              end ** 2 - n, "absolute error:", abs(n ** 0.5 - start), "relative error:",
              abs((n ** 0.5 - start)) / n ** 0.5)
        if (start ** 2 - n) * (middle ** 2 - n) > 0:
            start = middle
        else:
            end = middle
        i += 1
    print("step:", i, "Interval: [{}, {}]".format(start, end), "f_start:", start ** 2 - n, "f_end:",
          end ** 2 - n, "absolute error:", abs(n ** 0.5 - start), "relative error:",
          abs((n ** 0.5 - start)) / n ** 0.5)
    return start


def newton_raphson_method(n):
    """
    finds sqrt(n) with newton_raphson method.
    :param n: a number that we need to find is root.
    :return: sqrt(n).
    """
    a = n / 2
    b = a - (a - n / a) / 2
    i = 1
    while abs(b - a) > 10 ** -6:
        print("step:", i, "x_i:", a, "f(x_i):", a ** 2 - n, "f'(x_i):", 2 * a, "absolute error:", abs(n ** 0.5 - a),
              "relative error:", abs((n ** 0.5 - a)) / n ** 0.5)
        temp = b
        b = b - (b - n / b) / 2
        a = temp
        i += 1
    print("step:", i, "x_i:", b, "f(x_i):", b ** 2 - n, "f'(x_i):", 2 * b, "absolute error:", abs(n ** 0.5 - b),
          "relative error:", abs((n ** 0.5 - b)) / n ** 0.5)
    return b


def regula_falsi_method(n):
    """
    finds sqrt(n) with regula falsi method.
    :param n: a number that we need to find is root.
    :return: sqrt(n).
    """
    a = 1
    b = c = n
    i = 1
    while b - a > 10 ** -6 and abs(c ** 2 - n) > 10 ** -6:
        c = a - (a ** 2 - n) * (b - a) / (b ** 2 - a ** 2)
        print("step:", i, "Interval: [{}, {}]".format(a, b), "f_a:", a ** 2 - n, "f_b:",
              b ** 2 - n, "absolute error:", abs(n ** 0.5 - a), "relative error:",
              abs((n ** 0.5 - a)) / n ** 0.5)
        if (a ** 2 - n) * (c ** 2 - n) > 0:
            a = c
        else:
            b = c
        i += 1
    print("step:", i, "Interval: [{}, {}]".format(a, b), "f_a:", a ** 2 - n, "f_b:",
          b ** 2 - n, "absolute error:", abs(n ** 0.5 - c), "relative error:",
          abs((n ** 0.5 - c)) / n ** 0.5)
    return c


print("Newton-Raphson method:")
newton_raphson_method(2202)
print("")
print("Bisection method:")
bisection_method(2202)
print("")
print("Regula falsi method")
regula_falsi_method(2202)
