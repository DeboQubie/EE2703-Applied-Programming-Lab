import math
from msilib.schema import Error
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import quad


def exp_x(x):
    return np.exp(x)


def cos_cos_x(x):
    return np.cos(np.cos(x))


# plotting the functions over the interval [-2pi,4pi]

x = np.linspace(-2*math.pi, 4*math.pi, 1500)
y_exp = exp_x(x)
y_cos_cos = cos_cos_x(x)

plt.title("Plot of exp(x) vs x")
plt.xlabel("x")
plt.ylabel("exp(x)")
plt.plot(x, y_exp, "r")
plt.show()

plt.title("Plot of cos(cos(x)) vs x")
plt.xlabel("x")
plt.ylabel("cos(cos(x))")
plt.plot(x, y_cos_cos, "b")
plt.show()

# plotting the functions with [0,2pi] interval repeating thrice

x_sub = np.linspace(0, 2*math.pi, 500)
y_exp_sub = exp_x(x_sub)
y_cos_cos_sub = cos_cos_x(x_sub)
y_exp_1 = np.tile(y_exp_sub, 3)
y_cos_cos_1 = np.tile(y_cos_cos_sub, 3)
x = np.linspace(-2*math.pi, 4*math.pi, 1500)

plt.title("plot of exp(x) vs x")
plt.xlabel("x")
plt.ylabel("exp(x)")
plt.plot(x, y_exp_1, "r1")
plt.plot(x, y_exp, "b")
plt.show()


plt.title("plot of cos(cos(x)) vs x")
plt.xlabel("x")
plt.ylabel("cos(cos(x))")
plt.plot(x, y_cos_cos_1, "r1")
plt.plot(x, y_cos_cos, "b")
plt.show()


def exp_x_cos_kx(x, k):
    return math.exp(x)*math.cos(k*x)


def exp_x_sin_kx(x, k):
    return math.exp(x)*math.sin(k*x)


def cos_cos_x_cos_kx(x, k):
    return math.cos(math.cos(x))*math.cos(k*x)


def cos_cos_x_sin_kx(x, k):
    return math.cos(math.cos(x))*math.sin(k*x)


def return_coeff_of_function(func):
    return_array = []
    if func == "exp(x)":
        for i in range(0, 26):
            a = quad(exp_x_cos_kx, 0, 2*math.pi, args=(i))[0]
            return_array += [a/math.pi]
        for i in range(1, 26):
            b = quad(exp_x_sin_kx, 0, 2*math.pi, args=(i))[0]
            return_array += [b/math.pi]
        return_array[0] = return_array[0]/2
    elif func == "cos(cos(x))":
        for i in range(0, 26):
            a = quad(cos_cos_x_cos_kx, 0, 2*math.pi, args=(i))[0]
            return_array += [a/math.pi]
        for i in range(1, 26):
            b = quad(cos_cos_x_sin_kx, 0, 2*math.pi, args=(i))[0]
            return_array += [b/math.pi]
        return_array[0] = return_array[0]/2
    else:
        return Error("Fourier series for this function is not defined in this python function")
    return np.array(return_array)


fs_coef_exp_x = return_coeff_of_function("exp(x)")
fs_coef_cos_cos_x = return_coeff_of_function("cos(cos(x))")
n = np.arange(0, 26, 1)


# plotting the coefficient values

# exp(x) coefficient plotting in semilogy plot
plt.title("plot of coefficients of cos(cos(x)) in semilogy")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.semilogy(n, np.abs(fs_coef_cos_cos_x[:26]), "bo", label="a")
plt.semilogy(n[1:], np.abs(fs_coef_cos_cos_x[26:]), "ro", label="b")
plt.legend(loc="upper right")
plt.show()

# cos(cos(x)) coefficient plotting in semilogy plot
plt.title("plot of coefficients of exp(x) in semilogy")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.semilogy(n, np.abs(fs_coef_exp_x[:26]), "bo", label="a")
plt.semilogy(n[1:], np.abs(fs_coef_exp_x[26:]), "ro", label="b")
plt.legend(loc="upper right")
plt.show()


# exp(x) coefficients plotting in loglog plot
plt.title("plot of coefficients of cos(cos(x)) in loglog")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.loglog(n, np.abs(fs_coef_cos_cos_x[:26]), "bo", label="a")
plt.loglog(n[1:], np.abs(fs_coef_cos_cos_x[26:]), "ro", label="b")
plt.legend(loc="upper right")
plt.show()

# cos(cos(x)) coefficients plotting in loglog plot
plt.title("plot of coefficients of exp(x) in loglog")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.loglog(n, np.abs(fs_coef_exp_x[:26]), "bo", label="a")
plt.loglog(n[1:], np.abs(fs_coef_exp_x[26:]), "ro", label="b")
plt.legend(loc="upper right")
plt.show()


# Estimating the function using least square method

def populate_A():
    x = np.linspace(0, 2*math.pi, 401)
    x = x[:-1]
    A = np.zeros((400, 51))  # allocate space for A
    A[:, 0] = 1  # col 1 is all ones
    for k in range(1, 26):
        A[:, 2*k-1] = np.cos(k*x)  # cos(kx) column
        A[:, 2*k] = np.sin(k*x)  # sin(kx) column
    return A


def converting_to_original_form(c):
    return_list = [c[0]]
    for i in range(1, len(c)):
        if i % 2 != 0:
            return_list += [c[i]]
    for i in range(1, len(c)):
        if i % 2 == 0:
            return_list += [c[i]]
    return np.array(return_list)


x = np.linspace(0, 2*math.pi, 401)
x = x[:-1]  # drop last term to have a proper periodic integral
b_exp_x = np.exp(x)  # exp(x) has been written to take a vector
# cos(cos(x)) has been written to take a vector
b_cos_cos_x = np.cos(np.cos(x))
A = populate_A()
c_exp_x = converting_to_original_form(np.linalg.lstsq(A, b_exp_x)[0])
c_cos_cos_x = converting_to_original_form(np.linalg.lstsq(A, b_cos_cos_x)[0])

# exp(x) both coefficient plotting in semilogy plot
plt.title("plot of coefficients of cos(cos(x)) in semilogy")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.semilogy(n, np.abs(fs_coef_cos_cos_x[:26]), "ro", label="a")
plt.semilogy(n[1:], np.abs(fs_coef_cos_cos_x[26:]), "ro", label="b")
plt.semilogy(n, np.abs(c_cos_cos_x[:26]), "go", label="a_new")
plt.semilogy(n[1:], np.abs(c_cos_cos_x[26:]), "go", label="b_new")
plt.legend(loc="upper right")
plt.show()

# cos(cos(x)) both coefficient plotting in semilogy plot
plt.title("plot of coefficients of exp(x) in semilogy")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.semilogy(n, np.abs(fs_coef_exp_x[:26]), "ro", label="a")
plt.semilogy(n[1:], np.abs(fs_coef_exp_x[26:]), "ro", label="b")
plt.semilogy(n, np.abs(c_exp_x[:26]), "go", label="a_new")
plt.semilogy(n[1:], np.abs(c_exp_x[26:]), "go", label="b_new")
plt.legend(loc="upper right")
plt.show()


# exp(x) coefficients plotting in loglog plot
plt.title("plot of coefficients of cos(cos(x)) in loglog")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.loglog(n, np.abs(fs_coef_cos_cos_x[:26]), "ro", label="a")
plt.loglog(n[1:], np.abs(fs_coef_cos_cos_x[26:]), "ro", label="b")
plt.loglog(n, np.abs(c_cos_cos_x[:26]), "go", label="a_new")
plt.loglog(n[1:], np.abs(c_cos_cos_x[26:]), "go", label="b_new")
plt.legend(loc="upper right")
plt.show()

# cos(cos(x)) coefficients plotting in loglog plot
plt.title("plot of coefficients of exp(x) in loglog")
plt.xlabel("n")
plt.ylabel("a and b")
plt.grid()
plt.loglog(n, np.abs(fs_coef_exp_x[:26]), "ro", label="a")
plt.loglog(n[1:], np.abs(fs_coef_exp_x[26:]), "ro", label="b")
plt.loglog(n, np.abs(c_exp_x[:26]), "go", label="a_new")
plt.loglog(n[1:], np.abs(c_exp_x[26:]), "go", label="b_new")
plt.legend(loc="upper right")
plt.show()

# To find the maximum absolute difference between the two sets of coefficients
print("The maximum error in a of exp(x) is : ",
      np.max(np.abs(c_exp_x[:26]-fs_coef_exp_x[:26])))
print("The maximum error in b of exp(x) is : ",
      np.max(np.abs(c_exp_x[26:]-fs_coef_exp_x[26:])))
print("The maximum error in a of cos(cos(x)) is : ", np.max(
    np.abs(c_cos_cos_x[:26]-fs_coef_cos_cos_x[:26])))
print("The maximum error in b of cos(cos(x)) is : ", np.max(
    np.abs(c_cos_cos_x[26:]-fs_coef_cos_cos_x[26:])))


# Plotting the new function generated by c and the function itself
def converting_to_required_format(fs_coef):
    return_list = np.zeros(51)
    return_list[0] = fs_coef[0]
    for i in range(1, 26):
        return_list[2*i-1] = fs_coef[i]
    for i in range(26, len(fs_coef)):
        return_list[2*(i-25)] = fs_coef[i]
    return return_list


x = np.linspace(0, 2*math.pi, 401)
x = x[:-1]
y_cos_cos_x = np.cos(np.cos(x))
y_exp_x = np.exp(x)
A = populate_A()
c_cos_cos_x = converting_to_required_format(fs_coef_cos_cos_x)
c_exp_x = converting_to_required_format(fs_coef_exp_x)
y_cos_cos_x_from_fs = np.matmul(A, c_cos_cos_x)
y_exp_x_from_fs = np.matmul(A, c_exp_x)


plt.title("Plot of exp(x)")
plt.grid()
plt.semilogy(x, y_exp_x_from_fs, "ro", label="exp(x) from fourier series")
plt.semilogy(x, y_exp_x, "blue", label="exp(x) function")
plt.show()

plt.title("Plot of cos(cos(x))")
plt.grid()
plt.plot(x, y_cos_cos_x_from_fs, "ro", label="cos(cos(x)) from fourier series")
plt.plot(x, y_cos_cos_x, "blue", label="cos(cos(x)) function")
plt.show()
