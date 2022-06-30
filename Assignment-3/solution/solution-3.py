from cProfile import label
from matplotlib.markers import MarkerStyle
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.linalg import lstsq
from pylab import *


def Q1():
    N = 101                           # no of data points
    k = 9                             # no of sets of data with varying noise

    # generate the data points and add noise
    t = linspace(0, 10, N)              # t vector
    y = 1.05*sp.jn(2, t)-0.105*t       # f(t) vector
    Y = meshgrid(y, ones(k), indexing='ij')[0]  # make k copies
    scl = logspace(-1, -3, k)           # noise stdev
    n = dot(randn(N, k), diag(scl))     # generate k vectors
    yy = Y+n                          # add noise to signal

    # shadow plot
    plot(t, yy)
    xlabel(r'$t$', size=20)
    ylabel(r'$f(t)+n$', size=20)
    title(r'Plot of the data to be fitted')
    grid(True)
    savetxt("fitting.dat",
            c_[t, yy])  # write out matrix to file
    show()


def Q2():
    t, v2, v3, v4, v5, v6, v7, v8, v9, v10 = np.loadtxt(
        "D:\\2nd year\\4th sem\\EE2703 APL\\Assignment-3\\solution\\fitting.dat", unpack=True)

    return t, [v2, v3, v4, v5, v6, v7, v8, v9, v10]


def Q3():
    t, fn = Q2()
    sigma = np.logspace(-1, -3, 9)
    print("print len(fn) is :", len(fn))
    for i in range(9):
        plt.plot(t, fn[i], label=r'$\sigma$=%.3f' % sigma[i])
    plt.plot(t, fn[-1], label="True Value", color="black")
    plt.xlabel("t")
    plt.ylabel("f(t)+n(t)")
    plt.legend(loc=1)
    plt.show()


def Q4():
    t, fn = Q2()

    def g(A, B, t):
        g = A*sp.jn(2, t)+B*t
        return g

    plt.xlabel("t", size=20)
    plt.ylabel("g", size=20)
    g = g(1.05, -0.105, t)
    plt.plot(t, g)
    plt.title("Plot of g(t,A,B) and t")
    plt.legend()
    plt.show()


def Q5():
    def get_std_dev(f, g):
        sq_sum = 0
        for i in range(len(f)):
            sq_sum += (f-g)**2
        return (sq_sum/101)**0.5
    t, fn = Q2()
    stdev = get_std_dev(fn[-1], fn[0])
    print("len(t)=", len(t[::5]))
    print("len(fn[0])=", len(fn[0][::5]))
    plt.errorbar(t, fn[0], stdev, fmt="ro", label="ErrorBar")  # error
    plt.plot(t, fn[-1], color="black", label="f(t)")
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.legend()
    plt.show()


def Q6():
    def populate_M_p(A, B, t):
        # p = np.zeros([2, 1], dtype=float)
        p = np.array([A, B])
        # p[0][0] = A
        # p[1][0] = B
        x = sp.jn(2, t)
        M = np.c_[x, t]
        g = A*sp.jn(2, t)+B*t
        return M, p, g

    t, fn = Q2()
    M, p, g = populate_M_p(1.05, -0.105, t)
    print("g is :")
    print(g)

    # To check
    g_check = np.dot(M, p)
    print("g-check is :")
    print(g_check)
    if (g == g_check).all():
        print("it matchess!!")
    else:
        print("No it does not match!!")


def Q7():
    t, fn = Q2()
    A = np.array([i/100 for i in range(201)])
    B = np.array([(-0.2+i/100) for i in range(21)])
    for a in A:
        for b in B:
            sq_sum = 0
            g = a*sp.jn(2, t)+b*t
            for i in range(len(t)):
                sq_sum += (fn[0][i]-g[i])**2
            e = sq_sum/101  # mean squared error
            # print("for [A=", a, ",B=", b, "] mean squared error is :", e)


def Q8():
    t, fn = Q2()
    A = np.array([(i/100) for i in range(201)])
    B = np.array([(-0.2+i/100) for i in range(21)])
    E = np.zeros([201, 21], dtype=float)
    for i in range(len(A)):
        for j in range(len(B)):
            a = A[i]
            b = B[j]
            sq_sum = 0
            g = a*sp.jn(2, t)+b*t
            for k in range(len(t)):
                sq_sum += (fn[0][k]-g[k])**2
            e = sq_sum/101  # mean squared error
            E[i][j] = e
            if i == 200 and j == 20:
                print("E[-1][-1] in loop is :", E[-1][-1])
            # print("for [A=", a, ",B=", b, "] mean squared error is :", e)
    fig, ax = plt.subplots(1, 1)
    # print("E is :")
    # print(E)
    print("A is:")
    print(A)
    print("B is :")
    print(B)
    ax.contour(B, A, E)

    # To plot the exact point
    y = [1.05]
    x = [-0.105]
    ax.plot(x, y, marker="o", markerfacecolor="red")
    ax.annotate("Exact Location", (-0.105, 1.05))
    ax.set_title(r'Contour Plot of $\epsilon_{ij}$')
    ax.set_xlabel('B')
    ax.set_ylabel('A')
    plt.show()


def Q9(i):
    def populate_M(t):
        x = sp.jn(2, t)
        M = np.c_[x, t]
        return M

    t, fn = Q2()
    A = np.array([(i/100) for i in range(201)])
    B = np.array([(-0.2+i/100) for i in range(21)])
    M = populate_M(t)
    g = fn[i]
    p, res, rnk, s = lstsq(M, g)
    # print("p matrix using lstsp is :")
    # print(p)
    return p


def Q10():
    t, fn = Q2()
    sigma = np.logspace(-1, -3, 9)
    error_in_estimating_A = np.zeros(len(fn))
    error_in_estimating_B = np.zeros(len(fn))
    exact_p = np.array([1.05, -0.105])
    for i in range(len(fn)):
        p = Q9(i)
        error_in_estimating_A[i] = np.abs(exact_p-p)[0]
        error_in_estimating_B[i] = np.abs(exact_p-p)[1]
    plt.plot(sigma, error_in_estimating_A, label="A",
             markerfacecolor="red", marker="D", linestyle="dashed", color="red")
    plt.plot(sigma, error_in_estimating_B, label="B",
             markerfacecolor="blue", marker="D", linestyle="dashed", color="blue")
    plt.ylabel("Ms error")
    plt.xlabel("Noise standard Deviation")
    plt.legend()
    plt.show()


def Q11():
    t, fn = Q2()
    sigma = np.logspace(-1, -3, 9)
    error_in_estimating_A = np.zeros(len(fn))
    error_in_estimating_B = np.zeros(len(fn))
    exact_p = np.array([1.05, -0.105])
    for i in range(len(fn)):
        p = Q9(i)
        error_in_estimating_A[i] = np.abs(exact_p-p)[0]
        error_in_estimating_B[i] = np.abs(exact_p-p)[1]
    loglog(sigma, error_in_estimating_A, linestyle='',
           marker='o', color='r', label='Aerr', markerfacecolor="r")
    loglog(sigma, error_in_estimating_B, linestyle='',
           marker='o', color='b', label='Berr', markerfacecolor="b")
    stem(sigma, error_in_estimating_A, '-ro')
    stem(sigma, (error_in_estimating_B), '-bo')
    plt.xlabel("Ms error")
    plt.ylabel("Noise standard Deviation")
    plt.legend()
    plt.show()


# type the function of each question here
Q1()
