import numpy as np
import scipy.signal as sig
from pylab import *


def Q1(printing=False, filename="h.csv"):
    f = open(filename)

    b = f.readlines()
    for i in range(len(b)):
        b[i] = float(b[i].strip("\n"))
    f.close()

    b = np.array(b)

    if printing:
        print("The coefficients are :\n")
        print(b)

    return b


def Q2(plotting=False):

    # Obtaining and plotting the magnitude and the phase response of the digital filter
    w, h = sig.freqz(Q1())

    if plotting:
        figure()
        subplot(2, 1, 1)
        plot(w, abs(h), 'b', lw=2)
        ylabel(r"magnitude", size=16)
        title(r"Spectrum of Low Pass Filter")
        grid(True)
        subplot(2, 1, 2)
        plot(w, angle(h), 'r', lw=2)
        ylabel(r"Phase", size=16)
        xlabel(r"$\omega$", size=16)
        grid(True)
        show()

    return w, h


def Q3(plotting=False):
    n = np.arange(1, (2**10)+1)
    x = np.cos(0.2*pi*n) + np.cos(0.85*pi*n)

    if plotting:
        title("Input signal")
        plot(n, x, "b")
        xlabel("n")
        ylabel("Amplitude")
        show()

    return x, n


# linear convolution
def Q4(plotting=False):

    b = Q1()
    x, n = Q3()

    y = np.convolve(x, b, mode="same")

    if plotting:
        title("Output after linear convolution")
        plot(n, real(y), "b")
        xlabel("n")
        ylabel("Amplitude")
        show()

    return y


# circular convolution
def Q5(plotting=False):

    x, n = Q3()
    w, h = Q2()
    b = Q1()

    y = ifft(fft(x)*fft(concatenate((b, zeros(len(x)-len(b))))))

    if plotting:
        title("Output after circular convolution")
        plot(n, real(y), "b")
        xlabel("n")
        ylabel("Amplitude")
        show()

    return y


# circular convolution using linear convolution
def Q6(plotting=False):

    x, n = Q3()
    b = Q1()

    P = len(b)
    m = int(ceil(log2(P)))

    b_padded = np.concatenate((b, zeros((2**m)-P)))

    len_of_zero_array = (int(ceil(len(x)/2**m)))*(int(2**m))-len(x)
    x_padded = np.concatenate((x, np.zeros(len_of_zero_array)))

    y = []

    for i in range(int(len(x_padded)/(2**m))):
        x_i = np.concatenate((x_padded[i*(2**m):(i+1)*(2**m)], np.zeros(P-1)))
        x_i = x_padded[i*(2**m):(i+1)*(2**m)]
        y_i = ifft(
            fft(x_i)*fft(b_padded))
        y = np.concatenate((y, y_i))

    if plotting:
        title("circular convolution using linear convolution")
        plot(n, real(y), "b")
        xlabel("n")
        ylabel("Amplitude")
        show()

    return y


def file_reading_for_Q7(filename="x1.csv", printing=False):
    f = open(filename)

    b = f.readlines()
    for i in range(len(b)):
        b[i] = complex(b[i].strip(" ").strip("i\n")+"j")
    f.close()

    b = np.array(b)

    if printing:
        print("The coefficients are :\n")
        print(b)

    return b


# Circular correlation
def Q7():
    x = file_reading_for_Q7()

    x2 = np.roll(x, 5)

    cor = np.fft.ifftshift(np.correlate(x2, x, 'full'))
    n = linspace(0, len(cor)-1, len(cor))

    title("Correlation")
    plot(n, abs(cor), "b")
    xlabel("t")
    ylabel("correlation")
    xlim(0, 20)
    show()


Q6(plotting=True)
# Q7()
