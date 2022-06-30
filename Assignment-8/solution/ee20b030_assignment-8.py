from pylab import *


# Magnitude and phase plot of Fourier Transform of sin(5t)
def Q1_example_1():
    t = linspace(0, 2*pi, 129)
    t = t[:-1]
    y = sin(5*t)
    Y = fftshift(fft(y))/128.0
    w = linspace(-64, 63, 128)
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), lw=2)
    xlim([-10, 10])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\sin(5t)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), "ro", lw=2)
    ii = where(abs(Y) > 1e-3)
    plot(w[ii], angle(Y[ii]), "go", lw=2)
    xlim([-10, 10])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    savefig("fig9-2.png")
    show()


# Magnitude and phase plot of Fourier Transform of (1+0.1cos(t))cos(10t)
def Q1_example_2():
    t = linspace(-4*pi, 4*pi, 513)
    t = t[:-1]
    y = (1+0.1*cos(t))*cos(10*t)
    Y = fftshift(fft(y))/512.0
    w = linspace(-64, 64, 513)
    w = w[:-1]
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), lw=2)
    xlim([-15, 15])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\left(1+0.1\cos\left(t\right)\right)\cos\left(10t\right)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), "ro", lw=2)
    xlim([-15, 15])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    savefig("fig9-4.png")
    show()


# Plotting the spectrum of sin(t)**3
def Q2_part_1():
    t = linspace(0, 2*pi, 129)
    t = t[:-1]
    y = sin(t)**3
    Y = fftshift(fft(y))/128.0
    w = linspace(-64, 63, 128)
    figure()
    plot(w, abs(Y), lw=2)
    xlim([-10, 10])
    ylabel(r"$|Y|$", size=16)
    xlabel(r"$\omega$", size=16)
    title(r"Spectrum of $\sin^3(t)$")
    grid(True)
    show()


# Plotting the spectrum of cos(t)**3
def Q2_part_2():
    t = linspace(0, 2*pi, 129)
    t = t[:-1]
    y = cos(t)**3
    Y = fftshift(fft(y))/128.0
    w = linspace(-64, 63, 128)
    figure()
    plot(w, abs(Y), lw=2)
    xlim([-10, 10])
    ylabel(r"$|Y|$", size=16)
    xlabel(r"$\omega$", size=16)
    title(r"Spectrum of $\cos^3(t)$")
    grid(True)
    show()


# Plotting the fourier transform of cos(20t+5cos(t))
def Q3():
    t = linspace(-4*pi, 4*pi, 513)
    t = t[:-1]
    y = cos(20*t+5*cos(t))
    Y = fftshift(fft(y))/512.0
    w = linspace(-64, 64, 513)
    w = w[:-1]
    figure()
    plot(w, abs(Y), lw=2)
    xlim([-40, 40])
    ylabel(r"$|Y|$", size=16)
    xlabel(r"$\omega$", size=16)
    title(r"Spectrum of $\cos(20t+5\cos(t))$")
    grid(True)
    show()
    ii = where(abs(Y) > 1e-3)
    plot(w[ii], angle(Y[ii]), "go", lw=2)
    xlim([-40, 40])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    show()


# plotting the Fourier Transform of gaussian exp(-0.5*(t^2))
def Q4():
    T = 2*pi
    N = 128
    tolerance = 1e-15
    error = tolerance+1

    while error > tolerance:
        t = linspace(-T/2, T/2, N+1)[:-1]
        w = N/T * linspace(-pi, pi, N+1)[:-1]
        y = exp(-0.5*t**2)

        Y = fftshift(fft(y))*T/(2*pi*N)
        Y_actual = (1/(sqrt(2*pi)))*exp(-0.5*w**2)
        error = mean(abs(abs(Y)-Y_actual))

        T = T*2
        N = N*2

    print("The Error in the fourier transform calculated is : ", error*100, " %")
    print("The frequency domain is : ", T/pi)
    figure()
    plot(w, abs(Y), lw=2)
    xlim([-10, 10])
    ylabel(r"$|Y|$", size=16)
    xlabel(r"$\omega$", size=16)
    title(r"Spectrum of $\exp(-t^2/2)$")
    grid(True)
    show()


Q4()
