from pylab import *


# with y = sin(sqrt(2)*t)
def Q1_Example_1():

    t = linspace(-4*pi, 4*pi, 257)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt
    n = arange(256)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/256))
    y = sin(sqrt(2)*t)
    y = y*wnd
    y[0] = 0  # the sample corresponding to -tmax should be set zeroo
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/256.0
    w = linspace(-pi*fmax, pi*fmax, 257)
    w = w[:-1]
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), 'b', w, abs(Y), 'bo', lw=2)
    xlim([-4, 4])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\sin\left(\sqrt{2}t\right)\times w(t)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), 'ro', lw=2)
    xlim([-4, 4])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    show()


# with y = sin(1.25*t)
def Q1_Example_2():

    t = linspace(-4*pi, 4*pi, 257)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt
    n = arange(256)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/256))
    y = sin(1.25*t)
    y = y*wnd
    y[0] = 0  # the sample corresponding to -tmax should be set zeroo
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/256.0
    w = linspace(-pi*fmax, pi*fmax, 257)
    w = w[:-1]
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), 'b', w, abs(Y), 'bo', lw=2)
    xlim([-4, 4])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\sin\left(1.25t\right)\times w(t)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), 'ro', lw=2)
    xlim([-4, 4])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    show()


def Q2_Without_Hamming_Window():

    t = linspace(-pi, pi, 65)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt
    omega0 = 0.86
    y = (cos(omega0*t))**3
    y[0] = 0  # the sample corresponding to -tmax should be set zeroo
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/64.0
    w = linspace(-pi*fmax, pi*fmax, 65)
    w = w[:-1]
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), lw=2)
    xlim([-10, 10])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\cos^{3}\left(\omega_{0}t\right)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), 'ro', lw=2)
    xlim([-10, 10])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    savefig("fig10-1.png")
    show()


def Q2_With_Hamming_Window():
    t = linspace(-pi, pi, 65)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt
    n = arange(64)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/63))
    omega0 = 0.86
    y = ((cos(omega0*t))**3)*wnd
    y[0] = 0  # the sample corresponding to -tmax should be set zeroo
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/64.0
    w = linspace(-pi*fmax, pi*fmax, 65)
    w = w[:-1]
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), lw=2)
    xlim([-8, 8])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\cos^{3}\left(\omega_{0}t\right)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), 'ro', lw=2)
    xlim([-8, 8])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    show()


# generates cos(omega0t + delta)
def gen_cos(omega0, delta):
    t = linspace(-pi, pi, 129)
    t = t[:-1]
    y = cos((omega0*t) + delta)
    y[0] = 0

    return y


def omega0_and_delta_estimator(y, plotting=False):
    t = linspace(-pi, pi, 129)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt

    # generating the window function
    n = arange(128)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/127))
    y = y*wnd

    # generating Y
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/128.0
    w = linspace(-pi*fmax, pi*fmax, 129)
    w = w[:-1]

    # Plotting spectrum of y
    if plotting == True:
        figure()
        subplot(2, 1, 1)
        plot(w, abs(Y), lw=2)
        xlim([-8, 8])
        ylabel(r"$|Y|$", size=16)
        title(r"Spectrum of $\cos\left(\omega_{0}t + \delta\right)$")
        grid(True)
        subplot(2, 1, 2)
        plot(w, angle(Y), 'ro', lw=2)
        xlim([-8, 8])
        ylabel(r"Phase of $Y$", size=16)
        xlabel(r"$\omega$", size=16)
        grid(True)
        show()

    # Estimating omega
    ii = where(w > 0)
    omega = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))  # weighted average
    print("omega = ", omega)

    # Estimating delta
    ii_1 = np.where(np.logical_and(np.abs(Y) > 1e-4, w > 0))[0]
    np.sort(ii_1)
    points = ii_1[1:3]
    # weighted average for first 2 points
    print("delta = ", np.sum(np.angle(Y[points]))/len(points))


def Q3(omega0, delta, plotting=False):
    y = gen_cos(omega0, delta)

    print("omega and delta value without noise :\n")
    omega0_and_delta_estimator(y, plotting)


# generates cos(omega0t + delta) + noise
def gen_cos_with_noise(omega0, delta):
    t = linspace(-pi, pi, 129)
    t = t[:-1]
    y = cos((omega0*t) + delta)
    y[0] = 0

    # generating white gaussian noise
    n = 0.1*randn(128)

    return y+n


def Q4(omega0, delta, plotting=False):
    y = gen_cos_with_noise(omega0, delta)

    print("omega and delta value with noise is :\n")
    omega0_and_delta_estimator(y, plotting)


# generates the "chirped" signal
def gen_chirp_signal(N):
    t = linspace(-pi, pi, N+1)
    t = t[:-1]
    y = cos(16*t*(1.5+(t/(2*pi))))
    y[0] = 0

    return y


def plotting_spectrum_of_chirp_signal(y, N):

    t = linspace(-pi, pi, N+1)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt

    # generating the window function
    n = arange(N)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/(N-1)))
    y = y*wnd

    # generating Y
    y = fftshift(y)  # make y start with y(t=0)
    Y = fftshift(fft(y))/N
    w = linspace(-pi*fmax, pi*fmax, N+1)
    w = w[:-1]

    # Plotting spectrum of y
    figure()
    subplot(2, 1, 1)
    plot(w, abs(Y), lw=2)
    xlim([-60, 60])
    ylabel(r"$|Y|$", size=16)
    title(r"Spectrum of $\cos\left(16t\left(1.5+\frac{t}{2\pi}\right)\right)$")
    grid(True)
    subplot(2, 1, 2)
    plot(w, angle(Y), 'ro', lw=2)
    xlim([-8, 8])
    ylabel(r"Phase of $Y$", size=16)
    xlabel(r"$\omega$", size=16)
    grid(True)
    show()


def Q5():
    N = 1024
    plotting_spectrum_of_chirp_signal(gen_chirp_signal(N), N)


def gen_Y_from_t(t):
    N = 64
    n = arange(N)
    wnd = fftshift(0.54+0.46*cos(2*pi*n/(N-1)))
    y = cos(16*t*(1.5 + t/(2*pi)))*wnd
    y[0] = 0
    y = fftshift(y)
    Y = fftshift(fft(y))/N

    return Y


def Q6():

    N = 1024
    t = linspace(-pi, pi, N+1)
    t = t[:-1]
    dt = t[1]-t[0]
    fmax = 1/dt

    # splitting t into 64 samples wide
    N_array = 64
    t_array = split(t, (N/N_array))

    # finding Y for each t_array
    Y = zeros((16, 64), dtype=complex)
    for i in range(len(t_array)):
        Y[i] = gen_Y_from_t(t_array[i])

    # plotting "time-frequency" plot of Y
    t = t[::N_array]
    w = linspace(-pi*fmax, pi*fmax, N_array+1)
    w = w[:-1]
    t, w = meshgrid(t, w)

    fig1 = figure()
    ax = fig1.add_subplot(111, projection='3d')
    surf = ax.plot_surface(w, t, angle(Y).T, cmap='viridis',
                           linewidth=0, antialiased=False)
    fig1.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_title('Phase plot')
    ylabel(r"$\omega\rightarrow$")
    xlabel(r"$t\rightarrow$")
    show()


Q6()
