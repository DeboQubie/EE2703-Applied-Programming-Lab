import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt

# generates transfer function for forced damped spring system


def transfer_func_generator(decay, freq):
    p = np.polymul([1, 0, 2.25], [1, 2*decay, decay*decay+freq*freq])
    return sp.lti([1, decay], p)


# question-1
t = np.linspace(0, 50, 10001)
x_low_decay = sp.impulse(transfer_func_generator(0.5, 1.5), None, t)[1]
plt.title("plot of spring funtion with decay=0.5/s vs t")
plt.plot(t, x_low_decay, "r")
plt.xlabel("time")
plt.ylabel("x(t)")
plt.show()

# question-2
x_high_decay = sp.impulse(transfer_func_generator(0.05, 1.5), None, t)[1]
plt.title("plot of spring funtion with decay=0.05/s vs t")
plt.plot(t, x_high_decay, "r")
plt.xlabel("time")
plt.ylabel("x(t)")
plt.show()

# question-3
H = sp.lti(1, [1, 0, 2.25])
freq = np.linspace(1.4, 1.6, 5)
t = np.linspace(0, 200, 100001)
for f in freq:
    f_t = np.cos(f*t)*np.exp(-0.05*t)
    x = sp.lsim(H, f_t, t)[1]
    plt.title("plot of x(t) vs t for freq="+str(f))
    plt.xlabel("t")
    plt.ylabel("x(t)")
    plt.plot(t, x, "r")
    plt.show()


# question-4
X = sp.lti([1, 0, 2], [1, 0, 3, 0])
Y = sp.lti([2], [1, 0, 3, 0])

t = np.linspace(0, 20, 5001)
x = sp.impulse(X, None, t)[1]
y = sp.impulse(Y, None, t)[1]

plt.title("x(t) and y(t) vs time")
plt.plot(t, x, "r", label="x(t)")
plt.plot(t, y, "b", label="y(t)")
plt.xlabel("t")
plt.ylabel("f(t)")
plt.legend(loc="lower right")
plt.show()


# question-5
H = sp.lti([10**12], [1, 10**8, 10**12])

# plotting the bode plot
w, mag, phase = sp.bode(H)
plt.figure()
plt.title("Magnitude response of H(s)")
plt.semilogx(w, mag)    # Bode magnitude plot
plt.figure()
plt.title("Phase response of H(s)")
plt.semilogx(w, phase)  # Bode phase plot
plt.show()


# question-6
t = np.linspace(0, 30e-6, 10000)
vi_t = np.cos((10**3)*t)-np.cos((10**6)*t)
vo_t = sp.lsim(H, vi_t, t)[1]

# plotting
plt.title("vo(t) vs time for t<30\u03BCsec")
plt.xlabel("t")
plt.ylabel("vo(t)")
plt.plot(t, vo_t, "r", label="vo(t)")
plt.show()

t = np.linspace(0, 30e-3, 10000)
vo_t = sp.lsim(H, vi_t, t)[1]

# plotting
plt.title("vo(t) vs time for t<30msec")
plt.xlabel("t")
plt.ylabel("vo(t)")
plt.plot(t, vo_t, "r", label="vo(t)")
plt.show()
