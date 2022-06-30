'''
    NAME:Debojyoti Mazumdar
    ROLL NO.:EE20B030
    ASSIGNMENT:7
'''

# importing the required modules
import sympy as sm
import pylab as p
import scipy.signal as sp
import numpy as np


# returns the unit step function
def u(t):
    return np.where(t > 0, 1, 0)


# returns the sinusoidal function
def vi(t):
    return (np.sin(2e3*np.pi*t) + np.cos(2e6*np.pi*t))*(t > 0)


# returns the damped sinusoidal function
def vi_ds(t, w: float = 1e6, a: float = 500):
    return np.where(t < 0, 0, np.sin(w*t)*np.exp(-1*a*t))


# returns the Voltage solution of the given lowpass filter
def lowpass(R1, R2, C1, C2, G, Vi):
    s = sm.symbols("s")
    A = sm.Matrix([[0, 0, 1, -1/G], [-1/(1+s*R2*C2), 1, 0, 0],
                   [0, -G, G, 1], [-(1/R1)-(1/R2)-s*C1, 1/R2, 0, s*C1]])
    b = sm.Matrix([0, 0, 0, -Vi/R1])
    V = A.inv()*b
    return (A, b, V)


# returns the Voltage solution of the given highpass filter
def highpass(R1, R3, C1, C2, G, Vi):
    s = sm.symbols("s")
    A = sm.Matrix([[0, -1, 0, 1/G],
                   [s*C2*R3/(s*C2*R3+1), 0, -1, 0],
                   [0, G, -G, 0],
                   [-s*C2-(1/R1)-s*C1, 0, s*C2, 1/R1]])
    b = sm.Matrix([0, 0, 0, -Vi*s*C1])
    V = A.inv()*b
    return (A, b, V)


# converts the sympy expression to an lti response
def symExpToLTI(exp, s=sm.symbols('s')):
    exp_num, exp_den = exp.as_numer_denom()
    num_coeffs = sm.Poly(exp_num, s).all_coeffs()
    den_coeffs = sm.Poly(exp_den, s).all_coeffs()
    num = p.poly1d([float(n) for n in num_coeffs])
    den = p.poly1d([float(n) for n in den_coeffs])
    return sp.lti(num, den)


# returns the transfer function of the given lowpass filter
def transfer_function_of_LPF(Vi=1, R1: float = 10e3, R2: float = 10e3, C1: float = 1e-9, C2: float = 1e-9, G: float = 1.586):
    A, b, V = lowpass(R1, R2, C1, C2, G, Vi)
    V = A.inv()*b
    Vo = V[3]
    print(Vo)
    ww = p.logspace(0, 8, 801)
    ss = 1j*ww
    s = sm.symbols("s")
    hf = sm.lambdify(s, Vo, "numpy")
    v = hf(ss)
    p.loglog(ww, abs(v), lw=2)
    p.title(r'$\mid H(j{\omega}\mid$ of given LPF')
    p.xlabel(r'$\omega (log) \longrightarrow$')
    p.ylabel(r'${\mid H(j{\omega}\mid} (dB) \longrightarrow$')
    p.grid(True)
    p.show()
    return Vo


# returns the transfer function of the given highpass filter
def transfer_function_of_HPF(Vi=1, R1: float = 10e3, R3: float = 10e3, C1: float = 1e-9, C2: float = 1e-9, G: float = 1.586):
    A, b, V = highpass(R1, R3, C1, C2, G, Vi)
    V = A.inv()*b
    Vo = V[3]
    print(Vo)
    ww = p.logspace(0, 8, 801)
    ss = 1j*ww
    s = sm.symbols("s")
    hf = sm.lambdify(s, Vo, "numpy")
    v = hf(ss)
    p.loglog(ww, abs(v), lw=2)
    p.title(r'$\mid H(j{\omega}\mid$ of given HPF')
    p.xlabel(r'$\omega (log) \longrightarrow$')
    p.ylabel(r'${\mid H(j{\omega}\mid} (dB) \longrightarrow$')
    p.grid(True)
    p.show()
    return Vo


def Q1(t=np.linspace(0, 1e-3, 1001)):
    s = sm.symbols("s")
    Vo = transfer_function_of_LPF()
    H = symExpToLTI(Vo)
    u_ = u(t)
    _, u_resp, _ = sp.lsim(H, u_, t)
    p.plot(t, u_resp)
    p.title('Step response of given LPF')
    p.xlabel(r'time $\longrightarrow$')
    p.ylabel("output")
    p.show()


def Q2(t=np.linspace(0, 1e-3, 10001)):
    s = sm.symbols("s")
    Vo = transfer_function_of_LPF()
    H = symExpToLTI(Vo)
    _, vi_resp, _ = sp.lsim(H, vi(t), t)
    p.plot(t, vi_resp)
    p.title('Input response of the given LPF')
    p.xlabel(r'time $\longrightarrow$')
    p.ylabel("output")
    p.show()


def Q3():
    Vo = transfer_function_of_HPF()


def Q4(t=np.linspace(0, 1e-2, 100001)):
    s = sm.symbols("s")

    # output from HPF
    Vo = transfer_function_of_HPF()
    H = symExpToLTI(Vo)
    _, vi_ds_resp, _ = sp.lsim(H, vi_ds(t), t)
    p.plot(t, vi_ds(t))
    p.title('damped sinosoid')
    p.show()
    p.plot(t, vi_ds_resp)
    p.title('damped sinosoid response of the given HPF')
    p.xlabel(r'time $\longrightarrow$')
    p.ylabel("output")
    p.show()

    # output from LPF
    Vo = transfer_function_of_LPF()
    H = symExpToLTI(Vo)
    _, vi_ds_resp, _ = sp.lsim(H, vi_ds(t), t)
    p.plot(t, vi_ds(t))
    p.title('damped sinosoid')
    p.xlabel(r'time $\longrightarrow$')
    p.show()
    p.plot(t, vi_ds_resp)
    p.title('damped sinosoid response of the given LPF')
    p.xlabel(r'time $\longrightarrow$')
    p.ylabel("output")
    p.show()


def Q5(t=np.linspace(0, 1e-3, 1001)):
    s = sm.symbols("s")
    Vo = transfer_function_of_HPF()
    H = symExpToLTI(Vo)
    _, u_resp, _ = sp.lsim(H, u(t), t)
    p.plot(t, u_resp)
    p.title('Step response of the given HPF')
    p.xlabel(r'time $\longrightarrow$')
    p.ylabel("output")
    p.show()


Q1()
