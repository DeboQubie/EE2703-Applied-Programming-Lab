'''
    NAME:     Debojyoti Mazumdar
    ROLL NO.: EE20B030
    DATE:     12-05-2022
'''

# Importing modules
from pylab import *
import sympy as sp
import sys

# Initialising variables

# Independent variables
l = 0.5                       # quarter of a wavelength
c = 2.9979e8                  # speed of light
mu0 = (4*pi)*(10**-7)         # permeability of free space
N = 100                       # later make it N=100
Im = 1.0                      # current injected into the antenna
a = 0.01                      # radius of wire

# Taking input for the value of N if given
if len(sys.argv) != 1:
    N = int(sys.argv[1].lstrip())


# Dependent variables
wavelength = 4*l
f = c/wavelength              # frequency
k = 2*pi*(1/wavelength)       # wavenumber
dz = l/N                      # spacing of current samples


#--------- PART-1 ---------------------------------------------#
'''
    Pseudo-code:

    z = array from -N*dz to N*dz
    u = array z without -N*dz , 0 and N*dz
    I = array of length 2*N+1 with I[0]=I[-1]=0 and I[N]=Im
    J = array of length 2*N-1

'''


# Locations in the antenna
z = np.arange(-N, N+1, 1)*dz
# Unknown locations in the antenna
u = np.concatenate((z[1:N], z[N+1:-1]), axis=None)
# Current values
I = np.zeros((2*N)+1)
I[0] = 0
I[-1] = 0
I[N] = Im
# Current values in Unknown locations
J = np.zeros((2*N)-2)


#--------- PART-2 ---------------------------------------------#
'''
    Pseudo-code:

    function(input N):
        M = (2*N-2,2*N-2) identity matrix multiplied with (1/2*pi*a)
        return M
        
'''


# Function to compute M
def compute_M(N):
    M = np.identity((2*N)-2, dtype=float)*(1/(2*pi*a))
    return M


#--------- PART-3 ----------------------------------------------#
'''
    pseudo-code:
    
    Rz = sqrt(a^2 + (z-zj)^2)
    Ru = sqrt(a^2 + (u-uj)^2)
    R = Matrix with rows as Ru for all j
    RiN = sqrt(a^2 + u^2)

    P = Matrix obtained by the given expression with R
    Pb = Vector obtained by the given expression with Rb 

'''

# Finding the vectors Ru and Rz
j = sp.symbols("j")
Rz = (np.square(a)+np.square(z-j*dz))**0.5
Ru = (np.square(a)+np.square(u-j*dz))**0.5


# Initialising R ad RiN
x = np.linspace(-N+1, N-1, num=(2*N)-1, dtype=float)    # column index range
y = np.linspace(-N+1, N-1, num=(2*N)-1, dtype=float)    # row index range
# column indexes without 0
x = np.delete(x, N-1, 0)
y = np.delete(y, N-1, 0)                                # row indexes without 0
Y, X = np.meshgrid(y, x, sparse=False)
R = np.sqrt(np.square(a)+np.square((X-Y)*dz))
RiN = (a**2 + u**2)**0.5

# Matrix P calculated from matrix R
P = (mu0/(4*pi))*(np.exp(-1j*k*R))*(dz/R)
# Matrix PB calculated from matrix RiN
PB = (mu0/(4*pi))*(np.exp(-1j*k*RiN))*(dz/RiN)


#--------- PART-4 ---------------------------------------------#
'''
    pseudo-code:

    Q = Matrix obtained by the given expression using R and P
    QB = Vector obtained by the given expression using RiN and PB

'''

# Matrix Q calculated from P
Q = P*(a/mu0)*((1j*k/R)+(1/(R**2)))

# Matrix QB calculated from PB
QB = np.transpose(PB*(a/mu0)*((1j*k/RiN)+(1/(RiN**2))))


#--------- PART-5 ---------------------------------------------#
'''
    pseudo-code:

    M = function call to compute M
    J = inv(M-Q) x Qb x Im
    I[unknown current values indexes] = J
    I_approx = Im*sin(k*(l - abs(z)))

    Plot the calculated current and approximate current expression

'''

# Computing J and then followed by I
M = compute_M(N)
J = np.dot(np.linalg.inv(M-Q), QB)*Im

# Taking only the real part of the unknown currents and filling up I vector
I[1:N] = J[:N-1]
I[N+1:-1] = J[N-1:]

# Values of I according to the formula given
I_approx = Im*np.sin(k*(l - abs(z)))

# Plotting the I vector and I_approx vector
plot(z, I, "r", label="Actual current values")
plot(z, I_approx, "b", label="Calculated current values")
ylabel("Current(in A)")
xlabel("z (in meters)")
legend(loc="center")
title("Plot of calculated current and approximate current expression")
show()
