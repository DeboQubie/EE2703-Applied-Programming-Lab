##############################################################################
##############                  End Semester Exam               ##############
##############              Aniruddha S     EE20B010            ##############
##############################################################################

# Import commands

from pylab import *
import numpy as np
import sympy as sp

# Initialisation of constants and dependant variables

l = 0.5                    # Half length of dipole antenna
c = 2.9979e8               # Speed of light in vacuum
mu0 = 4e-7*np.pi           # Permeability of free space
N = 4                    # Number of sections in each half
Im = 1                     # Current injected into the antenna
a = 0.01                   # Radius of wire

wl = 4*l                   # Wavelength
freq = c/l                 # Frequency
k = 2*np.pi/wl             # Wave number
dz = l/N                   # Length of each piece of the wire


### Part 1 ###

# Pseudo-code

# z = array from -N to N
# u = remove -N, 0, N from z
# I = array of 2N+1 zeros
# Set I[0], I[N], I[2N] to 0, Im, 0
# J = remove 0th, Nth and 2Nth element of I

# Code

z = np.linspace(-N, N, 2*N+1)*dz          # Locations in the antenna
u = np.concatenate((z[1:N], z[N+1:-1]))   # Unknown locations
I = np.zeros(2*N+1)                       # Vector I
I[0], I[N], I[2*N] = 0, Im, 0             # Setting known currents
J = np.concatenate((I[1:N], I[N+1:-1]))   # Vector J

# Printing obtained vector arrays
print('Vector z: ', z)
print('Vector u: ', u)
print('Vector I: ', I)
print('Vector J: ', J)
print('')


### Part 2 ###

# Pseudo-code

# Create function M(N)
# Return 1/(2*pi*a) times identity matrix of order 2N-2

# Code

# Function to find matrix M
def M(N):
    return 1/(2*np.pi*a)*np.identity(2*N-2)


### Part 3 ###

# Pseudo-code

# Rz = sqrt(a^2 + (z-zj)^2)
# Ru = sqrt(a^2 + (u-uj)^2)
# Rb = sqrt(a^2 + u^2)
# R = Matrix with rows as Ru for all j
# P = Matrix obtained by the given expression with R
# Pb = Vector obtained by the given expression with Rb

# Code

j = sp.symbols('j')               # Sympy symbol j for denoting index
Rz = (a**2 + (z-j*dz)**2)**0.5    # Rz as a function of j
Ru = (a**2 + (u-j*dz)**2)**0.5    # Ru as a function of j
Rb = (a**2 + u**2)**0.5           # R vector at middle position

# Calculation of Matrix R containing all observer-source distances
R = (np.ones((2*N-2, 1))*u)
R = R - np.transpose(np.ones((2*N-2, 1))*u)
R = (R**2 + a**2)**0.5

P = (mu0/(4*np.pi))*np.exp(-1j*k*R)*dz/R      # Matrix P from Matrix R
Pb = (mu0/(4*np.pi))*np.exp(-1j*k*Rb)*dz/Rb   # Array Pb from Array Rb

# Printing obtained vector arrays and matrices
print('Vector Rz: ', Rz)
print('Vector Ru: ', Ru)
# print('Vector Rb:\n\n\n\n ', Rb)
# print('Matrix R:\n\n ', R)
# print('Matrix P:\n\n\n ', P)
# print('Vector Pb: \n\n\n', Pb)
print('')

### Part 4 ###

# Pseudo-code

# Q = Matrix obtained by the given expression with R
# Qb = Vector obtained by the given expression with Rb

# Code

Q = P*(a/mu0)*(1/R**2 + 1j*k/R)          # Matrix Q from Matrix R
Qb = Pb*(a/mu0)*(1/Rb**2 + 1j*k/Rb)      # Array Qb from Array Rb

# Printing obtained vector arrays and matrices
print('Matrix Q:\n\n\n\n ', Q)
print('Vector Qb:\n\n\n\n ', Qb)
print('')

### Part 5 ###

# Pseudo-code

# J = inv(M-Q) x Qb x Im
# I = append zeros at the ends of J and Im in the middle
# Plot the calculated current and approximate current expression

# Code

# Calculation of J and hence I from
J = np.dot(np.linalg.inv(M(N)-Q), Qb*Im).real
I = np.array([0] + J[:N-1].tolist() + [Im] + J[N-1:].tolist() + [0])

# Printing obtained vector arrays
# print('Vector J: ', J)
# print('Vector I: ', I)
print('')

# Plot of the calculated current and approximate current expression
plot(I)
plot(Im*np.sin(k*(l - abs(z))))    # Approximate current expression
title('Plot of calculated current and approximate current expression')
xlabel('Element index')
ylabel('Current')
legend(['Calculated current', 'Approximate current expression'])
show()
