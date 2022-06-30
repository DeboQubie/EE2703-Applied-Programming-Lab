import math
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import sys

if len(sys.argv) == 5:
    Nx = sys.argv[1]  # size along x
    Ny = sys.argv[2]  # size along y
    radius = sys.argv[3]  # radius of central lead
    Niter = sys.argv[4]  # number of iterations to perform

else:
    Nx = 30  # size along x
    Ny = 30  # size along y
    radius = 8  # radius of central lead
    Niter = 1500  # number of iterations to perform

# initialize potential
phi = np.zeros((Nx, Ny), dtype=float)
x, y = np.linspace(-0.5, 0.5, num=Nx,
                   dtype=float), np.linspace(-0.5, 0.5, num=Ny, dtype=float)
Y, X = np.meshgrid(y, x, sparse=False)
ii = np.where(X**2+Y**2 < (0.35)**2)
phi[ii] = 1.0

# plot potential
plt.xlabel("X")
plt.ylabel("Y")
plt.contourf(X, Y, phi)
plt.colorbar()
plt.show()


def update_phi(phi):
    phi[1:-1, 1:-1] = 0.25 * \
        (phi[1:-1, 0:-2]+phi[1:-1, 2:]+phi[0:-2, 1:-1]+phi[2:, 1:-1])
    return phi


def assert_boundaries(phi):
    phi[1:-1, 0] = phi[1:-1, 1]  # left side normal component 0
    phi[1:-1, -1] = phi[1:-1, -2]  # right side normal component 0
    phi[0, 1:-1] = phi[1, 1:-1]  # top side normal component 0
    return phi


# performing the iteration and saving the error
errors = np.zeros(Niter)
for k in range(Niter):
    oldphi = phi.copy()
    phi = update_phi(phi)
    phi = assert_boundaries(phi)
    phi[ii] = 1.0  # To make the potential at the points of wire 1
    errors[k] = (abs(phi-oldphi)).max()

# Keeping the x value of error plots
x = np.arange(1, Niter+1)

# log(error) vs iter plot
plt.xlabel("Iteration number")
plt.ylabel("error value in log scale")
plt.semilogy(x, errors, "ro")
plt.show()

# fitting a line in the error graph


def populate_M():
    x = np.ones(Niter).T
    y = np.arange(1, Niter+1)
    y = y.T
    return c_[x, y]


def get_max_error(A, B, Niter):
    return -(A/B)*(math.e**(B*(Niter+0.5)))


M = populate_M()
logA, B = lstsq(M, np.log(errors))[0]
A = math.e**logA
fitting_error = A*np.exp(B*x)


# plotting the fitting points and error line in semilogy plot
plt.xlabel("Iteration number")
plt.ylabel("error value in log scale")
plt.semilogy(x, errors, "red")
plt.semilogy(x[::50], fitting_error[::50], "bo")
plt.show()

# To get the max error
print("Max error is :")
print(get_max_error(A, B, Niter))

# Surface plot of final potential
fig4 = figure(4)  # open a new figure
ax = p3.Axes3D(fig4)  # Axes3D is the means to do a surface plot
title("The 3-D surface plot of the potential")
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1, cmap=cm.jet)
plt.show()

# plotting 2d contour of final potential
plt.title("2D Contour plot of potential")
plt.xlabel("X")
plt.ylabel("Y")
x_c, y_c = ii
plt.plot((x_c-Nx/2)/Nx, (y_c-Ny/2)/Ny, 'ro')
plt.contourf(Y, X[::-1], phi)
plt.colorbar()
plt.show()

# To draw the vector plot of the currents

Jx, Jy = (1/2*(phi[1:-1, 0:-2]-phi[1:-1, 2:]),
          1/2*(phi[:-2, 1:-1]-phi[2:, 1:-1]))
plt.title("Vector plot of current flow")
plt.quiver(Y[1:-1, 1:-1], -X[1:-1, 1:-1], -Jx[:, ::-1], -Jy)
x_c, y_c = np.where(X**2+Y**2 < (0.35)**2)
plt.plot((x_c-Nx/2)/Nx, (y_c-Ny/2)/Ny, 'ro')
plt.show()
