
from helpers import *

# The grid
re = -1
rc = 1
n = 251                  # spatial grid points
x = np.linspace(re, rc, n)     # the grid
hgx = (rc-re)/(n-1)      # spatial grid spacing
hgt = 0.02               # time step

endtime = 100            # end time
tds = 2                 # downsampling in time
nsteps = int(endtime/tds/hgt)  # number of time steps

# Initial data
x0 = 0.5*(re+rc)
sigma = (re+rc)/40
initial = np.zeros((2, len(x)))
initial[1, :] = gaussian(x, x0, sigma)

# Coefficients
ell = 2  # the spherical harmonic mode
coefs = set_coefs(x, ell, re, rc)

# Solve the wave equation
data = wave1dsolve(initial, coefs, hgx, hgt, tds, nsteps)
