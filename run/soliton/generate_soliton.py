# Generates solitonic profiles for the SPE in d=1,3
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt


def f1d(r,y,a):
    return np.vstack((y[1,:],
                      2 * a *y[2,:] * y[0,:],
                      y[3,:],
                      y[0,:]**2))

def f3d(r,y,a):
    return np.vstack((y[1,:]*r,
                      2 * a*(y[2,:] * y[0,:] - y[1,:]),
                      y[3,:]*r,
                      4*np.pi*y[0,:]**2 - 2*y[3,:]))

def bc(ya,yb):
    return np.array([ya[0]-1, ya[1], ya[3], yb[0]])

def jac1d(r,y,a):
    J = np.zeros((y.shape[0],y.shape[0],r.shape[0]))
    J[0,1,:] = 1
    J[1,0,:] = 2 * a * y[2,:]
    J[1,1,:] = 0
    J[1,2,:] = 2 * a * y[0,:]
    J[2,3,:] = 1
    J[3,0,:] = y[0,:]
    J[3,3,:] = 0
    return J

def jac3d(r,y,a):
    J = np.zeros((y.shape[0],y.shape[0],r.shape[0]))
    J[0,1,:] = r
    J[1,0,:] = 2 * a * y[2,:]
    J[1,1,:] = -2
    J[1,2,:] = 2 * a * y[0,:]
    J[2,3,:] = r
    J[3,0,:] = 8 * np.pi * y[0,:]
    J[3,3,:] = -2
    return J

a = 1
L = 10
r = np.linspace(0,L,256)
y = np.zeros((4,256));

# Initial Guess is a gaussian
sigma = 1
y[0] = np.exp(-0.5 * r**2/sigma**2)
y[1] = -r/sigma**2 * np.exp(-0.5 * r**2/sigma**2)
y[2] = -np.exp(-0.5 * r**2/sigma**2 ) - 1
y[3] = r/sigma**2 * np.exp(-0.5 * r**2/sigma**2)

res = integrate.solve_bvp(lambda r,y: f1d(r,y,a), 
                          bc, 
                          r, y, 
                          fun_jac = lambda r,y: jac1d(r,y,a), 
                          verbose=2, 
                          tol=1e-10,
                          max_nodes = 1000000)

# x = np.linspace(-10*L,10*L,1000)
# plt.plot(x, res.sol(np.abs(x))[2])
# plt.plot(x, np.abs(x)*integrate.quad(lambda r: res.sol(r)[0]**2, 0, L)[0])
# plt.grid()
# plt.show()

alpha = 1.58
v = 0
N = 2**14
L_sim = 200

x = np.linspace(-L_sim, L_sim, N) 
x0 = 10
v = 0 
dphi = 0

f1 = alpha * res.sol(np.sqrt(alpha) * np.abs(x - x0))[0]
f2 = alpha * res.sol(np.sqrt(alpha) * np.abs(x + x0))[0]

f1 = np.clip(f1, 0, None)
f2 = np.clip(f2, 0, None)

phase1 = np.exp(1.0j*v*x)
phase2 = np.exp(-1.0j*(v*x - dphi))

# psi_x = f1 * phase1 
psi_x = f1 * phase1 + f2 * phase2

fig,ax = plt.subplots()
plt.loglog(np.linspace(0, L_sim, N//2), alpha * 
         res.sol(np.sqrt(alpha) * np.linspace(0,L_sim,N//2))[0])
plt.show()

np.savetxt("psi_constructive_%d_%d_%f_%d.txt"%(L_sim, alpha, a, N), psi_x, fmt= "%.15e %.15e")

