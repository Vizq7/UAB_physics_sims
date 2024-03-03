import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle

def psi0(x, y, x0, y0, sigma=0.5, k=15*np.pi):
    return np.exp(-1/2*((x-x0)**2 + (y-y0)**2)/sigma**2)*np.exp(-1j*k*(x-x0))


def double_walls(psi, j0, j1, i0, i1, i2, i3):
    psi = np.asarray(psi)
    psi[0:i3, j0:j1] = 0
    psi[i2:i1,j0:j1] = 0
    psi[i0:,  j0:j1] = 0
    
    return psi
    
def one_walls(psi, j0, j1, i0, i1, i2, i3):
    psi = np.asarray(psi)
    psi[0:i2, j0:j1] = 0
    psi[i1:,j0:j1] = 0
    
    return psi

L = 8
dy = 0.05
dt = dy**2/4
Nx = int(L/dy) + 1
Ny = int(L/dy) + 1
Nt = 200
ct1 = -dt/(2j*dy**2)
ct2 = -dt/(2j*dy**2)

x0 = L/5
y0 = L/2
w = 0.1
s = 0.8
a = 0.4

j0 = int(1/(2*dy)*(L-w))
j1 = int(1/(2*dy)*(L+w))
i0 = int(1/(2*dy)*(L+s) + a/dy) 
i1 = int(1/(2*dy)*(L+s))        
i2 = int(1/(2*dy)*(L-s))   
i3 = int(1/(2*dy)*(L-s) - a/dy)

v = np.zeros((Ny,Ny), complex)
Ni = (Nx-2)*(Ny-2)

A = np.zeros((Ni,Ni), complex)
M = np.zeros((Ni,Ni), complex)


for k in range(Ni):     
    i = 1 + k//(Ny-2)
    j = 1 + k%(Ny-2)
    A[k,k] = 1 + 2*ct1 + 2*ct2 + 1j*dt/2*v[i,j]
    M[k,k] = 1 - 2*ct1 - 2*ct2 - 1j*dt/2*v[i,j]
    
    if i != 1: 
        A[k,(i-2)*(Ny-2)+j-1] = -ct2 
        M[k,(i-2)*(Ny-2)+j-1] = ct2
        
    if i != Nx-2:
        A[k,i*(Ny-2)+j-1] = -ct2
        M[k,i*(Ny-2)+j-1] = ct2
    
    if j != 1:
        A[k,k-1] = -ct1 
        M[k,k-1] = ct1 

    if j != Ny-2:
        A[k,k+1] = -ct1
        M[k,k+1] = ct1


from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

Asp = csc_matrix(A)

x = np.linspace(0, L, Ny-2)
y = np.linspace(0, L, Ny-2)
x, y = np.meshgrid(x, y)
psis = []

psi = psi0(x, y, x0, y0)
psi[0,:] = psi[-1,:] = psi[:,0] = psi[:,-1] = 0 
psi = double_walls(psi, j0, j1, i0, i1, i2, i3) 
psis.append(np.copy(psi))

for i in range(1,Nt):
    psi_vect = psi.reshape((Ni)) 
    b = np.matmul(M,psi_vect) 
    psi_vect = spsolve(Asp,b) 
    psi = psi_vect.reshape((Nx-2,Ny-2)) 
    psi = double_walls(psi, j0, j1, i0, i1, i2, i3)
    psis.append(np.copy(psi)) 

mod_psis = [] 

for wavefunc in psis:
    re = np.real(wavefunc) 
    im = np.imag(wavefunc)
    mod = np.sqrt(re**2 + im**2)
    mod_psis.append(mod)
fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0,L), ylim=(0,L))

img = ax.imshow(mod_psis[0], extent=[0,L,0,L], cmap=plt.get_cmap("gray"), vmin=0, vmax=np.max(mod_psis), zorder=1, interpolation="none")


bot = Rectangle((j0*dy,0),     w, i3*dy,      color="r", zorder=50) 
nid = Rectangle((j0*dy,i2*dy), w, (i1-i2)*dy, color="r", zorder=50)
top    = Rectangle((j0*dy,i0*dy), w, i3*dy,      color="r", zorder=50)

ax.add_patch(bot)
ax.add_patch(nid)
ax.add_patch(top)


def animate(i):
    img.set_data(mod_psis[i]) 
    img.set_zorder(1)
    return img,

anim = FuncAnimation(fig, animate, interval=1, frames=np.arange(0,Nt,2), repeat=False, blit=0)
plt.show()
anim.save('./slit_double.mp4', writer="ffmpeg", fps=30)

#-----------

Asp = csc_matrix(A)

x = np.linspace(0, L, Ny-2)
y = np.linspace(0, L, Ny-2)
x, y = np.meshgrid(x, y)
psis = []

psi = psi0(x, y, x0, y0)
psi[0,:] = psi[-1,:] = psi[:,0] = psi[:,-1] = 0 
psi = one_walls(psi, j0, j1, i0, i1, i2, i3) 
psis.append(np.copy(psi))

for i in range(1,Nt):
    psi_vect = psi.reshape((Ni)) 
    b = np.matmul(M,psi_vect) 
    psi_vect = spsolve(Asp,b) 
    psi = psi_vect.reshape((Nx-2,Ny-2)) 
    psi = one_walls(psi, j0, j1, i0, i1, i2, i3)
    psis.append(np.copy(psi)) 

mod_psis = [] 

for wavefunc in psis:
    re = np.real(wavefunc) 
    im = np.imag(wavefunc)
    mod = np.sqrt(re**2 + im**2)
    mod_psis.append(mod)
fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0,L), ylim=(0,L))

img = ax.imshow(mod_psis[0], extent=[0,L,0,L], cmap=plt.get_cmap("gray"), vmin=0, vmax=np.max(mod_psis), zorder=1, interpolation="none")


bot = Rectangle((j0*dy,0),     w, i2*dy,      color="r", zorder=50) #(x,y, w,h)
top    = Rectangle((j0*dy,i1*dy), w, i1*dy,      color="r", zorder=50)

ax.add_patch(bot)
ax.add_patch(top)


def animate(i):
    img.set_data(mod_psis[i]) 
    img.set_zorder(1)
    return img,

anim = FuncAnimation(fig, animate, interval=1, frames=np.arange(0,Nt,2), repeat=False, blit=0)
plt.show()
anim.save('./slit_one.mp4', writer="ffmpeg", fps=30)