import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import os

x = np.linspace(-30,30,5000)
deltax = x[1]-x[0]


dirname = os.path.dirname(__file__)

def norm(phi):
    norm = np.sum(np.square(np.abs(phi)))*deltax
    return phi/np.sqrt(norm)

def complex_plot(x,y,prob=True,**kwargs):
    real = np.real(y)
    imag = np.imag(y)
    a,*_ = plt.plot(x,real,label='Re',linestyle='dashed',**kwargs)
    b,*_ = plt.plot(x,imag,label='Im',linestyle='dashed',**kwargs)
    plt.xlim(-2,2)
    if prob:
        p,*_ = plt.plot(x,np.abs(y),label='$\sqrt{P}$')
        return a,b,p
    else:
        return a,b
    
def wave_packet(pos=0,mom=0,sigma=0.2):
    a = (np.exp(-np.square(x-pos)/sigma/sigma,dtype=complex))
    return norm(np.exp(-1j*mom*x)*np.exp(-np.square(x-pos)/sigma/sigma,dtype=complex))

def d_dxdx(phi,x=x):
    dphi_dxdx = -2*phi
    dphi_dxdx[:-1] += phi[1:]
    dphi_dxdx[1:] += phi[:-1]
    return dphi_dxdx/deltax*deltax

def d_dt(phi,h=1,m=1,V=0):
    return 1j*h/2/m * d_dxdx(phi) - 1j*V*phi/h

def rk4(phi, dt, **kwargs):
    k1 = d_dt(phi, **kwargs)
    k2 = d_dt(phi+dt/2*k1, **kwargs)
    k3 = d_dt(phi+dt/2*k2, **kwargs)
    k4 = d_dt(phi+dt*k3, **kwargs)
    return phi + dt/6*(k1+2*k2+2*k3+k4)

def simulate(phi_sim,V=0, steps=100000, dt=1e-1, condition=None, normalize=True,save_every=100):
    simulation_steps = [np.copy(phi_sim)]
    for i in range(steps):
        phi_sim = rk4(phi_sim,dt,V=V)
        if condition:
            phi_sim = condition(phi_sim)
        if normalize:
            phi_sim = norm(phi_sim)
        if save_every is not None and (i+1) % save_every == 0:
            simulation_steps.append(np.copy(phi_sim))
    return simulation_steps

        
def animate(simulation_steps,init_func=None):
    fig = plt.figure()
    re,im,prob = complex_plot(x,simulation_steps[0])
    plt.xlim(-2,2)
    plt.ylim(-2,2)
    if init_func:
        init_func()
    plt.legend()

    def animate(frame):
        prob.set_data((x, np.abs(simulation_steps[frame])))
        re.set_data((x, np.real(simulation_steps[frame])))
        im.set_data((x, np.imag(simulation_steps[frame])))
        return prob,re,im
    
    anim = FuncAnimation(fig, animate, frames=int(len(simulation_steps)), interval=50)
    plt.close()

    return anim

def free_init():
        plt.xlim(-4,4)
        plt.ylim(-3,3)

def box_init():
        plt.gcf().axes[0].axvspan(3, 4, alpha=0.2, color='red')
        plt.gcf().axes[0].axvspan(-4, -3, alpha=0.2, color='red')
        plt.xlim(-4,4)
        plt.ylim(-3,3)

def barrier_init():
    plt.gcf().axes[0].axvspan(1.4, 1.6, alpha=0.2, color='orange')
    plt.xlim(-15,16)
    plt.ylim(-1.5,1.5)


#-----------------------------------------------------------------------------------------------
'''
sim_free = simulate(wave_packet(sigma=0.2),steps=30000,save_every=300)

anim = animate(sim_free)
writervideo = animation.FFMpegWriter(15)            
anim.save(os.path.join(dirname,"./free.mp4"), writer=writervideo)
'''

#-----------------------------------------------------------------------------------------------

sim_free = simulate(wave_packet(pos=-3,mom=-15,sigma=0.2),steps=30000,save_every=500)

anim = animate(sim_free, init_func=free_init)
writervideo = animation.FFMpegWriter(20)            
anim.save(os.path.join(dirname,"./free_p.mp4"), writer=writervideo)

#-----------------------------------------------------------------------------------------------
'''
box_potential = np.where((x>-3)&(x<3),0,1)
sim_box_mom = simulate(wave_packet(mom=15,sigma=0.5),V=box_potential,steps=70000,save_every=500)

anim = animate(sim_box_mom, init_func=box_init)
writervideo = animation.FFMpegWriter(18)            
anim.save(os.path.join(dirname,"./box.mp4"), writer=writervideo)
'''
#-----------------------------------------------------------------------------------------------
'''
barrier_weak_potential = np.where((x>1.4)&(x<1.6),2e-2,0)
sim_barrier_mom = simulate(wave_packet(mom=-15,sigma=0.5,pos=-5),V=barrier_weak_potential,steps=105050,save_every=500)

anim = animate(sim_barrier_mom,init_func=barrier_init)
writervideo = animation.FFMpegWriter(40)            
anim.save(os.path.join(dirname,"./barrier.mp4"), writer=writervideo)
'''
#-----------------------------------------------------------------------------------------------
'''
finite_well = np.where((x>1.5)&(x<2),-20e-3,0)
sim_finite_well = simulate(wave_packet(mom=-10),V=finite_well,steps=90000,save_every=500)

anim = animate(sim_finite_well,init_func=barrier_init)
writervideo = animation.FFMpegWriter(60)            
anim.save(os.path.join(dirname,"./finite_well.mp4"), writer=writervideo)
'''


