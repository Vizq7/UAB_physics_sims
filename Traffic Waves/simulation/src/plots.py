import matplotlib.pyplot as plt
import numpy as np
import os
plt.style.use('default')
dirname = os.path.dirname(__file__)

def column(matrix, i):
    return [row[i] for row in matrix]


N=21000
dt=0.001
'''
t = []
for i in range(N+1):
    t.append(i*dt)
    
+t2 = []
for i in range(int(1/dt),int(21/dt)):
    t2.append(i*dt)
'''

C1_tau0_x = np.genfromtxt(os.path.join(dirname,'./C1_tau0_x.txt'),dtype=float,delimiter=',')

t = []
for i in range(len(column(C1_tau0_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C1_tau0_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C1_tau0_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau0_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau0_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau0_x, 4), marker="x", color="lime",s=1, label="4")
plt.xlim([0,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("C1_tau0_x.jpg",dpi=fig.dpi)

C1_tau0_v = np.genfromtxt(os.path.join(dirname,'./C1_tau0_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C1_tau0_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C1_tau0_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau0_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau0_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau0_v, 4), marker="x", color="lime",s=1, label="4")
plt.xlim([0,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C1_tau0_v.jpg",dpi=fig.dpi)

C1_tau0_err = np.genfromtxt(os.path.join(dirname,'./C1_tau0_err.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),int(21/dt)):
    t.append(i*dt)

a = column(C1_tau0_err, 0)
fig = plt.figure()
plt.scatter(t, column(C1_tau0_err, 0), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau0_err, 1), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau0_err, 2), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau0_err, 3), marker="x", color="lime",s=1, label="4")
plt.xlim([1,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C1_tau0_err.jpg",dpi=fig.dpi)


C2_tau0_x = np.genfromtxt(os.path.join(dirname,'./C2_tau0_x.txt'),dtype=float,delimiter=',')

t = []
for i in range(len(column(C2_tau0_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C2_tau0_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C2_tau0_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C2_tau0_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C2_tau0_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C2_tau0_x, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("C2_tau0_x.jpg")

C2_tau0_v = np.genfromtxt(os.path.join(dirname,'./C2_tau0_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C2_tau0_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C2_tau0_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C2_tau0_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C2_tau0_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C2_tau0_v, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C2_tau0_v.jpg",dpi=fig.dpi)

C2_tau0_err = np.genfromtxt(os.path.join(dirname,'./c2_tau0_err_total.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),len((C2_tau0_err))+int(1/dt)):
    t.append(i*dt)

fig = plt.figure()
plt.scatter(t, (C2_tau0_err), marker="x", color="orange",s=1)
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C2_tau0_err.jpg",dpi=fig.dpi)


C3_tau0_x = np.genfromtxt(os.path.join(dirname,'./C3_tau0_x.txt'),dtype=float,delimiter=',')

t = []
for i in range(len(column(C3_tau0_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C3_tau0_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C3_tau0_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C3_tau0_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C3_tau0_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C3_tau0_x, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
fig.savefig("C3_tau0_x.jpg")

C3_tau0_v = np.genfromtxt(os.path.join(dirname,'./C3_tau0_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C3_tau0_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C3_tau0_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C3_tau0_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C3_tau0_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C3_tau0_v, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C3_tau0_v.jpg",dpi=fig.dpi)

C3_tau0_err = np.genfromtxt(os.path.join(dirname,'./c3_tau0_err_total.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),len((C3_tau0_err))+int(1/dt)):
    t.append(i*dt)

fig = plt.figure()
plt.scatter(t, (C3_tau0_err), marker="x", color="orange",s=1)
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C3_tau0_err.jpg",dpi=fig.dpi)

C1_tau_x = np.genfromtxt(os.path.join(dirname,'./C1_tau_x.txt'),dtype=float,delimiter=',')

t = []
for i in range(len(column(C1_tau_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C1_tau_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C1_tau_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau_x, 4), marker="x", color="lime",s=1, label="4")
plt.xlim([0,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("C1_tau_x.jpg")

C1_tau_v = np.genfromtxt(os.path.join(dirname,'./C1_tau_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C1_tau_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C1_tau_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau_v, 4), marker="x", color="lime",s=1, label="4")
plt.xlim([0,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C1_tau_v.jpg",dpi=fig.dpi)

C1_tau_err = np.genfromtxt(os.path.join(dirname,'./C1_tau0_err.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),len((column(C1_tau_err, 0)))+int(1/dt)):
    t.append(i*dt)

fig = plt.figure()
plt.scatter(t, column(C1_tau_err, 0), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C1_tau_err, 1), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C1_tau_err, 2), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C1_tau_err, 3), marker="x", color="lime",s=1, label="4")
plt.xlim([1,10])
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C1_tau_err.jpg",dpi=fig.dpi)

C2_tau_x = np.genfromtxt(os.path.join(dirname,'./C2_tau_x.txt'),dtype=float,delimiter=',')
t = []
for i in range(len(column(C2_tau_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C2_tau_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C2_tau_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C2_tau_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C2_tau_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C2_tau_x, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("C2_tau_x.jpg",dpi=fig.dpi)

C2_tau_v = np.genfromtxt(os.path.join(dirname,'./C2_tau_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C2_tau_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C2_tau_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C2_tau_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C2_tau_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C2_tau_v, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C2_tau_v.jpg",dpi=fig.dpi)

C2_tau_err = np.genfromtxt(os.path.join(dirname,'./c2_tau_err_total.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),len((C2_tau_err))+int(1/dt)):
    t.append(i*dt)

fig = plt.figure()
plt.scatter(t, (C2_tau_err), marker="x", color="orange",s=1)
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C2_tau_err.jpg",dpi=fig.dpi)


C3_tau_x = np.genfromtxt(os.path.join(dirname,'./C3_tau_x.txt'),dtype=float,delimiter=',')
t = []
for i in range(len(column(C3_tau_x, 0))):
    t.append(i*dt)


fig = plt.figure()
plt.scatter(t, column(C3_tau_x, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C3_tau_x, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C3_tau_x, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C3_tau_x, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C3_tau_x, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("C3_tau_x.jpg")

C3_tau_v = np.genfromtxt(os.path.join(dirname,'./C3_tau_v.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(C3_tau_v, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(C3_tau_v, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(C3_tau_v, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(C3_tau_v, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(C3_tau_v, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("C3_tau_v.jpg",dpi=fig.dpi)

C3_tau_err = np.genfromtxt(os.path.join(dirname,'./c3_tau_err_total.txt'),dtype=float,delimiter=',')
t = []
for i in range(int(1/dt),len((C3_tau_err))+int(1/dt)):
    t.append(i*dt)

fig = plt.figure()
plt.scatter(t, (C3_tau_err), marker="x", color="orange",s=1)
plt.xlabel('$t$ (s)')
plt.ylabel('$\epsilon_\mathrm{r}$')
fig.savefig("C3_tau_err.jpg",dpi=fig.dpi)


t = []
for i in range(N+1):
    t.append(i*dt)

x_IDM = np.genfromtxt(os.path.join(dirname,'./x_IDM.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(x_IDM, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(x_IDM, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(x_IDM, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(x_IDM, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(x_IDM, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$x$ (m)')
fig.savefig("x_IDM.jpg")

v_IDM = np.genfromtxt(os.path.join(dirname,'./v_IDM.txt'),dtype=float,delimiter=',')

fig = plt.figure()
plt.scatter(t, column(v_IDM, 0), marker="x", color="black",s=1, label="0")
plt.scatter(t, column(v_IDM, 1), marker="x", color="green",s=1, label="1")
plt.scatter(t, column(v_IDM, 2), marker="x", color="red",s=1, label="2")
plt.scatter(t, column(v_IDM, 3), marker="x", color="cyan",s=1, label="3")
plt.scatter(t, column(v_IDM, 4), marker="x", color="lime",s=1, label="4")
lgnd = plt.legend()
lgnd.legendHandles[0]._sizes = [30]
lgnd.legendHandles[1]._sizes = [30]
lgnd.legendHandles[2]._sizes = [30]
lgnd.legendHandles[3]._sizes = [30]
lgnd.legendHandles[4]._sizes = [30]
plt.xlabel('$t$ (s)')
plt.ylabel('$v$ (km/h)')
fig.savefig("v_IDM.jpg",dpi=fig.dpi)


# --- ERR------
'''
c2_err = np.genfromtxt(os.path.join(dirname,'./C2_tau0_err.txt'),dtype=float,delimiter=',')

print(len(t2))
print(len(c2_err))
print(int(1/dt))
print(int(21/dt)+1)

fig = plt.figure()
plt.scatter(t2, c2_err, marker="x", color="black",s=1, label="0")
fig.savefig("c2_err.jpg",dpi=fig.dpi)

'''