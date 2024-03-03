import matplotlib.pyplot as plt
import numpy as np
import os
plt.style.use('default')
dirname = os.path.dirname(__file__)
from importlib import reload
plt=reload(plt)
from celluloid import Camera

N = 99
l=2

x = []
for i in range(N+1):
    x.append(i*l/N)
    
f = open(os.path.join(dirname,"./src/T1.txt"), "r")
T1 = list(map(float,f.read().strip().split("\t")))


f = open(os.path.join(dirname,"./src/T1_an.txt"), "r")
T1_an = list(map(float,f.read().strip().split("\t")))


errT1 = []
for i in range(len(T1)):  
    errT1.append(abs(T1[i]-T1_an[i]))

"""
plt.scatter(x, T1, marker="x", color="black",s=40)
plt.savefig(os.path.join(dirname,"src/T1.jpg"))
plt.show()
plt.scatter(x, err, marker="x", color="black",s=40)
plt.savefig("src/T1_err.jpg")
plt.show()"""



##########################################################

f = open(os.path.join(dirname,"./src/T2.txt"), "r")
T2 = list(map(float,f.read().strip().split("\t")))
f.close()
f = open(os.path.join(dirname,"./src/T2_an.txt"), "r")
T2_an = list(map(float,f.read().strip().split("\t")))
f.close()

errT2 = []
for i in range(len(T2)):  
    errT2.append(abs(T2[i]-T2_an[i]))
    
"""plt.scatter(x, T2, marker="x", color="black",s=40)
plt.savefig("src/T2.jpg")
plt.show()
plt.scatter(x, err, marker="x", color="black",s=40)
plt.savefig("src/T2_err.jpg")
plt.show()"""

##########################################################

f = open(os.path.join(dirname,"./src/T3.txt"), "r")
T3 = list(map(float,f.read().strip().split("\t")))
f.close()

f = open(os.path.join(dirname,"./src/T3_an.txt"), "r")
T3_an = list(map(float,f.read().strip().split("\t")))
f.close()

errT3 = []
for i in range(len(T3)):  
    errT3.append(abs(T3[i]-T3_an[i]))

"""plt.scatter(x, T3, marker="x", color="black",s=40)
plt.savefig("src/T3.jpg")
plt.show()
plt.scatter(x, err, marker="x", color="black",s=40)
plt.savefig("src/T3_err.jpg")
plt.show()"""

#plot explicit
plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')


plt.scatter(x, T1, marker="x", color="black",s=30)
plt.savefig(os.path.join(dirname,"src/T1.jpg"))
plt.show()

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')


plt.scatter(x, T2, marker="s", color="black",s=30,label="$\Delta \hat{t}=0.49(\Delta \hat{x})^2$")
plt.scatter(x, T3, marker="D", color="cyan",s=10,label="$\Delta \hat{t}=0.25(\Delta \hat{x})^2$")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5)
plt.savefig(os.path.join(dirname,"src/T2_3.jpg"),bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, errT1, marker="x", color="black",s=30)
plt.savefig(os.path.join(dirname,"src/T1_err.jpg"))
plt.show()

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, errT2, marker="s", color="black",s=30,label="$\Delta \hat{t}=0.49(\Delta \hat{x})^2$")
plt.scatter(x, errT3, marker="D", color="cyan",s=30,label="$\Delta \hat{t}=0.25(\Delta \hat{x})^2$")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5)
plt.savefig(os.path.join(dirname,"src/errT2_3.jpg"),bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()



##########################################################

f = open(os.path.join(dirname,"./src/Ti1.txt"), "r")
Ti1 = list(map(float,f.read().strip().split("\t")))
f.close()

f = open(os.path.join(dirname,"./src/Ti1_an.txt"), "r")
Ti1_an = list(map(float,f.read().strip().split("\t")))
f.close()


erri1 = []
for i in range(len(Ti1)):  
    erri1.append(abs(Ti1[i]-Ti1_an[i]))

"""plt.scatter(x, Ti1, marker="x", color="black",s=40)
plt.savefig("src/Ti1.jpg")
plt.show()
plt.scatter(x, erri1, marker="x", color="black",s=40)
plt.savefig("src/Ti1_err.jpg")
plt.show()"""

##########################################################

f = open(os.path.join(dirname,"./src/Ti2.txt"), "r")
Ti2 = list(map(float,f.read().strip().split("\t")))
f.close()

f = open(os.path.join(dirname,"./src/Ti2_an.txt"), "r")
Ti2_an = list(map(float,f.read().strip().split("\t")))
f.close()

erri2 = []
for i in range(len(Ti2)):  
    erri2.append(abs(Ti2[i]-Ti2_an[i]))


"""plt.scatter(x, Ti2, marker="x", color="black",s=40)
plt.savefig("src/Ti2.jpg")
plt.show()
plt.scatter(x, erri2, marker="x", color="black",s=40)
plt.savefig("src/Ti2_err.jpg")
plt.show()"""

##########################################################
f = open(os.path.join(dirname,"./src/Ti3.txt"), "r")
Ti3 = list(map(float,f.read().strip().split("\t")))
f.close()

erri3 = []
for i in range(len(Ti2)):  
    erri3.append(abs(Ti3[i]-Ti2_an[i]))

#plot explicit 
plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, Ti1, marker="s", color="black",s=30,label="$\Delta \hat{t}=\Delta \hat{x}$")
plt.scatter(x, Ti2, marker="D", color="cyan",s=20,label="$\Delta \hat{t}=0.5\Delta \hat{x}$")
plt.scatter(x, Ti3, marker="v", color="orange",s=15,label="$\Delta \hat{t}=0.125\Delta \hat{x}$")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5)
plt.savefig(os.path.join(dirname,"src/Ti1_2_3.jpg"),bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, erri1, marker="s", color="black",s=30,label="$\Delta \hat{t}=\Delta \hat{x}$")
plt.scatter(x, erri2, marker="D", color="cyan",s=30,label="$\Delta \hat{t}=0.5\Delta \hat{x}$")
plt.scatter(x, erri3, marker="v", color="orange",s=30,label="$\Delta \hat{t}=0.125\Delta \hat{x}$")
lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), fancybox=True, shadow=True, ncol=5)
plt.savefig(os.path.join(dirname,"src/Tierr1_2_3.jpg"),bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()

dif_T3_Ti2=[]
for i in range(len(Ti2)):  
    dif_T3_Ti2.append(abs(Ti2[i]-T3[i]))

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, dif_T3_Ti2, marker="x", color="black",s=30)
plt.savefig(os.path.join(dirname,"src/dif_T3_Ti2.jpg"))
plt.show()

dif_Ti2_Ti3=[]
for i in range(len(Ti2)):  
    dif_Ti2_Ti3.append(abs(Ti2[i]-Ti3[i]))
    

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, dif_Ti2_Ti3, marker="x", color="black",s=30)
plt.savefig(os.path.join(dirname,"src/dif_Ti2_Ti3.jpg"))
plt.show()



##########################################################

f = open(os.path.join(dirname,"./src/TSol.txt"), "r")
TSol = list(map(float,f.read().strip().split("\t")))
f.close()

f = open(os.path.join(dirname,"./src/T_sol_an.txt"), "r")
T_sol_an = list(map(float,f.read().strip().split("\t")))
f.close()

f = open(os.path.join(dirname,"./src/TiSol.txt"), "r")
TiSol = list(map(float,f.read().strip().split("\t")))
f.close()

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, TSol, marker="x", color="black",s=30)
plt.axvline(x = 0.75, color = 'b')
plt.axhline(y = 50, color = 'r',linestyle='--')
plt.axvline(x = 1.25, color = 'b')
plt.axvspan(0, 0.75, alpha=0.1, color="red")
plt.axvspan(0.75, 1.25, alpha=0.1, color="blue")
plt.axvspan(1.25, 2, alpha=0.1, color="red")
plt.savefig(os.path.join(dirname,"src/TSol.jpg"))
plt.show()

TSol_err=[]
for i in range(len(Ti2)):  
    TSol_err.append(abs(TSol[i]-T_sol_an[i]))

plt.tick_params(axis="y", direction='in', length=5)
plt.tick_params(axis="x", direction='in', length=5)
plt.xlabel('$z$ (cm)')
plt.ylabel('$T$ (°C)')

plt.scatter(x, TSol_err, marker="x", color="black",s=30)
plt.savefig(os.path.join(dirname,"src/TSol_err.jpg"))
plt.show()


##########################################################

f = open(os.path.join(dirname,"./src/TSol_matrix.txt"), "r")
lines = f.readlines();
f.close()

TSol_m = np.zeros( (1537, 100) )
i=0
for line in lines:
   TSol_m[i]=(list(map(float,line.strip().split(",")[:-1])))
   i+=1

fig = plt.figure()
camera = Camera(fig)
for i in range(0,1537,16):
    p = plt.plot(x,TSol_m[i], 'bx')
    plt.tick_params(axis="y", direction='in', length=5)
    plt.tick_params(axis="x", direction='in', length=5)
    plt.xlabel('$z$ (cm)')
    plt.ylabel('$T$ (°C)')
    t=(i*0.13*(1/100)**2)*3.5e4/60
    plt.legend(p,['t={0:.3g} min'.format(t)])
    camera.snap()

k=0
while(k<7):
    p = plt.plot(x,TSol_m[i], 'bx')
    plt.tick_params(axis="y", direction='in', length=5)
    plt.tick_params(axis="x", direction='in', length=5)
    plt.xlabel('$z$ (cm)')
    plt.ylabel('$T$ (°C)')
    t=(1566*0.13*(1/100)**2)*3.5e4/60
    plt.legend(p,['t={0:.3g} min'.format(t)])
    camera.snap()
    k+=1


animation = camera.animate()
animation.save(os.path.join(dirname,"./src/animacio.gif"))


