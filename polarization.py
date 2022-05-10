import numpy as np
import matplotlib.pyplot as plt
import mm1
from scipy.optimize import curve_fit
f=mm1.readfile("polarization.txt",1)
f=mm1.transpose(f)
print(len(f[0])==len(f[1]))
#a,b=np.polyfit(f[0],f[1],1)
#print(a,b)
# fit=[]
# for i in range(len(f[0])):
# 	fit.append(a*f[0][i]+b)
def model(x,a,b,c,d):
	x=x/180*np.pi
	return a*np.cos(b*x+c)**2+d

p,pcov=curve_fit(model,f[0],f[1])
for i in range(len(p)):
	p[i]=round(p[i],2)
y=[]
for i in range(len(f[0])):
	y.append(model(f[0][i],*p))
plt.plot(f[0],f[1],'ro-',label='data')
plt.plot(f[0],y,'b-')#,label=str(p[0])+'cos('+str(p[1])+'x+'+str(p[2])+')^2+'+str(p[3]))
#plt.plot(f[0],fit,'r-',label=str(round(a,4))+'x+'+str(round(b,4)))
plt.xlabel("Polarization angle (deg)")
plt.ylabel("Amplitude (mV)")
plt.grid()
plt.legend()
plt.show()