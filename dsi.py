import numpy as np
import matplotlib.pyplot as plt
import mm1

f=mm1.readfile("dsi.txt",1)
f=mm1.transpose(f)
a,b=np.polyfit(f[0],f[1],1)
print(a,b)
# fit=[]
# for i in range(len(f[0])):
# 	fit.append(a*f[0][i]+b)
plt.plot(f[0],f[1],'o-',label='data')
#plt.plot(f[0],fit,'r-',label=str(round(a,4))+'x+'+str(round(b,4)))
plt.xlabel("Angle (deg)")
plt.ylabel("Amplitude (mV)")
plt.grid()
#plt.legend()
plt.show()