import numpy as np
import matplotlib.pyplot as plt
import mm1

f=mm1.readfile("bragg.txt",0)
f=mm1.transpose(f)

# plt.plot(f[0],f[1],'+-',label='plane=100')
# plt.plot(f[2],f[3],'o-',label='plane=110')
# plt.xlabel('Angle (deg)')
# plt.ylabel('Amplitude(mV)')
# plt.grid()
# plt.legend()
# plt.show()
#for fixed crystal orientation
plt.plot(f[4],f[5],'o-')
plt.title('Fixed Crystal at 100-plane')
plt.xlabel('Angle (deg)')
plt.ylabel('Amplitude(mV)')
plt.grid()
plt.show()
