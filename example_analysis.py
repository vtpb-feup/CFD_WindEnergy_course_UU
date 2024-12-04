import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt

a=loadmat('simulation_2600.mat')

XT=[500,1500,2500,3500,4500]
X=a['x_u'][0]
Y=a['y_u'][0]
# Z=a['z_u'][0]
U=(a['u'][31,31,:]+a['u'][32,31,:]+a['u'][31,32,:]+a['u'][32,32,:])/4 # average of the 4 points so simplify the interpolation
Umiddle=(a['u'][31,:,:]+a['u'][32,:,:])/2 

UTurbines=np.interp(XT,X,U)
CT=0.8
a_AD=(1-np.sqrt(1-CT))/2
PowerTurbines=0.5*1.225*np.pi*50**2*UTurbines**3/(1-a_AD)**2*CT/1e6
plt.subplot(211)
plt.plot(X,U)
plt.plot(XT,UTurbines,'ro')
plt.xlabel('x [m]')
plt.ylabel('u [m/s]')
plt.title('u vs x')
plt.grid()

plt.subplot(212)
plt.plot(XT,PowerTurbines,'ok')
plt.xlabel('x [m]')
plt.ylabel('Power [MW]')
plt.title('Power vs x')
plt.tight_layout()
plt.grid()

plt.figure()
plt.contourf(X,Y,Umiddle,30,cmap='jet')
plt.colorbar()  
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('u vs x,y')
plt.show()

# convergence analysis
tt=np.arange(1,27)*100
turbine=2
for t in tt:
    a=loadmat('simulation_'+str(t)+'.mat')
    U=(a['u'][31,31,:]+a['u'][32,31,:]+a['u'][31,32,:]+a['u'][32,32,:])/4
    UTurbines=np.interp(XT[turbine],X,U)
    PowerTurbines=0.5*1.225*np.pi*50**2*UTurbines**3/(1-a_AD)**2*CT/1e6
    plt.plot(t,PowerTurbines,'ok')
plt.xlabel('iteration')
plt.ylabel('Power [MW]')
plt.title('Power vs time')
plt.grid()
plt.show()

