import meep as mp 
import numpy as np 
import matplotlib.pyplot as plt 
import array as arr
#-------------------------------------------------------#
#Create gyrotropic medium:
#-------------------------------------------------------#
epsn, f0, gamma, sn, b0 = 1.5, 1.0, 1e-6, 0.1, 0.15
#susc=[mp.GyrotropicLorentzianSusceptibility(frequency=f0, gamma=gamma, sigma=sn, bias=mp.Vector3(0,0,b0))]

#mat=mp.Medium(epsilon=epsn, mu=1, E_susceptibilities=susc)

#-------------------------------------------------------#
# Set up simulation cell:
#-------------------------------------------------------#
tmax, L = 100, 20
#tmax: number of time steps over which sim will run

cell, fcen, src_z, pml_layers = mp.Vector3(1,1,1), 0.8, mp.Vector3(0,0,-10), [mp.PML(1.0,mp.Z)] 

df=2.0
#df: frequency width of source
#-------------------------------------------------------#
# Create light source:
#-------------------------------------------------------#
src=[mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ex, center=src_z), mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Ey, center=src_z, amplitude=1j)]
#-------------------------------------------------------#
# Set up k-values over which to plot w(k) frequencies:
#-------------------------------------------------------#
kz=mp.Vector3(0,0,0.18)
kzm, n = mp.Vector3(0,0,-0.18), 30
kpts=mp.interpolate(n, [kzm, kz])

#-------------------------------------------------------#
# Create and run sim:
#-------------------------------------------------------#
	#simulate for left-handed polarized light
sim=mp.Simulation(cell_size=cell, geometry=[], sources=src, default_material=mp.vacuum, resolution=20)
allfreqs=sim.run_k_points(tmax, kpts)

print(allfreqs, kpts)
print('len(allfreqs)==len(kpts):', len(allfreqs)==len(kpts))
freqs=[]

for i in range(len(allfreqs)):
   freqs+=[allfreqs[i][4]] #freqs: take out first freq of each kpt, put into a list. 
   
print('len(freqs)==len(kpts):',len(freqs)==len(kpts))


dfsq = (f0**2 - 1j*fcen*gamma - fcen**2)
eperp = epsn + sn*f0**2*dfsq / (dfsq**2 - (fcen*b0)**2)
eta = sn * f0**2 * fcen * b0/ (dfsq**2- (fcen*b0)**2)


k=mp.interpolate(n, [-0.18,0.18]) #set up range of k_{z} values to plot
print('len(k)==len(kpts)',len(k)==len(kpts))
print('len(k)==len(freqs):', len(k)==len(freqs))
print('type(k)==type(freqs):',type(k)==type(freqs))

#swap, translate +1.0 freq
for i in range(len(freqs)):
    freqs[i]=-1*freqs[i]+1.48 - 0.065-0.0250 #add +1.0 for translate +1.0 freq

#print('swapped_+1.0translated_freqs:',freqs)
#

freqs_arr=np.array(freqs)
k_arr=np.array(k)
print('np.shape(freqs_arr)==np.shape(k_arr):',np.shape(freqs_arr)==np.shape(k_arr))


kz=np.linspace(-0.18,0.18,30)
omegap=np.abs(kz)/np.sqrt(eperp+eta)
omegam=np.abs(kz)/np.sqrt(eperp-eta)

plt.figure() 
plt.plot(k_arr, freqs_arr )
plt.plot(kz, omegap, 'r', label='Theory: LHCP light')
plt.plot(kz, omegam, 'g', label='Theory: RHCP light')
plt.legend()
plt.show()
#-------------------------------------------------------#
# Extract frequencies
#-------------------------------------------------------#



#import matplotlib.pyplot as plt 
#plt.figure()
#plt.plot(kpts_arr, freqs_arr)
#plt.plot(k, omegap, 'r', label='Theory: LHCP light')
#plt.plot(k, omegam, 'g', label='Theory: RHCP light')
#plt.legend([lo],['Expt: LHCP light'])
#plt.show()
#freq_arr, kspace_arr =arr.array('f', freq), np.linspace(-np.pi, np.pi, num=50)
# convert list of frequencies and k points into arrays for plotting

#-------------------------------------------------------#
# Calculate theoretical values:
#-------------------------------------------------------#
#dfsq = (f0**2 - 1j*fcen*gamma - fcen**2)
#eperp = epsn + sn*f0**2*dfsq / (dfsq**2 - (fcen*b0)**2)
#eta = sn * f0**2 * fcen * b0/ (dfsq**2- (fcen*b0)**2)


#k=np.linspace(-np.pi,np.pi,80) #set up range of k_{z} values to plot
#omegap=np.abs(k)/np.sqrt(eperp+eta)
#omegam=np.abs(k)/np.sqrt(eperp-eta)

#-------------------------------------------------------#
# Plot theory and experimental results
#-------------------------------------------------------#

#plt.figure()
#lo=plt.scatter(kspace_arr,freq_arr, s=70, facecolors='none', edgecolors='r')
#plt.plot(k, omegap, 'r', label='Theory: LHCP light')
#plt.plot(k, omegam, 'g', label='Theory: RHCP light')
#plt.legend([lo],['Meep: LHCP light'])
#plt.show()

