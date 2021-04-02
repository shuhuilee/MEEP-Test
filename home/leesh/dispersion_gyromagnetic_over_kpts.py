import meep as mp 
import numpy as np 
import matplotlib.pyplot as plt 
import array as arr

#-------------------------------------------------------#
#Create gyrotropic medium:
#-------------------------------------------------------#
epsn, f0, gamma, sn, b0 = 1.5, 4.0, 4e-6, 0.4, 1.8
fcen, df, alpha= 0.2, 0.5, 1e-5
df1 = f0 - 1j*fcen*alpha
df2 = fcen + 1j*gamma
muperp = epsn + sn * df1/(df1**2 - df2**2)
xi = sn * df2 / (df1**2 - df2**2)

susc=[mp.GyrotropicSaturatedSusceptibility(frequency=f0, gamma=gamma, sigma=sn, bias=mp.Vector3(0,0,b0))]

mat=mp.Medium(epsilon = 1, mu = epsn, H_susceptibilities= susc) #gyroelectric medium; this is just a dot cell

#-------------------------------------------------------#
# Set up simulation cell:
#-------------------------------------------------------#

tmax, L = 100, 20
#tmax: number of time steps over which sim will run

cell, src_z, pml_layers = mp.Vector3(1,1,1), mp.Vector3(0,0,-10), [mp.PML(1.0,mp.Z)] 


#df: frequency width of source
#-------------------------------------------------------#
# Create light source:
#-------------------------------------------------------#
right_circ_src=[mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Hy, center=src_z), mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Hx, center=src_z, amplitude=1j)]
left_circ_src=[mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Hx, center=src_z), mp.Source(mp.GaussianSource(fcen, fwidth=df), component=mp.Hy, center=src_z, amplitude=1j)]
#-------------------------------------------------------#
# Set up k-values over which to plot w(k) frequencies:
#-------------------------------------------------------#
kz=mp.Vector3(0,0, 0.5)
kzm, n = mp.Vector3(0,0,0.02), 80
kpts=mp.interpolate(n, [kzm, kz])

#-------------------------------------------------------#
# Create and run sim:
#-------------------------------------------------------#
	#simulate for left-handed polarized light
sim_right_circ_src=mp.Simulation(cell_size=cell, geometry=[], sources=right_circ_src, default_material=mat, resolution=20)
allfreqs_right_circ=sim_right_circ_src.run_k_points(tmax, kpts)
print(allfreqs_right_circ, kpts)
print('len(allfreqs_right_circ)==len(kpts):', len(allfreqs_right_circ)==len(kpts))
freqs=[]


f_right_circ=[]
data_store=[]
for i in range(len(allfreqs_right_circ)):
    data_store+=[allfreqs_right_circ[i][0].real]
f_right_circ+=data_store




kpts_1=np.linspace(0.02, 0.5, len(f_right_circ))

fpts_right_circ=np.array(f_right_circ)


#######################

#simulate for left-handed polarized light
sim_left_circ_src=mp.Simulation(cell_size=cell, geometry=[], sources=left_circ_src, default_material=mat, resolution=20)
allfreqs_left_circ=sim_left_circ_src.run_k_points(tmax, kpts)
print(allfreqs_left_circ, kpts)
print('len(allfreqs_left_circ)==len(kpts):', len(allfreqs_left_circ)==len(kpts))
freqs=[]



f_left_circ=[]
data_store=[]
for i in range(len(allfreqs_left_circ)):
    data_store+=[allfreqs_left_circ[i][0].real]
f_left_circ+=data_store

print('f_right_circ:',f_right_circ)
print('f_left_circ:',f_left_circ)


#print('Checking... trueorfalse=[f_right_circ==f_left_circ]')
#trueorfalse=[f_left_circ[i]-0.1e-6<=f_right_circ[i]<=f_left_circ[i]+0.1e-6]
#print('Check: trueorfalse:', trueorfalse)
fpts_left_circ=np.array(f_left_circ)


#########


kz=np.linspace(0.02,0.5,80)
omegap=np.abs(kz)/np.sqrt(muperp+xi)
omegam=np.abs(kz)/np.sqrt(muperp-xi)

plt.figure()
ro=plt.scatter(kpts_1,fpts_right_circ, s=60, facecolors='None', edgecolors='g', label='Light in gyromedium: RCP')
lo=plt.scatter(kpts_1,fpts_left_circ, s=40, facecolors='None', edgecolors='r', label='Light in gyromedium: LCP')
plt.plot(kz, omegap, 'r', label='Theory: LHCP light')
plt.plot(kz, omegam, 'g', label='Theory: RHCP light')
plt.plot(kz, kz, 'b:', label='Theory: Light in vacuum')
plt.plot(kz, np.array([1/epsn**0.5 for i in range(len(kz))])*kz, color='y', label='Light in dielectric of index epsn**0.5')
plt.xlabel(r'k  ($ \frac{2 \pi}{a} $)')
plt.ylabel(r'$\omega$  ($\frac{2 \pi c}{a}$)')
plt.xlim(0.159, 0.19)
plt.ylim(0.125,0.145)
plt.legend()
plt.show()





