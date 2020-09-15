#straight waveguide example

import meep as mp
import numpy as np
import matplotlib.pyplot as plt

resolution = 50
sx=20
sy=15


cell = mp.Vector3(sx,sy,0)
dpml=0.5
pml_layers=[mp.PML(dpml)]
vertwaveguide=mp.Block(size=mp.Vector3(4,sy,0), center=mp.Vector3(0,0,0))
metal1=mp.Block(size=mp.Vector3(sx/2-2,sy,0), center=mp.Vector3(-((sx/2-2)/2+2), 0,0), material=mp.metal)
metal2=mp.Block(size=mp.Vector3(sx/2-2,sy,0), center=mp.Vector3((sx/2-2)/2+2, 0,0), material=mp.metal)
metal3=mp.Block(size=mp.Vector3(sx,3,0),center=mp.Vector3(0,sy/2-1.5,0), material=mp.metal)
horizwaveguide=mp.Block(size=mp.Vector3(sx,4,0), center=mp.Vector3(0,sy/2-5,0))

mu_r=2.8
f0, gamma, sn, b0 = 6 , 35e-3/(2*np.pi), 1, 10
fcen, df, alpha= 6.68, 8.268, 0.22



susc=[mp.GyrotropicSaturatedSusceptibility(frequency=f0, gamma=gamma, sigma=sn, alpha=alpha, bias=mp.Vector3(0,0,b0))]
mat=mp.Medium(epsilon=1, mu=mu_r, H_susceptibilities=susc)

#gyro=mp.Cylinder(radius=1, material=mat, center=mp.Vector3(0,sy/2-5,0))

sourcelocation=-6

#beam_x0, beam_kdir, beam_w0, beam_E0 = mp.Vector3(0,sourcelocation,0), mp.Vector3(0,1,0), 1e-15, mp.Vector3(1,0,0)
source=mp.Source(mp.GaussianSource(frequency=fcen,width=1/fcen), center=mp.Vector3(0,sourcelocation, 0), component=mp.Hz, size=mp.Vector3(1,0,0))
                                 
geometry=[vertwaveguide, metal1,metal2,metal3, horizwaveguide]
sim=mp.Simulation(cell_size=cell, boundary_layers=pml_layers, geometry=geometry, sources=[source], resolution=resolution )
nfreq=100

infrontofsrc_fr=mp.FluxRegion(center=mp.Vector3(0, sourcelocation+0.5,0), size=mp.Vector3(4*2,0,0))
infrontofsrc=sim.add_flux(fcen, df, nfreq, infrontofsrc_fr)

refl_fr=mp.FluxRegion(center=mp.Vector3(0,sourcelocation-0.5, 0), size=mp.Vector3(4*2,0,0))
refl_before=sim.add_flux(fcen, df, nfreq, refl_fr)


refl_fr2 = mp.FluxRegion(center=mp.Vector3(0,sourcelocation+0.5,0), size=mp.Vector3(4*2,0,0))
refl_after = sim.add_flux(fcen, df, nfreq, refl_fr2)
# transmitted flux
tran_fr = mp.FluxRegion(center=mp.Vector3(sx/2-dpml-2.5,sy/2-5, 0), size=mp.Vector3(0,2*4,0))
tran = sim.add_flux(fcen, df, nfreq, tran_fr)

tran_fr2 = mp.FluxRegion(center=mp.Vector3(-(sx/2-dpml-2.5),sy/2-5, 0), size=mp.Vector3(0,2*4,0))
tran2= sim.add_flux(fcen, df, nfreq, tran_fr)

pt = mp.Vector3(sx/2-dpml-0.5,sy/2-5,0)

sim.use_output_directory()

sim.run(until_after_sources=mp.stop_when_fields_decayed(10,mp.Ez,pt,1e-3))
#mp.at_beginning(mp.output_epsilon), mp.at_every(0.1, mp.output_png(mp.Ez, "-Zc dkbluered")), 

fluxfreqs=mp.get_flux_freqs(refl_before)
reflectedflux_after=mp.get_fluxes(refl_after)
reflectedflux_before=mp.get_fluxes(refl_before)
infrontofsrcflux=mp.get_fluxes(infrontofsrc)
tran_flux = mp.get_fluxes(tran)
tran2_flux = mp.get_fluxes(tran2)

refl, tran=[],[]
freq, refl_after_fl, refl_before_fl, tran2 = [],[],[],[]
infrontofsrc=[]
for i in range(nfreq):
    freq=np.append(freq, fluxfreqs[i])
    infrontofsrc=np.append(infrontofsrc, infrontofsrcflux[i])
    refl_after=np.append(refl_after_fl, reflectedflux_after[i])
    refl_before=np.append(refl_before_fl, reflectedflux_before[i])
    tran=np.append(tran, tran_flux[i])
    tran2 = np.append(tran2, tran2_flux[i])


#----------------------------------
delta_n=f0**2-freq-1j*freq*gamma
chi_perp=(f0**2)*delta_n*sn/(delta_n-(freq**2)*b0**2)
chi_par=(f0**2) * sn/ delta_n
eta= f0**2 * freq* sn /(delta_n**2-freq**2)
#----------------------------------

incidentflux=infrontofsrc-refl_before
plt.figure()
plt.plot(freq,-refl_before/incidentflux, label='Reflectance')
plt.plot(freq, tran/incidentflux, color='green',label='Transmittance (Right arm)')

plt.plot(freq, tran2/incidentflux, ':b', linewidth=3, label='Transmittance (Left arm)')
plt.plot(freq, 1-(-refl_before/incidentflux)-tran2/incidentflux-tran/incidentflux, ':y', linewidth=3.5
, label='Loss')

#plt.plot(freq,chi_perp+1, label=r'$\mu_{11}$')
#plt.plot(freq,eta, label=r'$\mu_{12}$')
plt.xlabel(r'Frequency (GHz)')
#plt.ylabel(r'(a.u.)')

plt.legend(loc='lower right',borderaxespad=0.)

plt.show()







#######################

df1 = f0 - 1j*fcen*alpha
df2 = fcen + 1j*gamma
muperp = mu_r + sn * df1/(df1**2 - df2**2)
xi = sn * df2 / (df1**2 - df2**2)

tmax = 1/f0

kz=mp.Vector3(0,0, 0.5)
kzm, n = mp.Vector3(0,0,0.02), 80
kpts=mp.interpolate(n, [kzm, kz])
#simulate for left-handed polarized light
sim_left_circ_src=mp.Simulation(cell_size=cell, geometry=geometry, sources=[source], default_material=mat, resolution=20)
allfreqs_left_circ=sim_left_circ_src.run_k_points(tmax, kpts)
print(allfreqs_left_circ, kpts)
print('len(allfreqs_left_circ)==len(kpts):', len(allfreqs_left_circ)==len(kpts))
freqs=[]



f=[]
data_store=[]
for i in range(len(allfreqs_left_circ)):
    data_store+=[allfreqs_left_circ[i][0].real]
f+=data_store


print('f:',f)


#print('Checking... trueorfalse=[f_right_circ==f_left_circ]')
#trueorfalse=[f_left_circ[i]-0.1e-6<=f_right_circ[i]<=f_left_circ[i]+0.1e-6]
#print('Check: trueorfalse:', trueorfalse)
fpts_left_circ=np.array(f)


#########


kz=np.linspace(0.02,0.5,80)
omegap=np.abs(kz)/np.sqrt(muperp+xi)
omegam=np.abs(kz)/np.sqrt(muperp-xi)

kpts_1=np.linspace(0.02, 0.5, len(f))

plt.figure()
ro=plt.scatter(kpts_1,f, s=60, facecolors='None', edgecolors='g', label='Light in gyromedium')
plt.plot(kz, omegap, 'r', label='Theory: LHCP light')
plt.plot(kz, omegam, 'g', label='Theory: RHCP light')
plt.plot(kz, kz, 'b:', label='Theory: Light in vacuum')
plt.plot(kz, np.array([1/mu_r**0.5 for i in range(len(kz))])*kz, color='y', label='Light in dielectric of index epsn**0.5')
plt.xlabel(r'k  ($ \frac{2 \pi}{a} $)')
plt.ylabel(r'$\omega$  ($\frac{2 \pi c}{a}$)')
plt.xlim(0.159, 0.19)
plt.ylim(0.125,0.145)
plt.legend()
plt.show()