from argparse import ArgumentParser
import os.path
from os import path
import numpy as np
import scipy as sp
import subprocess
import matplotlib.pyplot as plt


T=1. # temperature in Kelvin
beta=1./T # in K^-1
P=4 # number of beads
tau=beta/float(P)
B=1. # rotational constant
Ngrid=1500

delta_gamma=2./float(Ngrid)

nskip=100

MC_steps=100000
step_z=1.
step_phi=1.
linprop_command='./linear_prop/linden.x'

if (path.exists(linprop_command)==False):
	subprocess.run('make',cwd="./linear_prop")

#check if executable exists
#arguments=str(T)+' '+str(P)+' '+str(B)+' '+str(Ngrid)+' -1'
process_log=subprocess.run([linprop_command,str(T),str(P),str(B),str(Ngrid),' -1'], capture_output=True)
print(process_log)
density_file_path='linden.out'
density_file=open(density_file_path,'r')
cos_gamma=np.zeros(Ngrid,float)
rho=np.zeros(Ngrid,float)
rho_E=np.zeros(Ngrid,float)
rho_E2=np.zeros(Ngrid,float)
# read-in density matrix and energy estimnators
i=0
for line in density_file:
	data=line.split()
	cos_gamma[i]=float(data[0])
	rho[i]=np.fabs(float(data[1]))
	rho_E[i]=float(data[2])
	rho_E2[i]=float(data[3])
	i+=1

# define a path
path_angles=np.zeros((2,P),float)
path_xyz=np.zeros((3,P),float)
path_xyz_new=np.zeros((3,P),float)

N_accept=0
N_total=0

E_avg = 0.0
E2_avg = 0.0

paths_output=open('paths.xyz','w')

for step in range(MC_steps):

	#sample one bead at a time; can one improve this?
	for p in range(P):

		z=path_angles[0,p]
		phi=path_angles[1,p]

		#uniform MC move
		z+=step_z*2.*(np.random.random()-.5)
		phi+=step_phi*2.*np.pi*(np.random.random())

		if (z >  1.0):
			z = 2.0 - z
		if (z < -1.0):
			z = -2.0 - z
		sint = np.sqrt(1.0 - z*z)

		x=sint*np.cos(phi)
		y=sint*np.sin(phi)

		path_xyz_new[0,p] = x
		path_xyz_new[1,p] = y
		path_xyz_new[2,p] = z

# the old density

		dot0 = 0.0	
		dot1 = 0.0
		p_minus=p-1
		p_plus=p+1
		if (p_minus<0):
			p_minus+=P
		if (p_plus>=P):
			p_plus-=P
		for id in range(3):
			dot0 += (path_xyz[id,p_minus]*path_xyz[id,p])
			dot1 += (path_xyz[id,p]*path_xyz[id,p_plus])

		# find points on cos_gamma grid
		
		index_0=int((dot0+1.)/delta_gamma)
		index_1=int((dot1+1.)/delta_gamma)

		dens_old = rho[index_0]*rho[index_1]

		# potential
#
#		pot_old=0.

# the new density 

		dot0 = 0.0	
		dot1 = 0.0
		for id in range(3):
			dot0 += (path_xyz[id,p_minus]*path_xyz_new[id,p])
			dot1 += (path_xyz_new[id,p]*path_xyz[id,p_plus])
		index_0=int((dot0+1.)/delta_gamma)
		index_1=int((dot1+1.)/delta_gamma)
		dens_new=rho[index_0]*rho[index_1]

		rd=dens_new/dens_old

#		pot_new=0.
#		rd *= np.exp(- tau*(pot_new-pot_old));
   
		Accepted = False

		if (rd>1.0):
			Accepted = True;
		elif (rd>np.random.random()): 
			Accepted = True;

		N_total+=1

		if (Accepted==True):
			N_accept+=1
			path_angles[0,p] = z
			path_angles[1,p] = phi

			for id in range(3):
				path_xyz[id,p]=path_xyz_new[id,p]

	if (step%nskip==0):
		#estimators
		srot=0.
		srot2=0.
		for p in range(P):

			dot1 = 0.0
			p_plus=p+1
			if (p_plus>=P):
				p_plus-=P
			for id in range(3):
				dot1 += (path_xyz[id,p]*path_xyz[id,p_plus])

			index_1=int((dot1+1.)/delta_gamma)
			dens=rho[index_1]
			E=rho_E[index_1]/dens
			E2=(rho_E[index_1]/dens)**2
			if (dens > 1.e-9):
				srot +=E
				srot2 +=E2
			#print(dot1,E)

		E_avg+=srot
		E2_avg+=srot2
# output paths to a file
		for p in range(P):
			paths_output.write(str(path_xyz[0,p])+' '+str(path_xyz[1,p])+' '+str(path_xyz[2,p])+'\n')
		paths_output.write('\n')
paths_output.close()

print('accept ratio: ',float(N_accept)/float(N_total))
E_avg=E_avg/float(MC_steps)*float(nskip)
E2_avg=E2_avg/float(MC_steps)*float(nskip)
print('<K>: ',E_avg,np.sqrt((E2_avg-E_avg**2))/(float(MC_steps)/float(nskip)))


