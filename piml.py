from argparse import ArgumentParser
import os.path
from os import path
import sys
from rotor_pi import linrotstep,table_sample
import numpy as np
#import scipy as sp
import subprocess


if (len(sys.argv) <8):
        print('wrong number of argument')
        print('usage: python3 ',sys.argv[0],' <T> <P> <Nsteps> <Neq> <Nskip> <Nskip2> <sampling flag: 0 old method, 1 table sampling>')
        sys.exit()

RZERO=1e-10
T=float(sys.argv[1]) # temperature in Kelvin
beta=1./T # in K^-1
P=int(sys.argv[2]) # number of beads
tau=beta/float(P)
kB= 0.6950356  # 1/cm per Kelvin
B=1. # rotational constant in Kelvin
B*=kB # # rotational constant in 1/CM (required for linden.x)
mu=0.
Ngrid=1500
delta_gamma=2./float(Ngrid-1)
MC_steps=int(sys.argv[3])
n_equilibrate=int(sys.argv[4])
nskip=int(sys.argv[5])
nskip2=int(sys.argv[6])
sampling_flag=sys.argv[7]

step_z=1.
step_phi=1.5

linprop_command='./linear_prop/linden.x'

if (path.exists(linprop_command)==False):
	subprocess.run('make',cwd="./linear_prop")

#check if executable exists
#arguments=str(T)+' '+str(P)+' '+str(B)+' '+str(Ngrid)+' -1'
process_log=subprocess.run([linprop_command,str(T),str(P),str(B),str(Ngrid),' -1'])
#process_log=subprocess.run([linprop_command,str(T),str(P),str(B),str(Ngrid),' -1'], capture_output=True)
#print(process_log)
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

#initialize xyz coordinates
path_xyz=np.zeros((3,P),float)
path_xyz_new=np.zeros((3,P),float)
for p in range(P):
	z=path_angles[0,p]
	phi=path_angles[1,p]
	sint = np.sqrt(1.0 - z*z)
	x=sint*np.cos(phi)
	y=sint*np.sin(phi)
	path_xyz[0,p] = x
	path_xyz[1,p] = y
	path_xyz[2,p] = z
	path_xyz_new[0,p] = x
	path_xyz_new[1,p] = y
	path_xyz_new[2,p] = z

#create an off-diagonal rho
#Nz=100
#Nphi=100
Nz=100
Nphi=100
dz=2./float(Nz)
dphi=(2.*np.pi)/float(Nphi)

if (sampling_flag=='1'):
	step_z=int(.3/dz)
	step_phi=int(1./dphi)


path_angles_indices=np.zeros((2,P),int) # has to match angles for now
i=int((path_angles[0,p]+1.)/dz)
j=int((path_angles[0,p]+np.pi)/dphi)
for p in range(P):
	path_angles_indices[0,p]=i
	path_angles_indices[1,p]=j

rho_eep='None'
# Ngrid_zphi=Nz*Nphi
# rho_eep=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
# #rho_E_eep=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
# x_grid=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
# y_grid=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
# z_grid=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
# for i in range(Nz):
# 	z=-1.+dz*float(i)
# 	sint = np.sqrt(1.0 - z*z)
# 	for j in range(Nphi):
# 		phi=-np.pi+dphi*float(j)
# 		x=sint*np.cos(phi)
# 		y=sint*np.sin(phi)
# 		x_grid[i,j]=x
# 		y_grid[i,j]=y
# 		z_grid[i,j]=z

# for i in range(Nz):
# 	for j in range(Nphi):
# 		for ip in range(Nz):
# 			for jp in range(Nphi):
# 				dot=x_grid[i,j]*x_grid[ip,jp]+y_grid[i,j]*y_grid[ip,jp]+z_grid[i,j]*z_grid[ip,jp]
# 				index=int((dot+1.)/delta_gamma)
# 				if (index>=Ngrid):
# 					index=Ngrid-1
# 				rho_eep[i*Nphi+j,ip*Nphi+jp]=rho[index]

# print('building rho_eep done')

#Ctau=np.zeros

N_accept=0
N_total=0
N_averages=0

step_out=0

E_avg = 0.0
E2_avg = 0.0

paths_output=open('paths.xyz','w')
E_output=open('E.dat','w')

histo_out=open('histo.xyz','w')
nbins=20
dx=2./float(nbins)
xyz_histo=np.zeros((3,nbins),float)

P_list=[]
for p in range(0,P,2):
	P_list.append(p)
for p in range(1,P,2):
	P_list.append(p)

for step in range(MC_steps):

	#sample one bead at a time; can one improve this?
	for p in P_list:

		if (sampling_flag=='0'):
			N_total,N_accept=linrotstep(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,Ngrid,mu,tau,N_total,N_accept,p,RZERO)
		if (sampling_flag=='1'):
			N_total,N_accept=table_sample(P,path_angles_indices,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,Ngrid,rho_eep,dz,dphi,Nz,Nphi,mu,tau,N_total,N_accept,p,RZERO)

	if (step%nskip==0 and step >n_equilibrate):
		#estimators
		N_averages+=1
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

			if (index_1>=Ngrid):
				index_1=Ngrid-1

			dens=rho[index_1]

			E=rho_E[index_1]/dens
			E2=(rho_E[index_1]/dens)**2
			E_output.write(str(E)+' '+str(E2)+'\n')
			if (dens > RZERO):
				srot +=E
				srot2 +=E2
			#print(dot1,E)
			for i in range(3):
				index=int((path_xyz[i,p]+1.)/dx)
				if (index>=nbins):
					index=nbins-1
				xyz_histo[i,index]+=1.

		E_avg+=srot
		E2_avg+=srot2
# output paths to a file
	if (step%nskip2==0):
		for p in range(P):
			step_out+=1
			paths_output.write(str(step_out)+' '+str(path_xyz[0,p])+' '+str(path_xyz[1,p])+' '+str(path_xyz[2,p])+' '+str(path_angles[0,p])+' '+str(path_angles[1,p])+' '+str(path_angles_indices[0,p])+' '+str(path_angles_indices[1,p])+'\n')			
		paths_output.write('\n')
paths_output.close()

for index in range(nbins):
	histo_out.write(str(-1.+float(index)*dx)+' ')	
	for i in range(3):
		histo_out.write(str(xyz_histo[i,index]/float(N_total)/dx)+' ')
	histo_out.write('\n')
histo_out.close()
E_output.close()

print('accept ratio: ',float(N_accept)/float(N_total))
E_avg=E_avg/float(N_averages)
E2_avg=E2_avg/float(N_averages)
print('<E> (in K): ',E_avg,np.sqrt((E2_avg-E_avg**2))/(float(N_averages)))


