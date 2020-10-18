from argparse import ArgumentParser
import os.path
from os import path
from rotor_pi import linrotstep
import numpy as np
import scipy as sp
import subprocess

RZERO=1e-10
T=1. # temperature in Kelvin
beta=1./T # in K^-1
P=16 # number of beads
tau=beta/float(P)
kB= 0.6950356  # 1/cm per Kelvin
B=1. # rotational constant in Kelvin
B*=kB # # rotational constant in 1/CM (required for linden.x)
mu=0.
Ngrid=1500

delta_gamma=2./float(Ngrid)

nskip=1
nskip2=1000

MC_steps=10000
step_z=.5
step_phi=2.
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
path_xyz=np.zeros((3,P),float)
path_xyz_new=np.zeros((3,P),float)

#create an off-diagonal rho
Nz=30
Nphi=30
dz=2./float(Nz)
dphi=2.*np.pi/float(Nphi)

Ngrid_zphi=Nz*Nphi
rho_eep=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
rho_E_eep=np.zeros((Ngrid_zphi,Ngrid_zphi),float)
for i in range(Nz):
	z=-1.+dz*float(i)
	sint = np.sqrt(1.0 - z*z)
	for j in range(Nphi):
		phi=dphi*float(j)
		x=sint*np.cos(phi)
		y=sint*np.sin(phi)
		for ip in range(Nz):
			zp=-1.+dz*float(ip)
			sintp = np.sqrt(1.0 - zp*zp)
			for jp in range(Nphi):
				phip=dphi*float(jp)
				xp=sintp*np.cos(phip)
				yp=sintp*np.sin(phip)
				dot=x*xp+y*yp+z*zp
				index=int((dot+1.)/delta_gamma)
				if (index>=Ngrid):
					index=Ngrid-1
				rho_eep[i*Nphi+j,ip*Nphi+jp]=rho[index]

print('building rho_eep done')


Ctau=np.zeros

N_accept=0
N_total=0
N_averages=0

E_avg = 0.0
E2_avg = 0.0

paths_output=open('paths.xyz','w')

P_list=[]
for p in range(0,P,2):
	P_list.append(p)
for p in range(1,P,2):
	P_list.append(p)

for step in range(MC_steps):

	#sample one bead at a time; can one improve this?
	for p in P_list:

		N_total,N_accept=linrotstep(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,mu,tau,N_total,N_accept,p,RZERO)
	
	if (step%nskip==0):
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
			dens=rho[index_1]

			E=rho_E[index_1]/dens
			E2=(rho_E[index_1]/dens)**2
			if (dens > RZERO):
				srot +=E
				srot2 +=E2
			#print(dot1,E)

		E_avg+=srot
		E2_avg+=srot2
# output paths to a file
	if (step%nskip2==0):
		for p in range(P):
			paths_output.write(str(path_xyz[0,p])+' '+str(path_xyz[1,p])+' '+str(path_xyz[2,p])+'\n')
		paths_output.write('\n')
paths_output.close()

print('accept ratio: ',float(N_accept)/float(N_total))
E_avg=E_avg/float(N_averages)
E2_avg=E2_avg/float(N_averages)
print('<K>: ',E_avg,np.sqrt((E2_avg-E_avg**2))/(float(N_averages)))


