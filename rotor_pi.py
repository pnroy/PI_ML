import numpy as np
def linrotstep(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,Ngrid,mu,tau,N_total,N_accept,p,RZERO):
   p_minus=p-1
   p_plus=p+1
   if (p_minus<0):
      p_minus+=P
   if (p_plus>=P):
      p_plus-=P

   z=path_angles[0,p]
   phi=path_angles[1,p]

   #uniform MC move
   z+=step_z*(np.random.random()-.5)
   phi+=step_phi*(np.random.random()-.5)


   if (z >  1.0):
      z = 2.0 - z
   if (z < -1.0):
      z = -2.0 - z
   sint = np.sqrt(1.0 - z*z)
   x=sint*np.cos(phi)
   y=sint*np.sin(phi)

# adjust pji in -Pi..Pi range
   phi=(np.arctan2(y,x))

   path_xyz_new[0,p] = x
   path_xyz_new[1,p] = y
   path_xyz_new[2,p] = z

# the old density
   dot0 = 0.0  
   dot1 = 0.0
   for id in range(3):
      dot0 += (path_xyz[id,p_minus]*path_xyz[id,p])
      dot1 += (path_xyz[id,p]*path_xyz[id,p_plus])

   # find points on cos_gamma grid  
   index_0=int((dot0+1.)/delta_gamma)
   index_1=int((dot1+1.)/delta_gamma)

   if (index_0>=Ngrid):
      print(index_0)
      index_0=Ngrid-1
   if (index_1>=Ngrid):
      print(index_1)
      index_1=Ngrid-1

   dens_old = rho[index_0]*rho[index_1]

   if (dens_old<RZERO): 
      dens_old = 0.0

   # potential
   pot_old=mu*path_xyz[2,p]

# the new density 

   dot0 = 0.0  
   dot1 = 0.0
   for id in range(3):
      dot0 += (path_xyz[id,p_minus]*path_xyz_new[id,p])
      dot1 += (path_xyz_new[id,p]*path_xyz[id,p_plus])
   index_0=int((dot0+1.)/delta_gamma)
   index_1=int((dot1+1.)/delta_gamma)

   if (index_0>=Ngrid):
      print(index_0)
      index_0=Ngrid-1
   if (index_1>=Ngrid):
      print(index_1)
      index_1=Ngrid-1

   dens_new=rho[index_0]*rho[index_1]

   if (dens_new<RZERO): 
      dens_new = 0.0
   pot_new=mu*path_xyz_new[2,p]


   if (dens_old>RZERO):
      rd = dens_new/dens_old
   else: 
      rd = 1.0

   rd *= np.exp(- tau*(pot_new-pot_old));

   Accepted = False

   if (rd>1.0):
      Accepted = True
   elif (rd>np.random.random()): 
      Accepted = True

   N_total+=1

   if (Accepted==True):
      N_accept+=1
      path_angles[0,p] = z
      path_angles[1,p] = phi

      for id in range(3):
         path_xyz[id,p]=path_xyz_new[id,p]
   return N_total,N_accept
   
def table_sample(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,Ngrid,rho_eep,dz,dphi,Nz,Nphi,mu,tau,N_total,N_accept,p,RZERO):
   p_minus=p-1
   p_plus=p+1
   if (p_minus<0):
      p_minus+=P
   if (p_plus>=P):
      p_plus-=P

# the old density

   i0=int((path_angles[0,p_minus]+1.)/dz)
   j0=int((path_angles[1,p_minus]+np.pi)/dphi)
   i1=int((path_angles[0,p]+1.)/dz)
   j1=int((path_angles[1,p]+np.pi)/dphi)
   i2=int((path_angles[0,p_plus]+1.)/dz)
   j2=int((path_angles[1,p_plus]+np.pi)/dphi)

   if (i0>=Nz):
      i0=Nz-1
   if (i1>=Nz):
      i1=Nz-1
   if (i2>=Nz):
      i2=Nz-1
   if (j0>=Nphi):
      j0=Nphi-1
   if (j1>=Nphi):
      j1=Nphi-1
   if (j2>=Nphi):
      j2=Nphi-1

   dens_old = rho_eep[i0*Nphi+j0,i1*Nphi+j1]*rho_eep[i1*Nphi+j1,i2*Nphi+j2]

# test below
# 100 grid points for z and phi each is good enough
   # dot0 = 0.0     
   # dot1 = 0.0
   # for id in range(3):
   #    dot0 += (path_xyz[id,p_minus]*path_xyz[id,p])
   #    dot1 += (path_xyz[id,p]*path_xyz[id,p_plus])

   # # find points on cos_gamma grid  
   # index_0=int((dot0+1.)/delta_gamma)
   # index_1=int((dot1+1.)/delta_gamma)

   # if (index_0>=Ngrid):
   #    index_0=Ngrid-1
   # if (index_1>=Ngrid):
   #    index_1=Ngrid-1

   # dens_old2 = rho[index_0]*rho[index_1]

   # print(rho_eep[i0*Nphi+j0,i1*Nphi+j1],rho[index_0])

   if (dens_old<RZERO): 
      dens_old = 0.0

   # potential
   pot_old=mu*path_xyz[2,p]

   z=path_angles[0,p]
   phi=path_angles[1,p]
   #uniform MC move
   z+=step_z*(np.random.random()-.5)
   phi+=step_phi*(np.random.random()-.5)

   # try integer sampling
   index_z=np.random.randint(-1,1)
   index_phi=np.random.randint(-1,1)

   if (z >  1.0):
      z = 2.0 - z
   if (z < -1.0):
      z = -2.0 - z

   sint = np.sqrt(1.0 - z*z)
   x=sint*np.cos(phi)
   y=sint*np.sin(phi)

   # adjust pji in -Pi..Pi range
   phi=(np.arctan2(y,x))

   i1_new=int((z+1.)/dz)
   j1_new=int((phi+np.pi)/dphi)
  # print(j1_new)

   if (i1_new>=Nz):
      i1_new=Nz-1
   if (j1_new>=Nphi):
      j1_new=Nphi-1

   path_xyz_new[0,p] = x
   path_xyz_new[1,p] = y
   path_xyz_new[2,p] = z

# the new density 

   dens_new=rho_eep[i0*Nphi+j0,i1_new*Nphi+j1_new]*rho_eep[i1_new*Nphi+j1_new,i2*Nphi+j2]

   if (dens_new<RZERO): 
      dens_new = 0.0
   pot_new=mu*path_xyz_new[2,p]

   if (dens_old>RZERO):
      rd = dens_new/dens_old
   else: 
      rd = 1.0

   rd *= np.exp(- tau*(pot_new-pot_old));

   Accepted = False

   if (rd>1.0):
      Accepted = True
   elif (rd>np.random.random()): 
      Accepted = True

   N_total+=1

   if (Accepted==True):
      N_accept+=1
      path_angles[0,p] = z
      path_angles[1,p] = phi

      for id in range(3):
         path_xyz[id,p]=path_xyz_new[id,p]
      i1=i1_new
      j1=j1_new
   return N_total,N_accept

def cosine_sample(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,Ngrid,mu,tau,N_total,N_accept,p,RZERO):
# under development
   p_minus=p-1
   p_plus=p+1
   if (p_minus<0):
      p_minus+=P
   if (p_plus>=P):
      p_plus-=P

   z=path_angles[0,p]
   phi=path_angles[1,p]

   #uniform MC move
   z+=step_z*(np.random.random()-.5)
   phi+=step_phi*(np.random.random()-.5)

   if (z >  1.0):
      z = 2.0 - z
   if (z < -1.0):
      z = -2.0 - z
   sint = np.sqrt(1.0 - z*z)
   x=sint*np.cos(phi)
   y=sint*np.sin(phi)

# adjust pji in -Pi..Pi range
#   phi=(np.arctan2(y,x))

   path_xyz_new[0,p] = x
   path_xyz_new[1,p] = y
   path_xyz_new[2,p] = z

# the old density
   dot0 = 0.0  
   dot1 = 0.0
   for id in range(3):
      dot0 += (path_xyz[id,p_minus]*path_xyz[id,p])
      dot1 += (path_xyz[id,p]*path_xyz[id,p_plus])

   # find points on cos_gamma grid  
   index_0=int((dot0+1.)/delta_gamma)
   index_1=int((dot1+1.)/delta_gamma)

   if (index_0>=Ngrid):
      print(index_0)
      index_0=Ngrid-1
   if (index_1>=Ngrid):
      print(index_1)
      index_1=Ngrid-1

   dens_old = rho[index_0]*rho[index_1]

   if (dens_old<RZERO): 
      dens_old = 0.0

   # potential
   pot_old=mu*path_xyz[2,p]

# the new density 

   dot0 = 0.0  
   dot1 = 0.0
   for id in range(3):
      dot0 += (path_xyz[id,p_minus]*path_xyz_new[id,p])
      dot1 += (path_xyz_new[id,p]*path_xyz[id,p_plus])
   index_0=int((dot0+1.)/delta_gamma)
   index_1=int((dot1+1.)/delta_gamma)

   if (index_0>=Ngrid):
      print(index_0)
      index_0=Ngrid-1
   if (index_1>=Ngrid):
      print(index_1)
      index_1=Ngrid-1

   dens_new=rho[index_0]*rho[index_1]

   if (dens_new<RZERO): 
      dens_new = 0.0
   pot_new=mu*path_xyz_new[2,p]


   if (dens_old>RZERO):
      rd = dens_new/dens_old
   else: 
      rd = 1.0

   rd *= np.exp(- tau*(pot_new-pot_old));

   Accepted = False

   if (rd>1.0):
      Accepted = True
   elif (rd>np.random.random()): 
      Accepted = True

   N_total+=1

   if (Accepted==True):
      N_accept+=1
      path_angles[0,p] = z
      path_angles[1,p] = phi

      for id in range(3):
         path_xyz[id,p]=path_xyz_new[id,p]
   return N_total,N_accept