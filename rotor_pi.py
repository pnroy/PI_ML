import numpy as np
def linrotstep(P,path_angles,step_z,step_phi,path_xyz,path_xyz_new,delta_gamma,rho,mu,tau,N_total,N_accept,p,RZERO):
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