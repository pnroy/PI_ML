void MCRotLinStep(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp)
{

   int it0 = (it1 - 1);
   int it2 = (it1 + 1);

   if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   int t0 = offset + it0;
   int t1 = offset + it1;
   int t2 = offset + it2;

   double n1[NDIM];

   double cost = MCAngles[CTH][t1];
   double phi  = MCAngles[PHI][t1];

// cost += (step*(rnd1()-0.5));
// phi  += (step*(rnd1()-0.5));
   cost += (step*(rand1-0.5));
   phi  += (step*(rand2-0.5));

   if (cost >  1.0)
   {
      cost = 2.0 - cost;
//    phi  = phi + M_PI;
   }

   if (cost < -1.0)
   {
       cost = -2.0 - cost;
//     phi  = phi  + M_PI;
   }

   double sint = sqrt(1.0 - cost*cost);

   newcoords[AXIS_X][t1] = sint*cos(phi);
   newcoords[AXIS_Y][t1] = sint*sin(phi);
   newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// the old density

   double p0 = 0.0;
   double p1 = 0.0;

   for (int id=0;id<NDIM;id++)
   {
      p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   }

   double dens_old;
   double rho1,rho2,erot;

   if(RotDenType == 0)
   {
      dens_old = SRotDens(p0,type)*SRotDens(p1,type);
   }
   else if(RotDenType == 1)
   {
      rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
      rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
      dens_old = rho1+rho2;
   }

   if (fabs(dens_old)<RZERO) dens_old = 0.0;
   if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");

   double pot_old  = 0.0;

   int itr0 = it1  * RotRatio;     // interval to average over
   int itr1 = itr0 + RotRatio;     // translational time slices

   for (int it=itr0;it<itr1;it++)  // average over tr time slices
   pot_old  += (PotRotEnergy(gatom,MCCosine,it));

// the new density 

   p0 = 0.0;
   p1 = 0.0;


   for (int id=0;id<NDIM;id++)
   {
       p0 += (MCCosine [id][t0]*newcoords[id][t1]);
       p1 += (newcoords[id][t1]*MCCosine [id][t2]);
   }

   double dens_new;

   if(RotDenType == 0)
   {
      dens_new = SRotDens(p0,type)*SRotDens(p1,type);
   }
   else if(RotDenType == 1)
   {
      rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
      rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
      dens_new = rho1 + rho2;
   }

   if (fabs(dens_new)<RZERO) dens_new = 0.0;
   if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");

   double pot_new  = 0.0;

   for (int it=itr0;it<itr1;it++)  // average over tr time slices
   pot_new  += (PotRotEnergy(gatom,newcoords,it));

   double rd;

   if(RotDenType == 0)
   {
      if (dens_old>RZERO)
       rd = dens_new/dens_old;
      else rd = 1.0;

      rd *= exp(- MCTau*(pot_new-pot_old));
   }
   else if(RotDenType == 1)
   {
      rd = dens_new - dens_old - MCTau*(pot_new-pot_old);
//    rd = exp(rd);
   }

   bool Accepted = false;
   if(RotDenType == 0)
   {
      if (rd>1.0)         Accepted = true;
//    else if (rd>rnd7()) Accepted = true;
      else if (rd>rand3) Accepted = true;
   }
   else if (RotDenType == 1)
   {
      if (rd > 0.0)   Accepted = true;
//    else if (rd > log(rnd7())) Accepted = true;
      else if (rd > log(rand3)) Accepted = true;
   }

   MCRotChunkTot += 1.0;

   if (Accepted)
   {
      MCRotChunkAcp += 1.0;

      MCAngles[CTH][t1] = cost;
      MCAngles[PHI][t1] = phi;

      for (int id=0;id<NDIM;id++)
      MCCosine [id][t1] = newcoords[id][t1];
   }

}
void MCRotationsMove(int type) // update all time slices for rotational degrees of freedom
{
#ifdef DEBUG_PIMC
   const char *_proc_=__func__;    //  MCRotationsMove() 
   if (type != IMTYPE)
   nrerror(_proc_,"Wrong impurity type");

   if (NDIM != 3)
   nrerror(_proc_,"Rotational sampling for 3D systems only");
#endif

   double step   = MCAtom[type].rtstep; 
   int    offset = MCAtom[type].offset;
 
   int atom0  = 0;                   // only one molecular impurtiy
   offset    += (NumbTimes*atom0);   // the same offset for rotational
   int gatom  = offset/NumbTimes;    // and translational degrees of freedom

   double MCRotChunkTot = 0.0;
   double MCRotChunkAcp = 0.0;

   RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"
   double rand1,rand2,rand3;

/*
   for (int it1=0;it1<NumbRotTimes;it1++)
   {
      rand1=runif(Rng);
      rand2=runif(Rng);
      rand3=runif(Rng);
      MCRotLinStep(it1,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
   }
*/

/*
   #pragma omp parallel reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3)
   {
      int tid=omp_get_thread_num();
      int itini=chunksize*tid;
      int itfnl=itini+chunksize;
      for (int itrot=itini;itrot<itfnl-1;itrot++)
      {
         rand1=runif(Rng);
         rand2=runif(Rng);
         rand3=runif(Rng);
         MCRotLinStep(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
      }
   }  // end omp parallel

   for (int itrot=chunksize-1;itrot<NumbRotTimes;itrot=itrot+chunksize)
   {
      rand1=runif(Rng);
      rand2=runif(Rng);
      rand3=runif(Rng);
      MCRotLinStep(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
   }

   for (int itrot=NThreads*chunksize;itrot<NumbRotTimes;itrot++)
   {
      rand1=runif(Rng);
      rand2=runif(Rng);
      rand3=runif(Rng);
      MCRotLinStep(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
   }

   MCTotal[type][MCROTAT] += MCRotChunkTot;
   MCAccep[type][MCROTAT] += MCRotChunkAcp;
*/

   #pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3)
   for (int itrot=0;itrot<NumbRotTimes;itrot=itrot+2)
   {
      rand1=runif(Rng);
      rand2=runif(Rng);
      rand3=runif(Rng);
      MCRotLinStep(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
   }

   MCTotal[type][MCROTAT] += MCRotChunkTot;
   MCAccep[type][MCROTAT] += MCRotChunkAcp;

   MCRotChunkTot = 0;
   MCRotChunkAcp = 0;

   #pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3)
   for (int itrot=1;itrot<NumbRotTimes;itrot=itrot+2)
   {
      rand1=runif(Rng);
      rand2=runif(Rng);
      rand3=runif(Rng);
      MCRotLinStep(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
   }

   MCTotal[type][MCROTAT] += MCRotChunkTot;
   MCAccep[type][MCROTAT] += MCRotChunkAcp;

}
