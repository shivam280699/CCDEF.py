c ABFMCCDEF.  CCDEF, A SIMPLIFIED COUPLED CHANNEL CODE FOR FUSION CROSS   
c 1   SECTIONS INCLUDING STATIC NUCLEAR DEFORMATIONS.                     
c 2   J. FERNANDEZ-NIELLO, C.H. DASSO, S. LANDOWNE.                       
c REF. IN COMP. PHYS. COMMUN. 54 (1989) 409                               
c     CCDEF :  a simplified coupled-channels code for calculating fusion  
c     cross sections and angular momentum distributions including static  
c     deformations in projectile and target                             
c                                                                       
c     coded by j. fernandez niello, c.h. dasso and s. landowne          
c     ------------------------------------------------------------------
c                                                                       
c     input data:                                                       
c                                                                       
c     read *,wma,wza,wmb,wzb,dv,fcc,be2a,be4a,be2b,be4b                 
c         wma = projectile mass                                         
c         wza = projectile charge                                       
c         wmb = target mass                                             
c         wzb = target charge                                           
c         dv  = parameter which can be used to adjust the barrier       
c               parameters somewhat.  a value dv=0 is recommended       
c               to start.  typical ranges for changes:  -10<dv<10       
c         fcc = 0 :diagonalization at barrier position rb               
c               1 :(recommended for strong couplings) exponential       
c               character of the form factors taken into account.       
c               second order estimation of the effective barriers       
c               within a one-fermi distance from rb                     
c         be2a = static  quadrupole  deformation  for  projectile       
c         be4a = static hexadecapole deformation  for  projectile       
c         be2b = static  quadrupole  deformation  for the target        
c         be4b = static hexadecapole deformation  for the target        
c                                                                       
c                                                                       
c     read *, emin,emax,de                                              
c         emin = minimum value of the energy for cross section          
c         emax = maximum value of the energy for cross section          
c         de   = interval in the energy scale                           
c                                                                       
c     read *,ns,na                                                      
c         ns = number of surface inelastic channels to be included      
c              in the calculation                                       
c         na = number of additional channels to be included in the      
c              calculation                                              
c                                                                       
c     read *,beta,flam,q        (only if ns.ne.0, read ns times)        
c         beta = deformation parameter for the mode                     
c         flam = multipolarity of the mode (it should be entered        
c                negative for projectile  modes and positive for        
c                target modes)                                          
c         q    = q-value for the channel (negative)                     
c                                                                       
c     read *,f,q                (only if na.ne.0, read na times)        
c         f = strength of the coupling                                  
c         q = q-value for the channel  (with corresponding sign)        
c                                                                       
c                                                                       
c     note: * entering a single value of the energy (i.e. emin.eq.emax, 
c             de.ne.0)  the angular momentum cross section distribution 
c             is constructed                                            
c           * a  static  deformation  characterized  by  quadrupole and 
c             hexadecapole  parameters  be2, be4,  can be attributed to 
c             either the projectile or the target  by entering the mass 
c             of the given nucleus as a negative number. if both masses 
c             are positive, the values of be2, be4 are ignored          
c           * maximum number of energies is 51                          
c           * maximum number of channels (na+ns) is 10                  
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      implicit real*8(a-h,o-z)                                          
      dimension beta (10),fslam(10),qs(10),fs(10)                       
      dimension fla(10,2),pa(10,2),dfla(10,2),dfla2(10,2),              
     *          sum(100),sum0(100),sigl(100),sigl0(100),np(10),n1(10)   
      common/one/v0r,a0r,rab,redm,au,hc,etak                            
      common/two/ra,rb,b2y2a,b2y2b,b4y4a,b4y4b                          
      alog(xxx)=dlog(xxx)                                               
      sqrt(xxx)=dsqrt(xxx)                                              
      cbrt(xx)=xx**(1./3.)                                              
      cos(xxx)=dcos(xxx)                                                
      sin(xxx)=dsin(xxx)                                                
      float(iii)=dfloat(iii)
      
      open (unit=19,file='input.dat',status='unknown')
      
   20 format(/,'  Projectile:  mass=',f5.0,'          charge=',f3.0)    
   21 format(/,'  Parameters for DV=',f6.2,//,                          
     *'      Vb=',f6.1,'       Rb=',f6.2,'       h-Omega=',f5.2,/)      
   22 format('  Target:      mass=',f5.0,'          charge=',f3.0)      
   23 format(//,2x,35('*'),/,'  *    beta    * lamda *   Q(MeV)   *',/  
     *,2x,35('*'))
   24 format('  *  ',f6.2,'    *  ',f 3.0,'  *  ',f6.2,'    *')         
   25 format(2x,35('*'),/)                                              
   26 format(/,'  Additional strength:  enter F, Q',/)                  
   27 format(/,2x,35('*'),/,'  *     F(MeV)     *     Q(MeV)     *',/,  
     *2x,35('*'))                                                       
   28 format('  *    ',f6.2,'      *    ',f6.2,'      *')               
   29 format(2x,35('*'),/)                                              
   30 format(//,'  Fusion cross sections vs. energy   (in mb)')         
   31 format(/,'  ',42('*'),/,'  *  E(MeV)  *    Coupled   *',          
     *       '   Uncoupled  *',/,'  ',42('*'))                          
   32 format('  *  ',f6.1,'  *  ',1pe10.3e2,'  *  ',1pe10.3e2,'  *')    
   33 format('  ',42('*'))                                              
   34 format(/,'  Cross sections for E =',1f6.1,' MeV are:',/,          
     *'  Coupled =',1pe10.3e2,' mb       Uncoupled =',1pe10.3e2,' mb',/)
   35 format(/,'  Partial wave cross sections (in mb/h-bar)',/)         
   36 format(/,'  ',42('*'),/,'  *     l    *    Coupled   *',          
     *       '   Uncoupled  *',/,'  ',42('*'))                          
   37 format('  *   ',i3,'    *  ',1pe10.3e2,'  *  ',1pe10.3e2,'  *')   
   38 format(/,'  Moments of the partial wave distribution:',/)         
   39 format('  Coupled:       <l> =',1f5.1,'     sigma =',1f5.1)       
   40 format('  Uncoupled:     <l> =',1f5.1,'     sigma =',1f5.1,/)     
   41 format(/,'  ',27('*'),/,'  *  E(MeV)  *',                         
     *       '   Uncoupled  *',/,'  ',27('*'))                          
   42 format('  *  ',f6.1,'  *  ',1pe10.3e2,'  *')                      
   43 format('  ',27('*'))                                              
   44 format(/,'  Cross section for E =',1f6.1,' MeV is:',1pe10.3e2,    
     *         ' mb',/)                                                 
   46 format(/,'  ',27('*'),/,'  *     l    *',                         
     *       '   Uncoupled  *',/,'  ',27('*'))                          
   47 format('  *   ',i3,'    *  ',1pe10.3e2,'  *  ',1pe10.3e2,'  *')   
   48 format('  Uncoupled:     <l> =',1f5.1,'     sigma =',1f5.1,/)     
   49 format('               beta2=',f4.2,'          beta4=',f4.2)      
      do 50 k=1,100                                                     
      sigl0(k)=0.                                                       
      sigl(k)=0.                                                        
      sum0(k)=0.                                                        
   50 sum(k)=0.                                                         
      read (19,*)wma,wza,wmb,wzb,dv,fcc,be2a,be4a,be2b,be4b                 
      idefa=0                                                           
      idefb=0                                                           
      if(wma.lt.0.) idefa=1                                         
      if(wmb.lt.0.) idefb=1                                             
      wma=abs(wma)            
      wmb=abs(wmb)                                                      
      print 20,wma,wza                                                  
      if(idefa.ne.0) print 49,be2a,be4a                                 
      print 22,wmb,wzb                                                  
      if(idefb.ne.0) print 49,be2b,be4b                                 
      read (19,*)emin,emax,de                                               
      read (19,*)ns,na                                                      
      nmax=ns+na                                                        
      ns1=ns+1                                                          
      ne=(emax-emin)/de+1.001                                           
      if(ne.gt.51) stop '  Too many energies...'                        
      ra=1.233*cbrt(wma)-0.978/cbrt(wma)                            
      rb=1.233*cbrt(wmb)-0.978/cbrt(wmb)                                
      rred=ra*rb/(ra+rb)                                                
      rab=ra+rb+0.29                                                    
      redm=wma*wmb/(wma+wmb)                                            
      pi=4.*datan(1.d0)                                                 
      au=931.5016                                                       
      hc=197.3286                                                       
      h2m=hc*hc/au                                                      
      fpi=3.544908                                                      
      a0r=0.63                                                          
      etak=1.43997*wza*wzb                                              
      pdws=0.                                                           
      ithdma=1                                                          
      if(idefa.ne.0) ithdma=90                                          
      do 200 ithda=1,ithdma                                             
      theda=pi*ithda/180.                                               
      pdwa=sin(theda)*pi/180.                                           
      if(idefa.eq.0) then                                               
      pdwa=1.                                                           
      be2a=0.                                                           
      be4a=0.                                                           
      endif                                                             
      xsqa=cos(theda)**2                                                
      b2y2a=be2a*sqrt(5./4./pi)*.5*(3.*xsqa-1.)                         
      b4y4a=be4a*sqrt(9./4./pi)*(35.*xsqa*xsqa-30.*xsqa+3.)/8.          
      ithdmb=1                                                          
      if(idefb.ne.0) ithdmb=90                                          
      do 205 ithdb=1,ithdmb                                             
      thedb=pi*ithdb/180.                                               
      pdwb=sin(thedb)*pi/180.                                           
      if(idefb.eq.0) then                                               
      pdwb=1.                                                           
      be2b=0.                                                           
      be4b=0.                                                           
      endif                                                             
      xsqb=cos(thedb)**2                                                
      b2y2b=be2b*sqrt(5./4./pi)*.5*(3.*xsqb-1.)                         
      b4y4b=be4b*sqrt(9./4./pi)*(35.*xsqb*xsqb-30.*xsqb+3.)/8.          
      pdw=pdwa*pdwb                                                     
      v0rx=30.08*(1.-1.8*(1.-2.*wza/wma)*(1.-2.*wzb/wmb))               
      v0r=v0rx*rred+dv-20.                                              
      call bar(rbar,vb,homega)                                          
      if(idefa.eq.0.and.idefb.eq.0) print 21, dv,vb,rbar,homega         
      eps=homega/6.283185                                               
      rcal=rbar                                                         
      if(ns.eq.0.or.homega.lt.0.) go to 110                             
      if(ithda.eq.1.and.ithdb.eq.1) print 23                            
      do 100 n=1,ns                                                     
      if(ithda.eq.1.and.ithdb.eq.1) then                                
      read (19,*) beta(n),fslam(n),qs(n)                                    
      print 24,beta(n),fslam(n),qs(n)                                   
      endif                                                             
      rrt=rb                                                            
      if(fslam(n).lt.0.) rrt=ra                                         
      flam=abs(fslam(n))                                                
      q=qs(n)                                                           
      flam1=flam-1                                                      
      call pot(rcal,ur,dur,ddur)                                        
      fs(n)=beta(n)*rrt*(-dur+3.*etak*(rrt/rcal)**flam1/(2.*flam        
     1 +1.)/rcal**2)/fpi                                                
      f=fs(n)                                                           
      fkap=-ddur/dur                                                    
      sq=sqrt(q*q+4.*f*f)                                               
      fla(n,1)=.5*(-q-sq)                                               
      fla(n,2)=.5*(-q+sq)                                               
      dfla(n,1)=2.*fkap*f**2/sq                                         
      dfla(n,2)=-dfla(n,1)                                              
      dfla2(n,1)=-4.*(fkap*f)**2/sq+8.*(fkap*f*f)**2/sq**3              
      dfla2(n,2)=-dfla2(n,1)                                            
      pa(n,1)=f*f/(f*f+fla(n,1)**2)                                     
      pa(n,2)=f*f/(f*f+fla(n,2)**2)                                     
  100 continue                                                          
      if(ithda.eq.1.and.ithdb.eq.1) print 25                            
  110 continue                                                          
      if(na.eq.0.or.homega.lt.0.) go to 130                             
      if(ithda.eq.1.and.ithdb.eq.1) print 27                            
      do 120 n=ns1,nmax                                                 
      fkap=0.71                                                         
      if(ithda.eq.1.and.ithdb.eq.1) then                                
      read (19,*) fs(n),qs(n)                                               
      print 28,fs(n),qs(n)                                              
      endif                                                             
      q=qs(n)                                                           
      f=fs(n)                                                           
      sq=sqrt(q*q+4.*f*f)                                               
      fla(n,1)=.5*(-q-sq)                                               
      fla(n,2)=.5*(-q+sq)                                               
      pa(n,1)=f*f/(f*f+fla(n,1)**2)                                     
      pa(n,2)=f*f/(f*f+fla(n,2)**2)                                     
      dfla(n,1)=2.*fkap*f**2/sq                                         
      dfla(n,2)=-dfla(n,1)                                              
      dfla2(n,1)=-4.*(fkap*f)**2/sq+8.*(fkap*f*f)**2/sq**3              
      dfla2(n,2)=-dfla2(n,1)                                            
  120 continue                                                          
      if(ithda.eq.1.and.ithdb.eq.1) print 29                            
  130 continue                                                          
      if (homega.lt.0.) go to 201                                       
      ilim=nmax                                                         
      if(nmax.gt.10) stop '  Maximum number of channels exceeded'       
      nd=2**nmax                                                        
      np(1)=1                                                           
      do 140 n=2,nmax                                                   
      np(n)=np(n-1)*2                                                   
  140 continue                                                          
      do 190 i1=1,nd                                                    
      ic1=i1-1                                                          
      p=1.                                                              
      fl=0.                                                             
      dfl=0.                                                            
      dfl2=0.                                                           
      ick=0                                                             
      do 150 n=1,nmax                                                   
      n1(n)=ic1/np(nmax+1-n)                                            
      n1t=n1(n)+1                                                       
      ick=ick+n1(n)                                                     
      if(ick.gt.ilim) go to 190                                         
      p=p*pa(n,n1t)                                                     
      fl=fl+fla(n,n1t)                                                  
      dfl=dfl+dfla(n,n1t)                                               
      dfl2=dfl2+dfla2(n,n1t)                                            
      ic1=ic1-n1(n)*np(nmax+1-n)                                        
  150 continue                                                          
      delta=fcc*dfl/(redm*homega**2/h2m-dfl2)                           
      if(abs(delta).gt.0.99) delta=-.99*fl/abs(fl)                      
      vbw=vb-.5*redm*(homega*delta)**2/h2m+fl+dfl*delta+.5*dfl2*delta**2
      rcald=rcal+delta                                                  
      facw=31.416*rcal**2*eps                                           
      facwd=31.416*rcald**2*eps                                         
      do 180 ie=1,ne                                                    
      e=emin+de*float(ie-1)                                             
      aqa=(e-vbw)/eps                                                   
      if(aqa.lt.30.) aqa=alog(1.+exp(aqa))                              
      sum(ie)=sum(ie)+facwd*p*pdw*aqa/e                                 
      if(i1.eq.1) then                                                  
      aqa=(e-vb)/eps                                                    
      if(aqa.lt.30.) aqa=alog(1.+exp(aqa))                              
      sum0(ie)=sum0(ie)+pdw*facw*aqa/e                                  
      endif                                                             
      if(ne.ne.1) go to 170                                             
      do 160 il=1,100                                                   
      gl=float(il-1)                                                    
      vbwl=vbw+0.5*hc**2*gl*(gl+1.)/(redm*au*rcald**2)                  
      vbl=vb+0.5*hc**2*gl*(gl+1.)/(redm*au*rcal**2)                     
      factor=31.41592*(2.*gl+1.)*hc**2/(2.*redm*au*e)                   
      tl=1.0                                                            
      aux=exp((e-vbwl)/eps)                                             
      if(aux.lt.30.) tl=exp(aux)/(1.+exp(aux))                          
      sigl(il)=sigl(il)+pdw*factor*p*tl                                 
      tl0=1.0                                                           
      aux0=exp((e-vbl)/eps)                                             
      if(aux0.lt.30.) tl0=exp(aux0)/(1.+exp(aux0))                      
      if(i1.eq.1) sigl0(il)=sigl0(il)+pdw*factor*tl0                    
  160 continue                                                          
  170 continue                                                          
  180 continue                                                          
  190 continue                                                          
  201 continue                                                          
      pdws=pdws+pdw                                                     
  205 continue                                                          
  200 continue                                                          
      if(ne.eq.1) go to 220                                             
      print 30                                                          
      if(nmax.ne.0) then                                                
      print 31                                                          
      else                                                              
      print 41                                                          
      endif                                                             
      do 210 ie=1,ne                                                    
      sig0=sum0(ie)/pdws                                                
      sig=sum(ie)/pdws                                                  
      e=emin+de*float(ie-1)                                             
      write(60,*) e,sig                                                 
      if(nmax.ne.0) then                                                
      print 32,e,sig,sig0                                               
      else                                                              
      print 42,e,sig                                                    
      endif                                                             
  210 continue                                                          
      if(nmax.ne.0) then                                                
      print 33                                                          
      else                                                              
      print 43                                                          
      endif                                                             
      go to 260                                                         
  220 continue                                                          
      sig=sum(1)/pdws                                                   
      sig0=sum0(1)/pdws                                                 
      if(nmax.ne.0) then                                                
      print 34,emin,sig,sig0                                            
      else                                                              
      print 44,emin,sig                                                 
      endif                                                             
      print 35                                                          
      if(nmax.ne.0) then                                                
      print 36                                                          
      else                                                              
      print 46                                                          
      endif                                                             
      do 230 il=1,100                                                   
      l=il-1                                                            
      if(sigl(il).lt.1.2e-6) go to 240                                  
      xsl=sigl(il)/pdws                                                 
      xsl0=sigl0(il)/pdws                                               
      if(nmax.ne.0) then                                                
      print 37,l,xsl,xsl0                                               
      else                                                              
      print 47,l,xsl                                                    
      endif                                                             
  230 continue                                                          
  240 continue                                                          
      if(nmax.ne.0) then                                                
      print 33                                                          
      else                                                              
      print 43                                                          
      endif                                                             
      s0=0.                                                             
      s1=0.                                                             
      s2=0.                                                             
      su0=0.                                                            
      su1=0.                                                            
      su2=0.                                                            
      do 250 il=1,100                                                   
      flo=float(il-1)                                                   
      s0=s0+sigl(il)                                                    
      s1=s1+flo*sigl(il)                                                
      s2=s2+flo**2*sigl(il)                                             
      su0=su0+sigl0(il)                                                 
      su1=su1+flo*sigl0(il)                                             
      su2=su2+flo**2*sigl0(il)                                          
  250 continue                                                          
      avl=s1/s0                                                         
      sd=sqrt((s2/s0)-avl**2)                                           
      avl0=su1/su0                                                      
      sd0=sqrt((su2/su0)-avl0**2)                                       
      print 38                                                          
      if(nmax.ne.0) then                                                
      print 39,avl,sd                                                   
      print 40,avl0,sd0                                                 
      else                                                              
      print 48,avl,sd                                                   
      endif                                                             
  260 continue                                                          
      end                                                               
c                                                                       
c     ------------------------------------------------------------------
      subroutine bar(rbar,vb,homega)                                    
      implicit real*8(a-h,o-z)                                          
      common/one/v0r,a0r,rab,redm,au,hc,etak                            
c      common/two/ra,rb,b2y2a,b2y2b,b4y4a,b4y4b                         
      sqrt(xxx)=dsqrt(xxx)                                              
      rmax=20.                                                          
      dr=.5                                                             
      s=-1.                                                             
      y1=-1.                                                            
      rbz=rmax                                                          
  270 call potent(rbz,v0,v1,v2)                                         
      y2=v1                                                             
      if(y2*y1.gt.0.) go to 280                                         
      if(dr.lt.0.01) go to 300                                          
      dr=.5*dr                                                          
      s=-s                                                              
      y1=y2                                                             
  280 rbz=rbz+s*dr                                                      
      if(rbz.gt.rmax.or.rbz.lt.0.8) go to 290                           
      go to 270                                                         
  290 continue                                                          
      homega=-1.                                                        
      return                                                            
  300 continue                                                          
      homega=hc*sqrt(-v2/(redm*au))                                     
      vb=v0                                                             
      rbar=rbz                                                          
      return                                                            
      end                                                               
c                                                                       
c     ------------------------------------------------------------------
      subroutine potent(rx,v0,v1,v2)                                    
      implicit real*8(a-h,o-z)                                          
      common/one/v0r,a0r,rab,redm,au,hc,etak                            
      common/two/ra,rb,b2y2a,b2y2b,b4y4a,b4y4b                          
      x=rx                                                              
      x2=x*x                                                            
      x3=x2*x                                                           
      x4=x3*x                                                           
      r2a=ra*ra                                                         
      r4a=r2a**2                                                        
      r2b=rb*rb                                                         
      r4b=r2b**2                                                        
      v0=etak/x*(1.+.6*b2y2a*r2a/x2+b4y4a*r4a/x4/3.                     
     1   +.6*b2y2b*r2b/x2+b4y4b*r4b/x4/3.)                              
      v1=-etak/x2*(1.+1.8*b2y2a*r2a/x2+5.*b4y4a*r4a/x4/3.               
     1   +1.8*b2y2b*r2b/x2+5.*b4y4b*r4b/x4/3.)                          
      v2=+2*etak/x3*(1.+3.6*b2y2a*r2a/x2+15.*b4y4a*r4a/x4/3.            
     1   +3.6*b2y2b*r2b/x2+15.*b4y4b*r4b/x4/3.)                         
      call pot(x,vn,v1n,v2n)                                            
      v0=v0+vn                                                          
      v1=v1+v1n                                                         
      v2=v2+v2n                                                         
      return                                                            
      end                                                               
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
      subroutine pot(rr,ur,dur,ddur)                                    
      implicit real*8(a-h,o-z)                                          
      common/one/v0r,a0r,rab,redm,au,hc,etak                            
      common/two/ra,rb,b2y2a,b2y2b,b4y4a,b4y4b                          
      ss=rr-rab-ra*(b2y2a+b4y4a)-rb*(b2y2b+b4y4b)                       
      arg=exp(-ss/a0r)                                                  
      arg1=1.+arg                                                       
      ur=-v0r*arg/arg1                                                  
      dur=-ur/arg1/a0r                                                  
      ddur=dur*(1.-2./arg1)/a0r                                         
      return                                                            
      end                                                               
C      CCDEF test deck                                                  
c 24. 12. -232. 90. +20. 1 0.3 0.0 0.22 0.12                              
c 110. 130. 2.                                                            
c 2 0                                                                     
c 0.3 -2 -1.                                                              
c 0.3  2 -2.
