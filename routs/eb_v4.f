c   Progran eb_v4.f  subroutines: Version 1.0 
c   Written by: Randy Millerson


c---------------------------------------------------------------------------------------------------
c                                   BE Const. Read-In Routine                                      |
c---------------------------------------------------------------------------------------------------

       subroutine azn_update()
       implicit real*8(a-h,o-z) 
       common/azn/ta, tz, tn        
       common/setup/n1, n2, n3, x1, x2
       common/parz/n_den,fff    
       
       icmd = 1
         
       write(*,*) icmd
       read(*,*) n1,n2,n3,x1,x2,ta,tz,n_den,fff
        
       end subroutine 


       subroutine be_const_read()
       implicit real*8 (a-h,o-z)

       common/const_nnv_control/nskin_opt_, nskin_print, nden_opt_, 
     1                          nden_print_, nff_opt_, nff_print_
        
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2
       common/parz/n_den,fff

       icmd = 1

       write(*,*) icmd
       read(*,*) nskin_opt_,nskin_print,nden_opt_,
     1           nden_print_,nff_opt_,nff_print_
       write(*,*) icmd                          
       read(*,*) n1,n2,n3,x1,x2,ta,tz,n_den,fff

       call init(n1, n2, n3, x1, x2, ta, tz)
                      
       end subroutine
                     
c---------------------------------------------------------------------------------------------------
c                                      Initalization Routine                                       |
c---------------------------------------------------------------------------------------------------

       subroutine init(n1_, n2_, n3_, x1_, x2_, ta_, tz_) 
       implicit real*8 (a-h,o-z)                     
       
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact

       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)

       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/binding/totbe,binde1,binde2,binde3
       common/maineos/xdata(75),zdata(75),nxdata,
     3                 xsnm(75), ydata(75),nsnm, 
     1                 breakz(75),cscoefz(4,75),       
     2                 breaky(75),cscoefy(4,75),       
     4                 xdatas(100),xdatan(100),
     5                 sdata(75),breaks(75),cscoefs(75)

       common/main/fint1,fint2,fint3       
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2

c       open(000,file='dump.don')
 
c    number of integration points                             
       n1 = n1_
       n2 = n2_
       n3 = n3_
c     integration limits 
       x1 = x1_
       x2 = x2_
c     characteristics of the nucleus 
       ta = ta_
       tz = tz_
       tn=ta-tz 

c        number of points for the Nuclear Matter EoS
c       n=11
        
       nxdata = n
       nxnm = n 
              
       do i=1,n
          ydata(i) = e0(i)
          zdata(i) = e1(i)
          sdata(i) = esym(i)
          xdata(i) = den(i)
c          write(*,*) xdata(i), ydata(i), zdata(i), sdata(i)
       end do
                
c      set up interpolation
         
       call dcsakm(n,xdata,zdata,breakz,cscoefz)
       call dcsakm(n,xdata,ydata,breaky,cscoefy)
       call dcsakm(n,xdata,sdata,breaks,cscoefs)       
       return
       end
        
c---------------------------------------------------------------------------------------------------
c                            Total Energy Per Particle Calculuations                               |
c---------------------------------------------------------------------------------------------------
c     Main Function: Calculates the energy per particle for given rho() parameters      
       function energy(rp, cp, wp, rn, cn, wn)
       implicit real*8 (a-h,o-z)                     
       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/binding/totbe,binde1,binde2,binde3
       common/abc/xnorm  
       common/main/fint1,fint2,fint3         
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2
       
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       
c       open(000,file='dump.don')

       dt = 0.d0

c     normalize the proton function
       if(n_den .EQ. 2 .OR. n_den .EQ. 4) then
          call xnormalize(pi,rp,cp,tz,dt)
       else if(n_den .EQ. 3) then
          call xnormalize(pi,rp,cp,tz,wp)
       end if
       ap=xnorm
         
              
c     normalize the neutron function
       if(n_den .EQ. 2 .OR. n_den .EQ. 4) then
          call xnormalize(pi,rn,cn,tn,dt)
       else if(n_den .EQ. 3) then
          call xnormalize(pi,rn,cn,tn,wn)
       end if
       an=xnorm 
         
          
       if(n_den .EQ. 2 .OR. n_den .EQ. 4) then
          call eos(n1,pi,x1,x2,ap,rp,cp,dt,an,rn,cn,dt)
          call be2(n2,pi,x1,x2,ap,rp,cp,dt,an,rn,cn,dt)
          call be3(n3,pi,x1,x2,ap,rp,cp,dt)
       else if(n_den .EQ. 3) then                      
          call eos(n1,pi,x1,x2,ap,rp,cp,wp,an,rn,cn,wn)               
          call be2(n2,pi,x1,x2,ap,rp,cp,wp,an,rn,cn,wn)               
          call be3(n3,pi,x1,x2,ap,rp,cp,wp)      
       end if    

c       write(000,*) ap, an, binde1, binde2, binde3
       energy = (binde1+binde2+binde3)/ta
                       
       end 

c---------------------------------------------------------------------------------------------------
c                                 Normalization Routines                                           |
c---------------------------------------------------------------------------------------------------
       subroutine xnormalize(pi,b,c,xg,wi) 
c 
       implicit real*8 (a-h,o-z)                     
       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/abc/xnorm 
       common/parz/n_den,fff

c       open(000,file='dump.don')
       
       sum=0.d0
       x1=0.d0
       x2=20.d0 
       xinf=0.d0
       nnorm=90    
       xnorm=1.d0 
       call lgauss(nnorm)
       call papoi(x1,x2,1,nnorm,xinf,1)
       sum=0.d0 
       do 10 i=1,nnorm
          r=pas(i)
          ww=poi(i)
          if(n_den .EQ. 2) then
             funct=rho(r,xnorm,b,c)*r**2*ww
          else if(n_den .EQ. 3) then 
             funct=rho_3pf(r,xnorm,b,c,wi)*r**2*ww
          else if(n_den .EQ. 4) then
             funct=rho_fy(r,xnorm,b,c,xg)*r**2*ww
          end if    
10     sum=sum+funct
       sum=sum*4.d0*pi
       xnorm=xg/sum
       return
       end 


c---------------------------------------------------------------------------------------------------


c---------------------------------------------------------------------------------------------------
c                                         Density Functions                                        |
c---------------------------------------------------------------------------------------------------


       function rho(xr,xa,xb,xc)
       implicit real*8 (a-h,o-z)                     
          rho=xa/(1.d0+dexp((xr-xb)/xc))
       end 


       function rho_3pf(xr,xa,xb,xc,xw)
       implicit real*8 (a-h,o-z)
          rho_3pf=(1.d0+((xr**2)*dabs(xw))/(xb**2))*
     1        (xa/(1.d0+dexp((xr-xb)/xc)))
       end 


       function rho_fy(xr,xa,xb,xc,xg)
       implicit real*8(a-h,o-z)
       cof = (3.d0*xg)/(4.d0*(3.14159d0)*xb**3)
       erterm = (xc/xr)*dexp(-(xr/xc))
       ebrterm = (xc/xr)*dexp(-(xb/xc))
       sh = dsinh(xb/xc)
       ch = dcosh(xb/xc)
       if(xr .lt. xb) then
          rho_fy = cof*( 1.d0 - ebrterm*(1.d0+(xb/xc))*dsinh(xr/xc))
       else if(xr .gt. xb) then
          rho_fy = cof*erterm*((xb/xc)*ch - sh)
       end if
       end   

c-------------------------------------------------
c               Density Derivative               |
c-------------------------------------------------

       function devrho(xr,xa,xb,xc)
       implicit real*8 (a-h,o-z)       
          devrho=(-(xa/xc)*dexp((xr-xb)/xc))/
     1 ((1.d0+dexp((xr-xb)/xc))**2)    
       end 


       function devrho_3pf(xr,xa,xb,xc,xw)
       implicit real*8 (a-h,o-z)
          devrho_3pf = devrho(xr,xa,xb,xc) + 
     1    2.d0*(dabs(xw)/xb**2)*xr*rho(xr,xa,xb,xc)+
     1    (dabs(xw)/xb**2)*(xr**2)*devrho(xr,xa,xb,xc)
       end


       function devrho_fy(xr,xa,xb,xc,xg)
       implicit real*8(a-h,o-z)
       cof = (3.d0*xg)/(4.d0*(3.14159d0)*xb**3)
       erterm = (xc/xr)*dexp(-(xr/xc))
       ebrterm = (xc/xr)*dexp(-(xb/xc))
       sh = dsinh(xb/xc)
       ch = dcosh(xb/xc)
       shr = dsinh(xr/xc)
       chr = dcosh(xr/xc)
       if(xr .lt. xb) then
          devrho_fy = -cof*(1.d0+(xb/xc))*dexp(-xb/xc)*
     1                ((chr/xr)-(xc*shr/xr**2))
       else if(xr .gt. xb) then
          devrho_fy = -cof*((xb/xc)*ch - sh)*
     1                ((erterm/xr)+(erterm/xc))
       end if
       end


c---------------------------------------------------------------------------------------------------
c                                        SEMF Components                                           |
c---------------------------------------------------------------------------------------------------

c-------------------------------------------------
c           Volume Nuclear Energy Term           |  
c-------------------------------------------------

c      This subroutine calculates the energy as a function of 
c      density and the symmetric energy parameter.

       subroutine eos(m,pi,x1,x2,a,b,c,wt,a2,b2,c2,wt2)             

       implicit real*8 (a-h,o-z)                     

       common/paspoi/pas(200),poi(200),x(200),w(200)
       common/binding/totbe,binde1,binde2,binde3
       common/maineos/xdata(75),zdata(75),nxdata,
     3                 xsnm(75), ydata(75),nsnm, 
     1                 breakz(75),cscoefz(4,75),       
     2                 breaky(75),cscoefy(4,75),  
     4                 xdatas(100),xdatan(100),
     5                 sdata(75),breaks(75),cscoefs(75)
       common/azn/ta, tz, tn
       dimension ee(100),ee2(100) 
       dimension xx(100)            

       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff
        
c  Generates the interpolated neutron matter EoS
c  parametrization of empirical EoS for symmmetric nuclear matter

       nintv=nxdata-1
       nintv2=nsnm-1 
c
       call lgauss(m)
       call papoi(x1,x2,1,m,0.d0,1)       
 
       sum=0.d0
       sum2=0.d0
       do 20 i=1,m
          r=pas(i)
          ww=poi(i)
          if(n_den .EQ. 2) then
             datapt=rho(r,a,b,c)+rho(r,a2,b2,c2)                      
             alp=(rho(r,a2,b2,c2)-rho(r,a,b,c))/datapt
          else if(n_den .EQ. 3) then
             datapt=rho_3pf(r,a,b,c,wt)+rho_3pf(r,a2,b2,c2,wt2)
             alp=(rho_3pf(r,a2,b2,c2,wt2)-rho_3pf(r,a,b,c,wt))/datapt
          else if(n_den .EQ. 4) then
             datapt=rho_fy(r,a,b,c,tz)+rho_fy(r,a2,b2,c2,tn)
             alp=(rho_fy(r,a2,b2,c2,tn)-rho_fy(r,a,b,c,tz))/datapt
          end if
          alp2=alp*alp

          e0_inc = dcsval(datapt,nintv,breaky,cscoefy)
          e1_inc = dcsval(datapt,nintv,breakz,cscoefz)
          esym_inc = dcsval(datapt,nintv,breaks,cscoefs)          

          if(isym_emp.eq.0) then
              pt=e0_inc+alp2*(e1_inc-e0_inc) 
          else
              pt=e0_inc+alp2*esym_inc
          end if
                       
          pt=pt*datapt*r**2.d0*ww
          sum=sum+pt       
 20    continue 
        
       finalint=sum*4.d0*pi
       binde1=finalint
       return
       end 


c---------------------------------------------------------------------------------------------------

c-------------------------------------------------
c           Surface Nuclear Energy Term          |  
c-------------------------------------------------


       subroutine be2(nn,pi,x1,x2,a,b,c,wt,a2,b2,c2,wt2)  

       implicit real*8 (a-h,o-z)                     
       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/binding/totbe,binde1,binde2,binde3
       common/azn/ta, tz, tn

       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff

       dimension f1(200),f2(200) 
       dimension break1(200),cscoef1(4,200)
       dimension break2(200),cscoef2(4,200)
       
       beta=1.d0     
c       fff = 65.d0
      
       n2=nn-1 
       call lgauss(nn)
       call papoi(x1,x2,1,nn,0.d0,1)   
       
       do i=1,nn
          r=pas(i)
          if(n_den .eq. 2) then
             f1(i)=rho(r,a,b,c)
          else if(n_den .eq. 3) then
             f1(i)=rho_3pf(r,a,b,c,wt)
          else if(n_den .eq. 4) then
             f1(i)=rho_fy(r,a,b,c,tz)
          end if
          if(n_den .eq. 2) then
             f2(i)=rho(r,a2,b2,c2) 
          else if(n_den .eq. 3) then
             f2(i)=rho_3pf(r,a2,b2,c2,wt2)
          else if(n_den .eq. 4) then
             f2(i)=rho_fy(r,a2,b2,c2,tn)
          end if     
       end do
 
       sum=0.d0
       do 20 i=1,nn 
          r=pas(i)
          ww=poi(i)
          if(n_den .EQ. 2) then
             funct1=devrho(r,a,b,c)
             funct2=devrho(r,a2,b2,c2)
          else if(n_den .EQ. 3) then
             funct1=devrho_3pf(r,a,b,c,wt)
             funct2=devrho_3pf(r,a2,b2,c2,wt2) 
          else if(n_den .EQ. 4) then
             funct1=devrho_fy(r,a,b,c,tz)
             funct2=devrho_fy(r,a2,b2,c2,tn)
          end if
          funct=funct1+funct2 
          funct3=(funct1-funct2)**2 
          funct=funct**2 
          funct=funct*fff*r**2*ww                                  
          sum=sum+funct
  20   continue 
       
       binde2=sum*4.d0*pi
       return
       end 
       
       
c---------------------------------------------------------------------------------------------------

c-------------------------------------------------
c              Coulomb Energy Term               |
c-------------------------------------------------


       subroutine be3(nn,pi,x1,x2,a,b,c,wt)  
       implicit real*8 (a-h,o-z)                     
       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/binding/totbe,binde1,binde2,binde3
       common/azn/ta, tz, tn

       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff

       dimension rp(200),wrp(200)       
c      
       n1=nn
       n2=nn
       call lgauss(nn)
       call papoi(x1,x2,1,nn,0.d0,1)        
       
       do 1 i=1,nn
          rp(i)=pas(i)
          wrp(i)=poi(i)
 1     continue 
       
       call lgauss(nn)
       sum=0.d0
       
       do 20 j=1,nn  
          xlow=0.d0
          xup=rp(j)
          call papoi(xlow,xup,1,nn,0.d0,1)
          sum2=0.d0
          do 10 i=1,nn
             r=pas(i)
             ww=poi(i)
c            internal integral
             if(n_den .EQ. 2) then
                funct=rho(r,a,b,c)*r**2.d0*ww
             else if(n_den .EQ. 3) then
                funct=rho_3pf(r,a,b,c,wt)*r**2.d0*ww
             else if(n_den .EQ. 4) then
                funct=rho_fy(r,a,b,c,tz)*r**2.d0*ww
             end if
             sum2=sum2+funct
 10       continue 
c     external integral
          if(n_den .EQ. 2) then
             functex=rp(j)*rho(rp(j),a,b,c)*sum2*wrp(j)
          else if(n_den .EQ. 3) then
             functex=rp(j)*rho_3pf(rp(j),a,b,c,wt)*sum2*wrp(j)
          else if(n_den .EQ. 4) then
             functex=rp(j)*rho_fy(rp(j),a,b,c,tz)*sum2*wrp(j)
          end if
          sum=sum+functex
 20    continue
       
       binde3=sum*(4.d0*pi)**2*1.44d0
       return
       end
c---------------------------------------------------------------------------------------------------

c---------------------------------------------------------------------------------------------------
c                            Atomic Parameters Calculations                                        |  
c---------------------------------------------------------------------------------------------------

c-------------------------------------------------
c         Charge Radius (Root Mean Squared)      |
c-------------------------------------------------

       subroutine chrms(a1,b1,c1,w1,pi,ttz,n)
       implicit real*8(a-h,o-z)
       common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
       common/charge/chr 
       common/azn/ta, tz, tn

       common/parz/n_den,fff

       a=0.70 
c      a=0.87*dsqrt(2.d0/3.d0)
       x1=0.d0
       x2=20.d0 
       call lgauss(n)
       call papoi(x1,x2,1,n,xinf,1)

       sum2=0.d0
       do 106 j=1,n
          r=pas(j)
          ww=poi(j)
          fact3=dexp(-(r/a)**2.d0)
          sum1=0.d0
          do 105 i=1,n
             r1=pas(i)
             ww1=poi(i)
             fact1=dexp(-(r1/a)**2.d0)
             fact2=dsinh(2.d0*r*r1/(a**2.d0))
             if(fact2.gt.1.0d300) then
                funct1=0.d0
             else
                if(n_den .EQ. 2) then
                   funct1=r1*ww1*fact1*rho(r1,a1,b1,c1)*
     1             fact2*fact3/r
                else if(n_den .EQ. 3) then
                   funct1=r1*ww1*fact1*rho_3pf(r1,a1,b1,c1,w1)*
     1             fact2*fact3/r
                else if(n_den .EQ. 4) then
                   funct1=r1*ww1*fact1*rho_fy(r1,a1,b1,c1,tz)*
     1             fact2*fact3/r
                end if
             end if 
             sum1=sum1+funct1
  105     continue
          dinte=2.d0*sum1/(a*dsqrt(pi))
          funct2=4.d0*pi*r**4.d0*ww*dinte
          sum2=sum2+funct2
 106   continue 
       chr=dsqrt(sum2/ttz)
       end 


c-------------------------------------------------
c         Root Mean Squared Routine              |
c-------------------------------------------------

       subroutine rms(a1,b1,c1,w1,a2,b2,c2,w2,pi,tz,tn,ta,n)
       implicit real*8(a-h,o-z)
       common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
       common/main/fint1,fint2,fint3
       common/parz/n_den,fff
  
       sum1=0.d0
       sum2=0.d0
       sum3=0.d0
       sum4=0.d0 
       x1=0.d0
       x2=20.d0
 
       call lgauss(n)
       call papoi(x1,x2,1,n,xinf,1)
c
       do 10 i=1,n
          r=pas(i)
          ww=poi(i)
          if(n_den .EQ. 2) then
             funct1=rho(r,a1,b1,c1)
             funct2=rho(r,a2,b2,c2)
             funct3=rho(r,a1,b1,c1) + rho(r,a2,b2,c2)
          else if(n_den .EQ. 3) then
             funct1=rho_3pf(r,a1,b1,c1,w1)
             funct2=rho_3pf(r,a2,b2,c2,w2)
             funct3=rho_3pf(r,a1,b1,c1,w1) + rho_3pf(r,a2,b2,c2,w2)
          else if(n_den .EQ. 4) then
             funct1=rho_fy(r,a1,b1,c1,tz)
             funct2=rho_fy(r,a2,b2,c2,tn)
             funct3=rho_fy(r,a1,b1,c1,tz) + rho_fy(r,a2,b2,c2,tn)
          end if
          funct1=funct1*r**4.d0
          funct2=funct2*r**4.d0
          funct3=funct3*r**4.d0
          sum1=sum1+funct1*ww
          sum2=sum2+funct2*ww
          sum3=sum3+funct3*ww
  10   continue

       const1=(4.d0*pi)/tz
       const2=(4.d0*pi)/tn
       const3=(4.d0*pi)/ta
       fint1=dsqrt(sum1*const1)
       fint2=dsqrt(sum2*const2)
       fint3=dsqrt(sum3*const3)
       end 

c---------------------------------------------------------------------------------------------------
c                                        Skin radii routines                                       |
c---------------------------------------------------------------------------------------------------

       subroutine skin_calc()
       implicit real*8(a-h,o-z)

       common/main/fint1,fint2,fint3
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2       

       common/pars/rp,cp,wp,rn,cn,wn
       
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)
       common/deng/desym(100),de0(100),de1(100)
       
       common/paspoi/pas(200),poi(200),x(200),w(200) 
       common/binding/totbe,binde1,binde2,binde3
       common/maineos/xdata(75),zdata(75),nxdata,
     3                 xsnm(75), ydata(75),nsnm, 
     1                 breakz(75),cscoefz(4,75),       
     2                 breaky(75),cscoefy(4,75),       
     4                 xdatas(100),xdatan(100),
     5                 sdata(75),breaks(75),cscoefs(75)
       common/abc/xnorm 
       common/charge/chr  

       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/parz/n_den,fff
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/const_lvl_2/alpha,beta,sigma,alph,a1,b1,c1,d1,d2,ff1,ff2

       common/skins/en,rad_prot,rad_neut,rad_chrg,
     1              skin_neut,skin_prot,esym_coef,den_ref

       dimension :: frac(100), coef(100,4) 

       nintv=nxdata-1

       en = energy(rp, cp, wp, rn, cn, wn)

c     normalize the proton function
       call xnormalize(pi,rp,cp,tz,wp)
       ap=xnorm       

       call xnormalize(pi,rn,cn,tn,wn)
       an=xnorm
       
       call rms(ap,rp,cp,wp,an,rn,cn,wn,pi,tz,tn,ta,64)
       call chrms(ap,rp,cp,wp,pi,tz,64)
         
       rad_prot = fint1
       rad_neut = fint2
       rad_chrg = chr
       skin_neut = fint2-fint1
       skin_prot = fint1-fint2
                      
       sum=0.d0
       down_lim = 0.0d0
       upper_lim = 20.d0
       nrep = 90
 
       if(tn .ne. tz) then
          coef_int = ta/(tn-tz)**2
       else 
          coef_int = 0.d0
       end if         

       call lgauss(nrep)
       call papoi(down_lim,upper_lim,1,nrep,xinf,1)
         
       do i =1,nrep
          dr = pas(i)
          ww = poi(i)
          if(n_den .EQ. 2) then
             rho_t = rho(dr,an,rn,cn) + rho(dr,ap,rp,cp)
             delt = (rho(dr,an,rn,cn) - rho(dr,ap,rp,cp))/rho_t
          else if(n_den .EQ. 3) then
             rho_t = rho_3pf(dr,an,rn,cn,wn) + rho_3pf(dr,ap,rp,cp,wp)
             delt = (rho_3pf(dr,an,rn,cn,wn) - 
     1               rho_3pf(dr,ap,rp,cp,wp))/rho_t
          else if(n_den .EQ. 4) then
             rho_t = rho_fy(dr,an,rn,cn,tn) + rho_fy(dr,ap,rp,cp,tz)
             delt = (rho_fy(dr,an,rn,cn,tn) -
     1               rho_fy(dr,ap,rp,cp,tz))/rho_t
          end if
           
          e0_eval = dcsval(datapt,nintv,breaky,cscoefy)
          e1_eval = dcsval(datapt,nintv,breakz,cscoefz)
          esym_eval = dcsval(datapt,nintv,breaks,cscoefs)           
          
          if(isym_emp.eq.0) then
             pt=(delt*delt)*(e1_eval-e0_eval) 
          else 
             pt=(delt*delt)*esym_eval      
          end if        
          
          val = rho_t*pt*(dr*dr)*ww
          sum = sum + val
       end do
       esym_coef = sum*coef_int*4.d0*pi
                
       call dcsakm(n,esym,den,frac,coef)
       den_ref = dcsval(esym_coef,nintv,frac,coef)          
       return 
       end subroutine
            
            
            
       subroutine skin_print()
       implicit real*8(a-h,o-z)
       common/skins/en,rad_prot,rad_neut,rad_chrg,
     1              skin_neut,skin_prot,esym_coef,den_ref
            
       open(945,file='skin_vals.don')
          
       fnr = rad_neut
       fpr = rad_prot
       fns = skin_neut
       frc = rad_chrg
       esc = esym_coef
       drf = den_ref
          
       write(945,*) ' BE        RN       RP       NS       CR       SC
     1     DR'
       write(945,7755) en, fnr, fpr, fns, frc, esc, drf
               
7755   format(2x,F7.4,2x,F7.4,2x,F7.4,
     1        2x,F7.4,2x,F7.4,2x,F7.4,2x,F7.4)
         
       end subroutine 
               
c---------------------------------------------------------------------------------------------------
c                                    Nuclear Densities Routine                                     |
c---------------------------------------------------------------------------------------------------

       subroutine nucden(nstep,aplim)
       implicit real*8(a-h,o-z)
       common/main/fint1,fint2,fint3
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2
         
       common/pars/rp,cp,wp,rn,cn,wn
         
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)
       common/deng/desym(100),de0(100),de1(100)
       
       common/paspoi/pas(200),poi(200),x(200),w(200)
       common/binding/totbe,binde1,binde2,binde3
       common/maineos/xdata(75),zdata(75),nxdata,
     3                 xsnm(75), ydata(75),nsnm,
     1                 breakz(75),cscoefz(4,75),
     2                 breaky(75),cscoefy(4,75),
     4                 xdatas(100),xdatan(100),
     5                 sdata(75),breaks(75),cscoefs(75)
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/abc/xnorm
       common/parz/n_den,fff
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       
       common/den_list/nstep_val,rad_step(100),
     1                 rhop(100),rhon(100),rhot(100)
                
       nstep_val = nstep

       call xnormalize(pi,rp,cp,tz,wp)
       ap=xnorm
       call xnormalize(pi,rn,cn,tn,wn)
       an=xnorm
           
       do i=1,nstep
          rad_step(i) = aplim*(1.d0*(i-1)/(1.d0*(nstep-1)))
       end do       
              
       if(n_den .EQ. 2) then
          do i=1,nstep
             rhop(i) = rho(rad_step(i),ap,rp,cp)
             rhon(i) = rho(rad_step(i),an,rn,cn)
             rhot(i) = rhop(i) + rhon(i)
          end do
       else if(n_den .EQ. 4) then
          do i=1,nstep
             rhop(i) = rho_fy(rad_step(i),ap,rp,cp,tz)
             rhon(i) = rho_fy(rad_step(i),an,rn,cn,tn)
             rhot(i) = rhop(i) + rhon(i)
          end do
       else if(n_den .EQ. 3) then
          do i=1,nstep
             rhop(i) = rho_3pf(rad_step(i),ap,rp,cp,wp)
             rhon(i) = rho_3pf(rad_step(i),an,rn,cn,wn)
             rhot(i) = rhop(i) + rhon(i)
          end do
       end if 
       return 
       end subroutine         
       


       subroutine nucden_print()
       implicit real*8(a-h,o-z)
       common/den_list/nstep_val,rad_step(100),
     1                 rhop(100),rhon(100),rhot(100)
          
       nn = nstep_val
       open(909,file='nuc_dens.don') 
       write(909,2332) 
       do i=1,nn
          write(909,2300) rad_step(i), rhop(i), rhon(i), rhot(i) 
       end do 
         
2332   format(2x,"dens",5x"rhop",5x,"rhon"5x,"rhot")
1000   format("                     ")
2300   format(F7.3,2x,F7.3,2x,F7.3,2x,F7.3)       
          
       end subroutine 
        
c---------------------------------------------------------------------------------------------------
c                                       Form Factor routine                                        |
c---------------------------------------------------------------------------------------------------

       subroutine form_fac(qf,ac,rc,dc,wc,pi,tz,n)
       implicit real*8(a-h,o-z) 
       common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
       common/form_pars/ff,ff2,ff2_log
       common/parz/n_den,fff

       xinf = 0.d0
       x1 = 0.d0
       x2 = 15.d0
         
       call lgauss(n)
       call papoi(x1,x2,1,n,xinf,1)
          
       coef = 4.d0*pi/tz
          
       sum=0.d0
       do j=1,n
          r=pas(j)
          ww=poi(j)
          if(dabs(qf*r) .LT. 0.0000001) then
             xjo = 1.0d0
          else 
             xjo = dsin(qf*r)/(qf*r) 
          end if
      
          if(n_den .EQ. 2) then
             rhozi = rho(r,ac,rc,dc)
          else if(n_den .EQ. 3) then
             rhozi = rho_3pf(r,ac,rc,dc,wc)
          else if(n_den .EQ. 4) then
             rhozi = rho_fy(r,ac,rc,dc,tz)
          end if
          gran = rhozi*xjo*(r*r)*ww 
          sum = sum + gran
       end do
       ff = coef*sum       
       ff2 = ff*ff
       ff2_log = dlog10(ff2)
       return 
       end subroutine

       subroutine form_fac_gen(nrep,nup)
       implicit real*8(a-h,o-z) 
       common/form_pars/ff,ff2,ff2_log
       common/pars/rp,cp,wp,rn,cn,wn
       common/azn/ta, tz, tn
       common/abc/xnorm
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/form_gen/fkf_gp(500),ff_gp(500),ff2_gp(500),ff2log_gp(500)

       call xnormalize(pi,rp,cp,tz,wp)
       ap = xnorm

       xinc = (1.d0*nup)/((nrep-1)*1.d0)

       x = 0.d0
       do i=1,nrep
          call form_fac(x,ap,rp,cp,wp,pi,tz,90)
          fkf_gp(i) = x
          ff_gp(i) = ff
          ff2_gp(i) = ff2
          ff2log_gp(i) = ff2_log
          x = x + xinc
       end do        
       end subroutine

c---------------------------------------------------------------------
c                     Form Factor Print routine                      |
c---------------------------------------------------------------------

       subroutine form_fac_print(nrep)
       implicit real*8(a-h,o-z)  
       common/form_gen/fkf_gp(500),ff_gp(500),ff2_gp(500),ff2log_gp(500)
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
         
       open(919,file='form_fac.don')  
       write(919,*)'  kf               ff               ff2             
     1  fflog'
               
       do i=1,nrep  
          write(919,1010) fkf_gp(i), ff_gp(i), ff2_gp(i), ff2log_gp(i)
       end do

1010   format(F15.10,2x,F15.10,2x,F15.10,2x,F15.10)       
       end subroutine 

c---------------------------------------------------------------------------------------------------
c                            Numeric Integration Subroutines                                       |
c---------------------------------------------------------------------------------------------------

c-------------------------------------------------
c          Legendre-Gaussian Integration         |
c-------------------------------------------------
       subroutine lgauss(n)
       implicit real*8(a-h,o-z)
       dimension z(200),wz(200)
       common/paspoi/pas(200),poi(200),xfs(200),wfs(200)
       
       if(n-1) 1,2,3
 1        return
 2        z(1)=0.d0
          wz(1)=2.d0
          return
 3        r=dfloat(n)
       g=-1.d0
       
       do 147 i=1,n
          test=-2.d0
          ic=n+1-i
          if(ic.lt.i) go to 150
 4           s=g
             t=1.d0
             u=1.d0
             v=0.d0
          do 50 k=2,n
             a=dfloat(k)
             p=((2.d0*a-1.d0)*s*g-(a-1.d0)*t)/a
             dp=((2.d0*a-1.d0)*(s+g*u)-(a-1.d0)*v)/a
             v=u
             u=dp
             t=s
 50       s=p
          if(abs((test-g)/(test+g)).lt.0.5d-09) go to 100
             sum=0.d0
          if(i.eq.1) go to 52
             do 51 j=2,i
                sum=sum+1.d0/(g-xfs(j-1))
 51          continue
 52       test=g
          g=g-p/(dp-p*sum)
          go to 4
 100      xfs(ic)=-g
          xfs(i)=g
          wfs(i)=2.d0/(r*t*dp)
          wfs(ic)=wfs(i)
 147   g=g-r*t/((r+2.d0)*g*dp+r*v-2.d0*r*t*sum)
 150   do 160 i=1,n
          z(i)=xfs(i)
          wz(i)=wfs(i)
 160   continue
       return
       end
c-------------------------------------------------
c              Steps and Weights                 |
c-------------------------------------------------
       subroutine papoi(xi,xf,ni,nf ,xinf,k)
       implicit real*8(a-h,o-z)
       common/paspoi/pas(200),poi(200),tp(200),hp(200)

       coeff=1.d0           
       go to (10,20,30) k
 10       xs=(xi+xf)/2.d0
          xd=(xf-xi)/2.d0
          do 1 i=ni,nf
             pas(i)=xs+xd*tp(i)
 1           poi(i)=xd*hp(i)
          return
 20       do 2 i=ni,nf
             pas(i)=((1.d0+tp(i))/(1.d0-tp(i)))*coeff+xinf
             poi(i)=2.d0*hp(i)/((1.d0-tp(i))**2)
             poi(i)=poi(i)*coeff
c 2   write(15,*) pas(i),poi(i),i
  2       continue
          return
 30    do 3 i=ni,nf
          pas(i)=xi*(1.d0+tp(i))/(1.d0-tp(i))+xinf
          poi(i)=2.d0*xi*hp(i)/((1.d0-tp(i))**2)
 3     continue                       
       return
       end
c-------------------------------------------------

