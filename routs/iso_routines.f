      
c      Subroutines
c      
c      input_constants()
c      eos_constants()
c      read_in()
c      phen_eos()
c      calc_iso_val()
c      print_iso_val()
c      parab()
c      phen_eos_print()
         
c--------------------------------------------------------------------------
        
       subroutine input_constants()
       implicit real*8(a-h,o-z) 
       common/const_main/n,mic,isnm,isym_emp,k0,rho0       
       common/const_iso_control/n_control_, iphen_print_, 
     1                          iso_calc_, iso_print_, parab_print_
           
c       open(unit=000, file='dump.don')
       icmd = 1
       write(*,*) icmd              
       read(*,*) n_control_, iphen_print_, iso_calc_, iso_print_,
     1           parab_print_
       write(*,*) icmd
       read(*,*) n, mic, isnm, isym_emp, k0, rho0
        
       end subroutine
          
c--------------------------------------------------------------------------

       subroutine eos_constants()
       implicit real*8(a-h,o-z)
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/const_lvl_2/alpha,beta,sigma,alph,a1,b1,c1,d1,d2,ff1,ff2

       den_s = 0.0d0
       den_e = 1.8d0       
       pi = 3.14159265
       pi2 = pi*pi 

       fact=(3.d0*pi2/2.d0)**(2.d0/3.d0) 
       hbc=197.327d0
       hbc2=hbc**2
       xm=938.926d0 
       tfact=(3.d0*hbc2/10.d0/xm)
       totfact=fact*tfact
c
c      k = 220,260
c
       alpha=-29.47-46.74*(k0+44.21)/(k0-166.11)
       beta=23.37*(k0+254.53)/(k0-166.11)
       sigma=(k0+44.21)/210.32
       gam=0.72d0
       alph=0.2d0 
       a1=119.14d0
       b1=-816.95d0
       c1=724.51d0
       d1=-32.99d0
       d2=891.15d0
       ff1=a1*2.d0*(0.5d0)**(5.d0/3.d0)
       ff2=d1*2.d0*(0.5d0)**(5.d0/3.d0) + 
     1     d2*2.d0*(0.5d0)**(8.d0/3.d0) 

       end subroutine

c--------------------------------------------------------------------------

       subroutine eos_read_in()
       implicit real*8(a-h,o-z)
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)       

       open(unit=510, file='ex.don')

       do i = 1,n
          read(510,*) xkf(i), e0(i), e1(i)
          den(i) = (2.d0*xkf(i)**3)/(3.d0*pi2)
          esym(i) = e1(i) - e0(i)
       end do
       return
       end 

c-------------------------------------------------------------------------- 
      
       subroutine phen_eos()
       implicit real*8(a-h,o-z)

       common/const_main/n,mic,isnm,isym_emp,k0,rho0 
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)   
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/const_lvl_2/alpha,beta,sigma,alph,a1,b1,c1,d1,d2,ff1,ff2

       dimension :: ee(100), ee2(100)    

       do i=1,n 
          datapt = den(i)
          rat=datapt/rho0
          ee(i)=ff1*(datapt)**(2.d0/3.d0) + b1*datapt +
     1    c1*(datapt)**(alph+1.d0) + ff2*(datapt)**(5.d0/3.d0) 
                       
          ee2(i)=totfact*(datapt)**(2.d0/3.d0) + (alpha/2.d0)*(rat)+
     1    (beta/(sigma + 1.d0))*(rat)**(sigma) 
           
          
c         mic: microscopic EoS         
c         mic = 1 : EoS is completely microscopic
c         mic = 0 : EoS options for phenomenological EoS

          if(mic.eq.1) go to 3355 
          
c         isnm = 1: symmetric matter EoS is phenomenological choice 1
c         isnm = 0: symmetric matter EoS is phenomenological choice 2
 
             if(isnm.eq.1) then
                e0(i)=ee(i) 
             else 
                e0(i)=ee2(i) 
             end if 
          
c             if(datapt.le.0.0019.and.e0.gt.0.d0) then
c                e0(i)=0.d0 
c             end if
               
             pt=22.d0*rat**gam + 12.d0*rat**(2.d0/3.d0)      

c            isym_emp = 0: esym is partially phenomenological
c                          the neutron matter part is microcopic
c            isym_emp = 1: esym is completely phenomenological 

             if(isym_emp.eq.0) then
                esym(i)= e1(i)-e0(i) 
             else 
                esym(i)=pt      
             end if      
             go to 3355 
          
3355      continue
       end do

       end subroutine

c--------------------------------------------------------------------------

       subroutine calc_iso_val()
       implicit real*8(a-h,o-z)   
       common/const_main/n,mic,isnm,isym_emp,k0,rho0 
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)  
       common/deng/desym(100),de0(100),de1(100)
       common/pressure/prse0(100),prse1(100),prsesym(100)    
       common/isovals_e/rhoo,e0o,rho1,prs0,prs02,esym0,esym1
       common/isovals_de/bigL,bigK,bigKD,bigK0
                
       dimension der(100)
       dimension outpt1(100), outpt2(100), coeff1(4,100), coeff2(4,100)
       dimension outpt3(100), outpt4(100), coeff3(4,100), coeff4(4,100)
       dimension P(100),  rho(100), rhog(100)
       dimension e0d(100), e1d(100)
       dimension xkf0(100), xkf1(100), den0(100), den1(100)
       dimension de0dr(100), ee(100), ee2(100)
       dimension outpt5(100), coeff5(100)
        
       n2 = n-1       
        
       call dcsakm(n,den,esym,outpt1,coeff1)
       call dcsakm(n,den,e0,outpt2,coeff2)
       call dcsakm(n,den,e1,outpt3,coeff3)
       
       do i=1,n
         desym(i) = dcsder(1,den(i),n2,outpt1,coeff1)
         de0(i) = dcsder(1,den(i),n2,outpt2,coeff2)
         de1(i) = dcsder(1,den(i),n2,outpt3,coeff3)
         prse0(i) = den(i)*den(i)*de0(i)   
         prse1(i) = den(i)*den(i)*de1(i)
         prsesym(i) = den(i)*den(i)*desym(i) 
       end do
         
       call dcsakm(n,de0,den,outpt4,coeff4)
       call dcsakm(n,den,prse0,outpt5,coeff5)
          
       rhoo = dcsval(0.d0,n2,outpt4,coeff4)
       rhoo_2 = rhoo*2.d0
       e0o = dcsval(rhoo,n2,outpt2,coeff2)
       rho1 = 0.1d0
       prs0 = dcsval(rhoo,n2,outpt5,coeff5)
       prs02 = dcsval(rhoo_2,n2,outpt5,coeff5)
         
       esym0 = dcsval(rhoo,n2,outpt1,coeff1)
       esym1 = dcsval(rho1,n2,outpt1,coeff1)
          
       bigL = 3.d0*rhoo*dcsder(1,rhoo,n2,outpt1,coeff1) 
       bigK = 9.d0*rhoo*rhoo*dcsder(2,rhoo,n2,outpt1,coeff1)
       bigKD = 9.d0*rho1*rho1*dcsder(2,rho1,n2,outpt1,coeff1)
       bigK0 = 9.d0*rhoo*rhoo*dcsder(2,rhoo,n2,outpt2,coeff2)
        
       end subroutine

c--------------------------------------------------------------------------
       
       subroutine print_iso_val()
       implicit real*8(a-h,o-z)
       common/const_main/n,mic,isnm,isym_emp,k0,rho0 
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)   
       common/deng/desym(100),de0(100),de1(100)
       common/pressure/prse0(100),prse1(100),prsesym(100)    
       common/isovals_e/rhoo,e0o,rho1,prs0,prs02,esym0,esym1
       common/isovals_de/bigL,bigK,bigKD,bigK0

       open(unit=140, file='data.srt')         
       open(unit=150, file='eos.srt')
       
1400   format(3x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3,
     1        2x,F7.3,2x,F7.3,2x,F7.3,2x,F7.3)
1500   format(2x,F8.4,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3,2x,F8.3
     1        ,2x,F8.3,2x,F8.3)
       
       
1001   format('-----------------------------------------------|')
1000   format('                                         ')         
6734   format('----------------------------------Pressures--------------
     1---------------------|')
4545   format('-------L-Value-------|')
4646   format('-------K-Value-------|')
3636   format('------K-0-Value------|')
3737   format('------K-D-Value------|')
4747   format('-----E-Sat-Value-----|')
4848   format('-----R-Sat-Value-----|')
4949   format('-----P-Sat-Value-----|')

5050   format('----Esym-Sat-Value---|')
5151   format('----Esym-Den-Value---|')


8629   format('     Den      E0       E1',       
     1       '       Esym      Kf      P0       P1       Psym')   

0001   format('-------------------------------Isoscalar-Values----------
     1---------------------|')

1111   format('------------------------------Isovector-Values-----------
     1---------------------|')

2222   format('-----------------------------------Tabulated-Values------
     1-----------------------------|')
        
       if(n_opt .eq. 0) then
c          write(140,6734)
          write(140,1000)
          write(140,2222)
          write(140,1000)
          write(140,*) "   rho_o     E_o      K_o        L_o    ",
     1                 "   K         K_d       E_sym_o   E_sym_den "
          write(140,1500) rhoo, e0o, bigK0, bigL, bigK, bigKD,
     1                    esym0, esym1

          write(140,1000)
          write(150,8629)
          do i=1,n
             write(150,1400) den(i),e0(i),e1(i),esym(i),
     1                       xkf(i),prse0(i),prse1(i),prsesym(i) 
          end do


          write(140,1000)
          write(140,0001)
          write(140,1000)

          write(140,1000)
          write(140,4848)
          write(140,1000)
          write(140,*) rhoo
          write(140,1000)
          write(140,4747)
          write(140,1000)
          write(140,*) e0o
          write(140,1000)
          write(140,4949) 
          write(140,1000) 
          write(140,*) prs0
          write(140,1000)
          write(140,3636)
          write(140,1000)
          write(140,*) bigK0    
          write(140,1000)      

          write(140,1000)
          write(140,1111)
          write(140,1000)
       
          write(140,1000)        
          write(140,4545)
          write(140,1000)
          write(140,*) bigL       
          write(140,1000)
          write(140,4646)
          write(140,1000)
          write(140,*) bigK
          write(140,1000)
          write(140,3737)
          write(140,1000)
          write(140,*) bigKD
          write(140,1000)
          write(140,5050)
          write(140,1000)
          write(140,*) esym0
          write(140,1000)
          write(140,5151)
          write(140,1000)
          write(140,*) esym1
          write(140,1000) 
          
       end if       
       end subroutine

c--------------------------------------------------------------------------

       subroutine parab()
       implicit real*8(a-h,o-z)

       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)
        
       dimension :: ea(6,100), alp(6)

       open(777,file='parab.don')
       
       do i=1,6 
          alp(i) = 0.d0 + 0.2d0*(i-1)
       end do

       do j=1,6
          do i=1,n
             ea(j,i) = e0(i)+esym(i)*alp(j)**2 
          end do
       end do

       write(777,*) '   den      ea0       ea2       ea4       ea6
     1ea8       ea1'
       do i=1,n
          write(777,1414) den(i), ea(1,i), ea(2,i),
     1                    ea(3,i), ea(4,i), ea(5,i), ea(6,i)

       end do

1414   format(2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4,
     1        2x,F8.4,2x,F8.4)

       end subroutine

c--------------------------------------------------------------------------

       subroutine phen_eos_print()
       implicit real*8(a-h,o-z)
       common/const_main/n,mic,isnm,isym_emp,k0,rho0
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)
       
       open(999,file='e0_phens.don')
               
       pi = 3.14159265
       pi2 = pi*pi 
       
c       rho0 = 0.16d0      

c      parameatrization of empirical EoS for symmmetric nuclear matter
        
       xk0 = 220
       xk1 = 260
       
       fact=(3.d0*pi2/2.d0)**(2.d0/3.d0) 
       hbc=197.327d0
       hbc2=hbc**2
       xm=938.926d0 
       tfact=(3.d0*hbc2/10.d0/xm)
       totfact=fact*tfact
c      
       alpha=-29.47-46.74*(xk0+44.21)/(xk0-166.11)
       beta=23.37*(xk0+254.53)/(xk0-166.11)
       sigma=(xk0+44.21)/210.32
       
       alpha1=-29.47-46.74*(xk1+44.21)/(xk1-166.11)
       beta1=23.37*(xk1+254.53)/(xk1-166.11)
       sigma1=(xk1+44.21)/210.32
       
        
       gam=0.72d0
       alph=0.2d0 
       a1=119.14d0
       b1=-816.95d0
       c1=724.51d0
       d1=-32.99d0
       d2=891.15d0
       ff1=a1*2.d0*(0.5d0)**(5.d0/3.d0)
       ff2=d1*2.d0*(0.5d0)**(5.d0/3.d0) + 
     1     d2*2.d0*(0.5d0)**(8.d0/3.d0) 
       
c      phenom_eos_section
       write(999,*) '   kf       den     e0_1      e0_220    e0_260
     1esym'  
       do i=1,n 
          datapt = den(i)
          rat=datapt/rho0
          ee=ff1*(datapt)**(2.d0/3.d0) + b1*datapt +
     1    c1*(datapt)**(alph+1.d0)+ff2*(datapt)**(5.d0/3.d0)  
          ee2=totfact*(datapt)**(2.d0/3.d0)+(alpha/2.d0)*(rat)+
     1    (beta/(sigma+1.d0))*(rat)**(sigma)                
          ee3=totfact*(datapt)**(2.d0/3.d0)+(alpha1/2.d0)*(rat)+
     1    (beta1/(sigma1+1.d0))*(rat)**(sigma1)
          pt=22.d0*rat**gam+12.d0*rat**(2.d0/3.d0)                
          write(999,1010) xkf(i), den(i), ee, ee2, ee3, pt
       end do      
1010   format(2x,F7.3,2x,F7.3,2x,F8.4,2x,F8.4,2x,F8.4,2x,F8.4)
       end
       
