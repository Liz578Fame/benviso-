
       Program iso_benv_calc()
       implicit real*8(a-h,o-z)
       common/const_iso_control/n_control_, iphen_print_,
     1                          iso_calc_, iso_print_, parab_print_ 
       common/const_nnv_control/nskin_opt_, nskin_print, nden_opt_,
     1                          nden_print_, nff_opt_, nff_print_
                                            
       icmd = 1                               
       call iso_pipeline()       
       write(*,*) icmd
                                   
       n_benv = n_control_     
                                
       if(n_benv .EQ. 1) then         
          call be_const_read()
          call benv_pipeline() 
       end if  
       write(*,*) icmd                 
                                  
       end program                    

c----------------------------------------------
c                 Subroutines                 |
c----------------------------------------------

       subroutine iso_pipeline()
       implicit real*8(a-h,o-z)

       common/const_iso_control/n_control_, iphen_print_,
     1                          iso_calc_, iso_print_, parab_print_
       
       common/const_main/n,mic,isnm,isym_emp,k0,rho0    
       common/const_lvl_1/den_s,den_e,pi,pi2,fact,hbc,hbc2,xm,totfact
       common/const_lvl_2/alpha,beta,sigma,alph,a1,b1,c1,d1,d2,ff1,ff2
       common/eos_gp/xkf(100),den(100),e0(100),e1(100),esym(100)       
       
       common/deng/desym(100),de0(100),de1(100)
       common/pressure/prse0(100),prse1(100),prsesym(100)    
       common/isovals_e/rhoo,e0o,rho1,prs0,prs02,esym0,esym1
       common/isovals_de/bigL,bigK,bigKD,bigK0
       
c       open(000,file='dump.don')   
       icmd = 1    
        
c      const_main: read in constants from stdin
       call input_constants()
       n_control = n_control_
       iphen_print = iphen_print_
       iso_calc = iso_calc_
       iso_print = iso_print_
       parab_print = parab_print_  
        
c      const_lvl 1 & 2
       call eos_constants()
       call eos_read_in()
        
c      compute phenomonological EoS
       call phen_eos()
c      print phenom_EoS
       if(iphen_print .EQ. 1) then 
          call phen_eos_print()
       end if
       
c      Compute Isovector and Isoscalar properties
       if(iso_calc .EQ. 1) then 
          call calc_iso_val()
       end if 
        
c      Print Isovector and Isoscalar Properties
       if(iso_calc .EQ. 1 .AND. iso_print .EQ. 1) then
          call print_iso_val()
       end if
              
c      Print range of parabolic EoS
       if(parab_print .EQ. 1) then 
          call parab()
       end if
                    
       end subroutine


                                 
       subroutine benv_pipeline() 
       implicit real*8(a-h,o-z)
       common/const_nnv_control/nskin_opt_, nskin_print, nden_opt_, 
     1                          nden_print_, nff_opt_, nff_print_
                  
       iradius_step = 20
       radlim = 10.d0   
       nfac_num = 40    
       nfac_step = 5     
         
c      const_main: parse control constants       
       nskin_opt  = nskin_opt_
       nskin_print= nskin_opt_
       nden_opt   = nden_opt_
       nden_print = nden_print_
       nff_opt    = nff_opt_
       nff_print  = nff_print_
        
c      optimize density paramters upon the 
c      minimization of the binding energy        
       call be_server()        
           
       if(nskin_opt .EQ. 1) then
          call skin_calc()
       end if       

       if(nskin_opt .EQ. 1 .AND. nskin_print .EQ. 1) then
          call skin_print()
       end if
            
       if(nden_opt .EQ. 1) then      
          call nucden(iradius_step,radlim)
       end if
                 
       if(nden_print .EQ. 1 .AND. nden_opt .EQ. 1) then
          call nucden_print()
       end if

       if(nff_opt .EQ. 1) then 
          call form_fac_gen(nfac_num,nfac_step)
       end if
                         
       if(nff_print .EQ. 1 .AND. nff_opt .EQ. 1) then
          call form_fac_print(nfac_num)
       end if
            
       end subroutine
            
             
c---------------------------------------------------------------------------------------------------
c                                         Server Routine                                           |
c---------------------------------------------------------------------------------------------------

       subroutine be_server()
       implicit real*8 (a-h,o-z)
       common/azn/ta, tz, tn
       common/setup/n1, n2, n3, x1, x2
       common/pars/rp,cp,wp,rn,cn,wn  

c      SERVER INTIALIZATION
c       call init(n1, n2, n3, x1, x2, ta, tz)

c-------------------------------------------------
c             SERVER MAIN LOOP                   |
c-------------------------------------------------

c      read the command and parameters
c      if command is 0 (or rather "not 1"), evaluate energy and repeat
c      if command is 1, stop the server

       icmd = 0
       icmdpp = 1
       write(*,*) icmdpp
       
       do 100 while (icmd .ne. 1)
          
          read(*,*) icmd,rp,cp,wp,rn,cn,wn
                  
          if (icmd .ne. 1) then
             en = energy(rp, cp, wp, rn, cn, wn)
c         write the energy value to STDOUT
             write(*,*) en              
          end if
 100   continue  
        
       end subroutine

             
c--------------------------------------------------------------------------

c      Main Routines
c      
c      iso_pipeline()
c      
c      
c      Isovalue Subroutines dependencies:  
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
