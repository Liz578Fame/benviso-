import sys                                   #System Module
import os                                    #Operating Module
import io                                    #Pipeline Module
import subprocess                            #Subprocess Module
import time                                  #Timing Module
from py_libs import PMOD_26 as Parse         #Parsing Module
  
import scipy.optimize as sciopt              #import SciPy module with optimization routines
import shutil


class fort_opt_server():
     
    def __init__(self,n_opt_var):
        self.n_opt_var = n_opt_var    
        
    def com_server_orig(self,serv_name,write):    

        global n_opt
        global function_call_counter
         
        n_opt = self.n_opt_var
        
        server = subprocess.Popen(serv_name, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.STDOUT)
        function_call_counter = 0   #counts the number of energy() evaluations
        # definition of function to minimize
        def energy(v): # v is a list of parameters. Modified Woods-Saxon Distribution v = [rp, cp, wp, rn, cn, wn]
            global n_opt
            global function_call_counter
            #write the "EVALUATE ENERGY" command (icmd = 0)  and parameter values into the STDIN of the FORTRAN server
            if(n_opt != 3):
                server.stdin.writelines(["0\n" + str(v[0]) + "\n" + str(v[1]) + "\n" + str(0.0) 
                                            + "\n" + str(v[2]) + "\n" + str(v[3]) + "\n" + str(0.0) + "\n"])
            if(n_opt == 3):
                    server.stdin.writelines(["0\n" + str(v[0]) + "\n" + str(v[1]) + "\n" + str(v[2]) + "\n" 
                                            + str(v[3]) + "\n" + str(v[4]) + "\n" + str(v[5]) + "\n"])        
            #wait for FORTRAN server responce and read the energy value from from its STDOUT 
            en = float(server.stdout.readline())
             
            function_call_counter += 1  #increase the counter

            #print out parameters and energy every n iterations
            if True: #set to False to disable printing
                n = 10
                if n == 1 or (function_call_counter % n) == 1:
                    print v, en	
            return en    #return the value of the energy
        
        if(n_opt != 3):
            v = [7.0,0.3,7.0,0.3] #set up the initial guess for parameters
        if(n_opt == 3):
            v = [7.0,0.3,3.0,7.0,0.3,3.0]
        
        time_0 = time.time() #start timing
         
        res = sciopt.minimize(
                          fun = energy,  #target function to be minimized. energy() in our case
                          x0 = v,        #initial guess for solution
                          method = 'Nelder-Mead',
                          tol = 1e-6   #max relative error between consecutive iterations
                         )
        
        print res    #And the result is:...

        time_1 = time.time() #timing
        server.communicate("1\n0\n0\n0\n0\n0\n0\n") # issue a "STOP" command (icmd = 1) to the FORTRAN server, wait for it to finish

        nma = res.x

        if(write):
         
            outputfile = open("opt_par.etr", "w")
                 
            if(n_opt != 3):
                outputfile.writelines([str(nma[0]),"  ",str(nma[1]),"  ",str(0.0),"  ",str(nma[2]),"  ",str(nma[3]),"  ",str(0.0),"\n"])
            if(n_opt == 3):
                outputfile.writelines([str(nma[0]),"  ",str(nma[1]),"  ",str(nma[2]),"  ",str(nma[3]),"  ",str(nma[4]),"  ",str(nma[5]),"\n"])
            outputfile.close()

    def os_eb_opt_serv(self,server):    

        global n_opt
        global function_call_counter
        global s_server
          
        s_server = server                 
        n_opt = self.n_opt_var
         
        function_call_counter = 0   #counts the number of energy() evaluations
        # definition of function to minimize
        def energy(v): # v is a list of parameters. Modified Woods-Saxon Distribution v = [rp, cp, wp, rn, cn, wn]
            global n_opt
            global function_call_counter
            #write the "EVALUATE ENERGY" command (icmd = 0)  and parameter values into the STDIN of the FORTRAN server
                            
            if(n_opt != 3):
                out_str_sum = "0\n"+str(v[0])+"\n"+str(v[1])+"\n"+str(0.0)+"\n"
                out_str_sum = out_str_sum+"\n"+str(v[2])+"\n"+str(v[3])+"\n"+str(0.0)+"\n"
                s_server.stdin.writelines([out_str_sum])
            if(n_opt == 3):
                out_str_sum = "0\n"+str(v[0])+"\n"+str(v[1])+"\n"+str(v[2])
                out_str_sum = out_str_sum+"\n"+str(v[3])+"\n"+str(v[4])+"\n"+str(v[5])+"\n"
                s_server.stdin.writelines([out_str_sum])     
   
            #wait for FORTRAN server responce and read the energy value from from its STDOUT 
            en = float(server.stdout.readline())
             
            function_call_counter += 1  #increase the counter

            #print out parameters and energy every n iterations
            if True: #set to False to disable printing
                n = 10
                if n == 1 or (function_call_counter % n) == 1:
                    print v, en	
            return en    #return the value of the energy
        
        if(n_opt != 3):
            v = [7.0,0.3,7.0,0.3] #set up the initial guess for parameters
        if(n_opt == 3):
            v = [7.0,0.3,3.0,7.0,0.3,3.0]
        
        time_0 = time.time() #start timing
         
        res = sciopt.minimize(
                          fun = energy,  #target function to be minimized. energy() in our case
                          x0 = v,        #initial guess for solution
                          method = 'Nelder-Mead',
                          tol = 1e-6   #max relative error between consecutive iterations
                         )
        
#        print res    #And the result is:...

        time_1 = time.time() #timing        
        nma = res.x
         
        if(n_opt != 3):           
            s_server.stdin.writelines(["1\n"+str(nma[0])+"\n"+str(nma[1])+"\n"+str(0.0)+"\n"+str(nma[2])+"\n"+str(nma[3])+"\n"+str(0.0)+"\n"])
        elif(n_opt == 3):
            s_server.stdin.writelines(["1\n"+str(nma[0])+"\n"+str(nma[1])+"\n"+str(nma[2])+"\n"+str(nma[3])+"\n"+str(nma[4])+"\n"+str(nma[5])+"\n"])
         
#        server.communicate("1\n0\n0\n0\n0\n0\n0\n") # issue a "STOP" command (icmd = 1) to the FORTRAN server, wait for it to finish    

def par_read_in():

    blnk = ' '
    nl   = '\n'
    nul  = '0\n'
    par_iso = Parse.list_file_grab('par.don',[],False,True)

    control_seq = par_iso[1]
    control_iso = par_iso[4]
    control_nnv = par_iso[7]
    control_azn = par_iso[10]

    # Read-in Iso options
    n_control = control_seq[0]       #
    iphen_print = control_seq[1]     #
    iso_calc = control_seq[2]        #
    iso_print = control_seq[3]       #
    parab_print = control_seq[4]     #

    # Read-in Iso values
    n = control_iso[0]               #
    mic = control_iso[1]             #
    isnm = control_iso[2]            #
    isym_emp = control_iso[3]        #
    k0 = control_iso[4]              #
    rho0 = control_iso[5]            #
    
    # Read-in NNV options
    nskin_opt  = control_nnv[0]      #
    nskin_print= control_nnv[1]      #
    nden_opt   = control_nnv[2]      #
    nden_print = control_nnv[3]      #
    nff_opt    = control_nnv[4]      #
    nff_print  = control_nnv[5]      #

    isop = [iso_print,iphen_print,parab_print]
    nnvp = [nskin_print,nden_print,nff_print]      

    colp = [isop,nnvp]
       
    # Read-in NNV values      
    n1    = control_azn[0]
    n2    = control_azn[1]
    n3    = control_azn[2]
    x1    = control_azn[3]
    x2    = control_azn[4]
    ta    = control_azn[5]
    tz    = control_azn[6]
    nden  = control_azn[7]
    fff   = control_azn[8]
          
    control_seq_str = nl.join(control_seq)+nl
    control_iso_str = nl.join(control_iso)+nl
    control_nnv_str = nl.join(control_nnv)+nl
    control_azn_str = nl.join(control_azn)+nl    
          
    return control_seq_str,control_iso_str,control_nnv_str,control_azn_str,n_control,nden,colp
             


def iso_pipeline(server,con_seq_str,con_iso_str):
     
    continue_bool = 1
    continue_bool = server.stdout.readline()
    server.stdin.writelines([con_seq_str])
    continue_bool = server.stdout.readline()
    server.stdin.writelines([con_iso_str])
    continue_bool = server.stdout.readline()
    
    
def nnv_pipeline(server,con_nnv_str,con_azn_str):

    continue_bool = server.stdout.readline()
    server.stdin.writelines([con_nnv_str])
    continue_bool = server.stdout.readline()
    server.stdin.writelines([con_azn_str])

    continue_bool = server.stdout.readline()
    fos = fort_opt_server(nden_val)
    fos.os_eb_opt_serv(server)

    continue_bool = server.stdout.readline()


def collect(colp):
    isop = colp[0]
    nnvp = colp[1]

    cd = os.getcwd()
    print(cd)
    
    if(int(isop[0]) == 1):
       shutil.move('eos.srt',"/data")
       shutil.move('data.srt',"/data") 
    if(int(isop[1]) == 1): 
       shutil.move('eo_phens.don',"/data")  
    if(int(isop[2]) == 1): 
       shutil.move('parab.don',"/data")
    
    if(int(nnvp[0]) == 1): 
       shutil.move('skin_vals.don',"/data")
    if(int(nnvp[0]) == 1):
       shutil.move('nuc_dens.don',"/data")
    if(int(nnvp[0]) == 1):
       shutil.move('form_fac.don',"/data")


# MAIN

# server set-up 
serv_str = './server'
server = subprocess.Popen(serv_str, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.STDOUT)

#Read in Parameters
con_seq_str,con_iso_str,con_nnv_str,con_azn_str,benv_opt,nden_val,colp = par_read_in()
benv_opt = int(benv_opt)


iso_pipeline(server,con_seq_str,con_iso_str)
if(benv_opt == 1):
    nnv_pipeline(server,con_nnv_str,con_azn_str)

collect(colp)


'''
continue_bool = 1
continue_bool = server.stdout.readline()                      
server.stdin.writelines([con_seq_str])
continue_bool = server.stdout.readline()
server.stdin.writelines([con_iso_str])
continue_bool = server.stdout.readline()

if(int(benv_opt) == 1):
     
    continue_bool = server.stdout.readline()     
    server.stdin.writelines([con_nnv_str])  
    continue_bool = server.stdout.readline()        
    server.stdin.writelines([con_azn_str])  

    continue_bool = server.stdout.readline()            
    fos = fort_opt_server(nden_val)             
    fos.os_eb_opt_serv(server)

    continue_bool = server.stdout.readline()
'''
          
       
