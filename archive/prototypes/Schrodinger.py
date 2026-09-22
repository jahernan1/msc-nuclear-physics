from numpy import arange, array, exp, sqrt
from pylab import plot,xlabel,ylabel,show,legend,figure,xlim,suptitle,subplot
from scipy.integrate import simps

#-----------------------------$ Jesse Hernandez $------------------------------#
'''


#---------------------------------$ MODULES $----------------------------------#
                      

'''
#------------------------------------------------------------------------------#
                      # -> Isotropic Harmonic Oscillator Potential  
def iho_V(x):         # 
    return x**2       # -> b**2 = hbar / (m * w) (length scale)
                      # -> en = 2*E / (hbar * w) (dimensionless energy)
#------------------------------------------------------------------------------#
                                       # -> Woods-Saxon potential
def ws_V(x):                           # 
    A = 40                             # -> Units: fm(femptometers, fermis)
    b = sqrt(197.3**2/(2*931.5*50))    # -> b**2 = (2*m)/ (v*hbar**2) (b=0.64)
    cbar = ((1.25)*(A**(1.0/3.0)))/b   # -> cbar = c / b 
    abar = 0.5/b                       #    ( c = 1.25 * A^(1/3)-mass #) 
                                       #
    return -(1/(1+exp((x-cbar)/abar))) # -> abar = a / b ( a = 0.5 )
                                       # 
                                       # -> en = E / v (Unit : MeV, for E and v)
                                       # 
#------------------------------------------------------------------------------#

def finite_V(x):
    #V = 100.0
    #a = 9.23  #corresponds to 4.2 fm'''
    #a = 20.65
    a = 27.7
    
    if x <= a : return -V
    else : return 0

#-----------------------------------------------------------------------------# 

def f(r,x,en):         # -> Decoupled First Order Equations from schrodinger
    u = r[0]           #    equation
    y = r[1]           # 
    fu = y             #----------------------------# 

    if x == 0:  x = 1.0e-12                         # -> Deal with singularity 
      
    fy = ( (finite_V(x) / V ) + ( ( l*( l+1 ) ) /x**2 )-en )*u            # -> Plug in the potential   

    return array([fu,fy],float)                     #    from above modules
                                                    #  
#------------------------------------------------------------------------------#
                                   # 
def intsol(u0,y0,a,b,h,en):        # -> Runge-Kutta Solutions for the interior
                                   #    function
    r = array([u0,y0],float)       # -> Initial conditions   
                                   #  
    for x_int in arange(a,b,h):    # -> Runs over chosen range by step size, h
        upoints.append(r[0])       #
        ypoints.append(r[1])       #
                                   #---------#
        k1 = h*f(r,x_int,en)                 # -> Runge-Kutta Equations
        k2 = h*f(r+0.5*k1,x_int+0.5*h,en)    # 
        k3 = h*f(r+0.5*k2,x_int+0.5*h,en)    #
        k4 = h*f(r+k3,x_int+h,en)            #   
        r += (k1+2*k2+2*k3+k4)/6             #  
                                             #    
#------------------------------------------------------------------------------#
                                   #  
def extsol(u0,y0,a,b,h,en):        # -> Runge-Kutta solutions for the exterior 
    r = array([u0,y0],float)       #    function
                                   #  
    for t_ext in arange(b,a,h):    # -> Runs from back to front
        ext_upoints.append(r[0])   #           
        ext_ypoints.append(r[1])   #---------# NOTE: MUST RUN WITH NEGATIVE 
        k1 = h*f(r,t_ext,en)                 #       STEPSIZE, h
        k2 = h*f(r+0.5*k1,t_ext+0.5*h,en)    #
        k3 = h*f(r+0.5*k2,t_ext+0.5*h,en)    #
        k4 = h*f(r+k3,t_ext+h,en)            #
        r += (k1+2*k2+2*k3+k4)/6             #
                                             #
#------------------------------------------------------------------------------#
                                     #
def det(en,a,b,c,d,h):               # -> Finds the determinate at the matching
    intsol(0.0,0.05,a,b,h,en)        #    point
    extsol(0.0,1.0e-1,c,d,-h,en)     #
                                     #   
    u1 = upoints[-1]                 #  
    u2 = ext_upoints[-1]             # 
    y1 = ypoints[-1]                 #   
    y2 = ext_ypoints[-1]             #
    return -u1 * y2 + u2 * y1        # 
                                     #
#------------------------------------------------------------------------------#
                                                            #    
def switch(e1,e2,a,b,c,d,h):                                # -> Checks if sign
    if det(e1,a,b,c,d,h) > 0 and det(e2,a,b,c,d,h) > 0:     #    of determinate
        return True                                         #    changes
    elif det(e1,a,b,c,d,h) < 0 and det(e2,a,b,c,d,h) < 0:   #
        return True                                         #
    else: return False                                      #
                                                            #
#------------------------------------------------------------------------------#
'''                                                            

#--------------------------------$ MAIN CODE $---------------------------------#

'''
#------------------------------------------------------------------------------#


upoints = []
ypoints = []    
ext_upoints = []
ext_ypoints = []

nwfs = []
#en = -0.02
V = 900
a = 0.0
b = 27.8
c = 27.6
d = 35.0
h = 0.1



for n in range(1):
    l = n
    en = -1.02
    print ('Angular momentum number: {}'.format(l))
    incr = 1
    #for n in range(3): #for iho

    while (en < 0):
        target = 1.0e-12
        s = 0.01
        en += 2*s
        en2 = en + s
       
        while abs(en2-en) > target:
            upoints = []
            ypoints = []    
            ext_upoints = []
            ext_ypoints = []

            if switch(en,en2,a,b,c,d,h) == True:
                if en > 0: break
                en += s
                en2 = en + s
            
            else:
                s = -s/4
                en =  en2
                en2 = en + s
                                     
        if en > 0:
            print 'No more bound states.'
            break
         
        print ('Bound State Energy: {}'.format(en))  
        print ('Bound State Energy: {} MeV'.format(en*V))
        
        upoints = []
        ypoints = []    
        ext_upoints = []
        ext_ypoints = []
    
        intsol(0.0,0.05,a,b,h,en)
        extsol(0.0,-h,c,d,-h,en)
        
        coef = upoints[-1] / ext_upoints[-1]  
        ext_upoints1 = [coef * i for i in ext_upoints[::-1]]
        del ext_upoints1[0]
        ext_ypoints1 = [coef * i for i in ext_ypoints[::-1]]
        del ext_ypoints1[0]

        wavefunc = upoints + ext_upoints1
        sq_wavefunc = [i**2 for i in wavefunc] 
        
        norm = sqrt(1/simps(sq_wavefunc))
        norm_wavefunc = [norm*i for i in wavefunc]
        nwfs.append(norm_wavefunc)

        dwavefunc = ypoints + ext_ypoints1
        norm_dwavefunc = [norm*i for i in dwavefunc]
        
        #print('Normalization Check: {}' "\n".format(simps([i**2 for i in norm_wavefunc])))
       
        figure()
        #subplot(3,1,1+n)
        
        plot(arange(a,d+h,h), norm_wavefunc,'-b',label='u')
        plot(arange(a,d+h,h), norm_dwavefunc,'-r',label='du/dx')

        #xlim(0,10)
        suptitle('Harmonic Oscillator Wavefunction (l = {})'.format(int(l))) 
        legend(loc='upper right')
        ylabel('u(x)')
        xlabel('x')

        x = arange(a, d+h , h)
       
        udata = open( '/home/jesse/Documents/Data_Files/udata_schrod_{0}MeV_{1}.dat'.format( str(V), str(incr) ), 'w+' )
                
        for i in arange(0,len(x)):
            udata.write( '{0:<9}{1:^1} \n'.format( str(x[i]),str(norm_wavefunc[i]) ) )
        incr += 1
show()


#-------------------------------------------------------------------        
# wf01 = simps(array(nwfs[0],float)*array(nwfs[1],float))
# wf02 = simps(array(nwfs[0],float)*array(nwfs[2],float))
# wf12 = simps(array(nwfs[1],float)*array(nwfs[2],float))


# print ('Check Orthagonality (<n|m> = 0):' "\n")
# print ('<0|1> = <1|0> = {}'.format(wf01))
# print ('<0|2> = <2|0> = {}'.format(wf02))

# print ('<1|2> = <2|1> = {}'.format(wf12))


