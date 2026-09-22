from numpy import arange, array, exp, sqrt, empty, pi
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
    b = sqrt(197.3**2/(2*931.5*60))    # -> b**2 = (2*m)/ (v*hbar**2) (b=0.64)
    cbar = ((1.25)*(A**(1.0/3.0)))/b   # -> cbar = c / b 
    abar = 0.5/b                       #    ( c = 1.25 * A^(1/3)-mass #) 
                                       #
    return -(1.0/(1.0+exp((x-cbar)/abar))) # -> abar = a / b ( a = 0.5 )
                                       # 
                                       # -> en = E / v (Unit : MeV, for E and v)
                                       # 
#------------------------------------------------------------------------------#
                       # 
def f(r,x,en):         # -> Decoupled First Order Equations from schrodinger
    u = r[0]           #    equation
    y = r[1]           # 
    fu = y             #----------------------------# 
    if x == 0:                                      # -> Deal with singularity 
        fy = (ws_V(x)+((l*(l+1))/1.0e-20)-en)*u     # 
    else: fy = (ws_V(x)+((l*(l+1))/x**2)-en)*u      # -> Plug in the potential  
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

def _occupied(n,lst):

    lst = sorted(lst)
    newlist = []
    a = 0
    for elm in lst:
        if a >= n: return newlist
        a += 2 * ( 2 * elm[1][1] + 1 )
        newlist.append( elm )
        
def _density(n,lst):

    x = arange( 0, 20 + h, h )
    x[0] = 1.0e-12
    a = 0
    p = 0
    
    for elm in lst:
        u_wf = array( elm[2] )
        R_wf2 = (u_wf/x)**2
        a += 2 * ( 2 * elm[1][1] + 1 )
        n_occ =  2 * ( 2 * elm[1][1] + 1 )
              
        if a > n : n_occ = n - ( a - n_occ) 
        
        p += R_wf2 * n_occ

    return p / (4 * pi)   

'''

#--------------------------------$ MAIN CODE $---------------------------------#

'''

#------------------------------------------------------------------------------#


upoints = []
ypoints = []    
ext_upoints = []
ext_ypoints = []

sys_info = []
#en = -1.02
#l = 4.0
a = 0.0
b = 2.05
c = 1.95
d = 20.0
h = 0.05
runs = 0

for i in range(5):
    l = i
    node = 1
    en = -1.02
    print ('Angular momentum number: {}'"\n".format(l))
    
    #for n in range(3): #for iho

    while (en < 0): 
        target = 1.0e-10
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
         
        print ('Bound State Epsilon: {}'.format(en))  
        print ('Bound State Energy: {} MeV'.format(en*60))
        
        upoints = []
        ypoints = []    
        ext_upoints = []
        ext_ypoints = []
    
        intsol(0.0,h,a,b,h,en)
        extsol(0.0,-1.0e-1,c,d,-h,en)
        
        coef = upoints[-1] / ext_upoints[-1]  
        ext_upoints1 = [coef * i for i in ext_upoints[::-1]]
        del ext_upoints1[0]
        ext_ypoints1 = [coef * i for i in ext_ypoints[::-1]]
        del ext_ypoints1[0]

        wavefunc = upoints + ext_upoints1
        sq_wavefunc = [i**2 for i in wavefunc] 
        
        norm = sqrt(1/simps(sq_wavefunc))
        norm_wavefunc = [norm*i for i in wavefunc]       

        dwavefunc = ypoints + ext_ypoints1
        norm_dwavefunc = [norm*i for i in dwavefunc]
        
        print('Normalization Check: {}' "\n".format(simps([i**2 for i in norm_wavefunc])))

        sys_info.append([en,[node,l],norm_wavefunc])
        #print sys_info[runs][1]
        runs += 1
        node += 1
        
        figure()
        #subplot( 3,1,1+n )
        plot( arange(a,d+h,h), norm_wavefunc,'-b',label='u' )
        plot( arange(a,d+h,h), norm_dwavefunc,'-r',label='du/dx' )

        #xlim(0,10)
        suptitle( 'Woods-Saxon Wavefunction (l = {})'.format(int(l)) ) 
        legend( loc= 'upper right' )
        ylabel( 'u(x)' )
        xlabel( 'x' )

figure()
proton_density = _density( 20, _occupied( 20, sys_info ) )
#neutron_density = _density( 126, _occupied( 126, sys_info ) )
plot ( arange( a, d+h, h ), proton_density, '-b', label = 'proton' )  
#plot ( arange( a, d+h, h ), neutron_density, '-r', label = 'neutron' )
legend ( loc = 'upper right' ) 

x = arange ( 0, d + h, h )
x[0] = 1.0e-12
print simps( proton_density*4*pi*(x**2) )
#print simps( neutron_density*4*pi*(x**2) )

# wf01 = simps(array(nwfs[0],float)*array(nwfs[1],float))
# wf02 = simps(array(nwfs[0],float)*array(nwfs[2],float))
# wf12 = simps(array(nwfs[1],float)*array(nwfs[2],float))

# print ('Check Orthagonality (<n|m> = 0):' "\n")
# print ('<0|1> = <1|0> = {}'.format(wf01))
# print ('<0|2> = <2|0> = {}'.format(wf02))
# print ('<1|2> = <2|1> = {}'.format(wf12))

show()
