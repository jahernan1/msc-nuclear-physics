from numpy import arange, array, exp, sqrt, empty, pi, interp
from pylab import plot,xlabel,ylabel,show,legend,figure,xlim,suptitle,subplot
from scipy.integrate import simps
from scipy.interpolate import splrep,splev
from fractions import Fraction

#----------------------------------MODULES-----------------------------------

def f(r,x,en): 
    
    g = r[0]
    f = r[1]

    if nucleon == 'proton':
        v = splev( x,tck_PV ) / m    # ->interpolates value at x using
                                     #   cubic spline for the given potential
        s = splev( x,tck_PS ) / m    
    elif nucleon == 'neutron':
        v = splev( x,tck_NV ) / m
        s = splev( x,tck_NS ) / m
        
    if x == 0: x = 1.0e-5

    dg = ( en - v + s + 1 )*f - ( kappa/x )*g
    df = ( kappa/x )*f - ( en - v - s - 1 )*g

    return array([dg,df],float)

#------------------------------------------------------------------------

def intsol(a,b,h,en):

    if kappa < 0 : r = array([0.0001,0.0],float)
    else : r = array([0.0,0.0001],float)

    for x_int in arange(a,b,h):
        gint.append(r[0])
        fint.append(r[1])

        k1 = h*f( r,x_int,en )
        k2 = h*f( r+0.5*k1,x_int+0.5*h,en )
        k3 = h*f( r+0.5*k2,x_int+0.5*h,en )
        k4 = h*f( r+k3,x_int+h,en )
        r += ( k1 + 2*( k2+k3 ) + k4 )/6

#-------------------------------------------------------------------

def extsol(c,d,h,en):

    if kappa < 0 : r = array([0.0001,0.0],float)
    else : r = array([0.0,0.0001],float)
    
    for x_ext in arange(d,c,h):
        gext.append(r[0])
        fext.append(r[1])

        k1 = h*f( r,x_ext,en )
        k2 = h*f( r+0.5*k1,x_ext+0.5*h,en )
        k3 = h*f( r+0.5*k2,x_ext+0.5*h,en )
        k4 = h*f( r+k3,x_ext+h,en )
        r += ( k1 + 2*( k2+k3 ) + k4 )/6

#---------------------------------------------------------------------

def det(en):
    intsol(a,b,h,en)
    extsol(c,d,-h,en)

    g1 = gint[-1]
    g2 = gext[-1]
    f1 = fint[-1]
    f2 = fext[-1]

    return (g2*f1) - (g1*f2)

#----------------------------------------------------------------------

def switch(e1,e2):
        
    if det(e1) > 0 and det(e2) > 0: return True
    elif det(e1) < 0 and det(e2) < 0: return True
    else: return False 

#---------------------------------------------------------------------------

def _occupied(n,lst):

    sortedlst = sorted(lst)
    newlist = []
    a = 0

    for elm in sortedlst:
        if a >= n: break
        
        a += ( 2 * elm[1] + 1 )    # elm[1] is j-value
        newlist.append( elm )

    return newlist

#----------------------------------------------------------------------------

def _density( n,lst,typ ):   # typ is 'vector' or 'scalar'
    a = 0
    density = 0
    
    for elm in lst:
        g_wf = array( elm[2][0] )
        f_wf = array( elm[2][1] )

        n_occ = ( 2 * elm[1] + 1 )   # number of total occupiable nucleons
        a += n_occ

        # if not 'magic' takes only the nucleons occupying that level
        
        if a > n : n_occ = n - ( a - n_occ)  

        x = xaxis
        x[0] = 1.0e-12
        
        if typ == 'vector':
            density += ( ( n_occ / ( 4*pi*(x**2) ) )) * ( g_wf**2 + f_wf**2)

        elif typ == 'scalar':
            density += ( ( n_occ / ( 4*pi*(x**2) ) )) * ( g_wf**2 - f_wf**2)
            
    return density

#--------------------------------------CODE-----------------------------------


m = 939.0     # avg poton and neutron rest mass [MeV]
kappa_values = array([-4,-3,-2,-1,1,2],float)
#en = 0.0

h = 0.04*(m/197.326)      # dimensionless interval; rint = 0.04 

a = 0.0                     # b = 1/m
b = 25.0 + h                    # r = b*x -> x = r/b = r*b
c = 25.0
d = 20.0*(m/197.326)      # r = 20.0 fm

xaxis = arange(a,d+h,h)

# All potential are in [MeV] 

neutron_Sf = open( 'Potentials/Ca40_SN.dat','r' )
neutron_S = array( neutron_Sf.read().split(),float )
neutron_Sf.close()

neutron_Vf = open( 'Potentials/Ca40_VN.dat','r' )
neutron_V = array( neutron_Vf.read().split(),float )
neutron_Vf.close()

proton_Sf = open( 'Potentials/Ca40_SP.dat','r' )
proton_S = array( proton_Sf.read().split(),float )
proton_Sf.close()

proton_Vf = open( 'Potentials/Ca40_VP.dat','r' )
proton_V = array( proton_Vf.read().split(),float )
proton_Vf.close()

tck_PS = splrep( xaxis,proton_S )       # cubic spline of the potentials 
tck_PV = splrep( xaxis,proton_V )
tck_NS = splrep( xaxis,neutron_S )
tck_NV = splrep( xaxis,neutron_V )

nucleon = 'proton'
wavefunctions = []

print '{}\n'.format(nucleon)
    
for k in kappa_values:
    kappa = k
    en = 0.94
    
    if kappa > 0: l = kappa
    else: l = -1 - kappa
    j = abs(kappa) - 0.5
    
    print 'kappa = {}; j = {}; l = {} \n'.format( kappa,Fraction(j),l )
    
    while (en < 1.0): 
        target = 1.0e-6
        s = 0.02
        en2 = en + s

        print ( en, en2, s )
        
        while abs(en2-en) > target:
            gint = []
            gext = []    
            fint = []
            fext = []

            if switch(en,en2) == True:
                if en > 1.0: break
                en += s
                en2 = en + s
            
            else:
                s = -s/2
                en =  en2
                en2 = en + s

            print en, en2, s
            print det(en), det(en2)
            
        if en > 1.0:
            print 'No more bound states.'
            break

        print ('Bound State Epsilon: {}'.format(en))  
        print ('Bound State Energy: {} MeV'.format(en*m))
        print ('Binding Energy (E-m): {} MeV \n'.format(en*m-m))
        
        gint = []
        gext = []    
        fint = []
        fext = []

        intsol(a,b,h,en)
        extsol(c,d,-h,en)

        coef = gint[-1] / gext[-1]
                    
        gext1 = [coef * i for i in gext[::-1]]
        del gext1[0]
        fext1 = [coef * i for i in fext[::-1]]
        del fext1[0]    

        gwavefunc = array(gint + gext1)
        fwavefunc = array(fint + fext1)

        norm_const = sqrt(1 / (simps( ( gwavefunc**2+fwavefunc**2 ) ) ))
        
        gnorm_wavefunc = [norm_const*i for i in gwavefunc]       
        fnorm_wavefunc = [norm_const*i for i in fwavefunc]

        wavefunctions.append([en,j,[gnorm_wavefunc,fnorm_wavefunc]])
        
        figure()
        plot( arange(a,d+h,h), gnorm_wavefunc,'-b',label='g(x)' )
        plot( arange(a,d+h,h), fnorm_wavefunc,'-r',label='f(x)' )
        
        suptitle( 'Dirac Wavefunctions (k = {})'.format(kappa)  ) 
        legend( loc= 'upper right' )
        xlabel( 'x' )
            

vector_density = array( _density( 20, _occupied( 20,wavefunctions ), 'vector' ))
print 'Number of {}: {}'.format(nucleon,simps( 4*pi*(xaxis**2)*vector_density ))
print 'Number of protons: {}'.format(simps( 4*pi*(xaxis**2)*vector_density ))
            
scalar_density = array( _density( 20, _occupied( 20,wavefunctions ), 'scalar' ))

figure()

plot( xaxis, vector_density, '-b',label='vector_p' )
plot( xaxis, scalar_density, '-r',label='scalar_p' )
legend( loc='upper right' )
xlabel( 'x' )
        
show()
