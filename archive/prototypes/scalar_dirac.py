from numpy import arange, array, exp, sqrt, empty, pi
from pylab import plot,xlabel,ylabel,show,legend,figure,xlim,suptitle,subplot
from scipy.integrate import simps
from fractions import Fraction
from matplotlib import rc

#rc('font',**{'family':'sans-serif','san-serif':['Helvetica']})
#rc('text', usetex=True)


def f(r,x,en):

    g = r[0]
    f = r[1]
    s = finite_S(x)/m

    if x == 0: x = 1.0e-12

    dg = ( en + 1 + s )*f - ( kappa/x )*g
    df = ( kappa/x )*f - ( en - 1 - s )*g

    return array([dg,df],float)

#------------------------------------------------------------------------

def finite_S(x):
    a = 20.0
    
    if x <= a : return -S
    else : return 0
        
#-----------------------------------------------------------------------
 
#------------------------------------------------------------------

def intsol(a,b,h,en):
    
    if kappa < 0 : r = array([.0001,0.0],float)
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

def det(a,b,c,d,h,en):

    intsol(a,b,h,en)
    extsol(c,d,-h,en)

    g1 = gint[-1]
    g2 = gext[-1]
    f1 = fint[-1]
    f2 = fext[-1]

    return (g2*f1) - (g1*f2)

#----------------------------------------------------------------------

def switch(a,b,c,d,h,e1,e2):
    if det(a,b,c,d,h,e1) > 0 and det(a,b,c,d,h,e2) > 0:
        return True                                        
    elif det(a,b,c,d,h,e1) < 0 and det(a,b,c,d,h,e2) < 0:  
        return True
    else: return False 


def normalize( g,f ):
    g = array(g,float)
    f = array(f, float)
    
    return sqrt(1 / (simps((g**2+f**2))))
#-------------------------------------------------------------------------------

kappa = -1.0
m = 939.0     # rest mass [MeV]
potentials = [100.0,500.0,900.0]
#en = 0.0

a = 0.0
b = 20.05
c = 19.95
d = 30.0
h = 0.05

for s in potentials:
    
    en = 0.0
    S = s
    plot_incr = 1
    figure()
    
    while (en < 1.0): 
        target = 1.0e-12
        step = 0.01
        en2 = en + step
    
        while abs(en2-en) > target:
            gint = []
            gext = []    
            fint = []
            fext = []
            
            if switch(a,b,c,d,h,en,en2) == True:
                if en > 1.0: break
                en += step
                en2 = en + step
                      
            else:
                step = -step/4
                en =  en2
                en2 = en + step
            
        if en > 1.0:
           print 'No more bound states.'
           break

        print ('Bound State Epsilon: {}'.format(en))  
        print ('Bound State Energy: {} MeV'.format(en*m))
        print ('Binding Energy: {} MeV' '\n'.format(en*m-m))

        gint = []
        gext = []    
        fint = []
        fext = []
        
        intsol(a,b,h,en)
        extsol(c,d,-h,en)
        
        coef = gint[-1] / gext[-1]
        
        gext = [coef * i for i in gext[::-1]]
        del gext[0]
        fext = [coef * i for i in fext[::-1]]
        del fext[0]

        gwavefunc = gint + gext
        fwavefunc = fint + fext
        
        norm_const = normalize( gwavefunc, fwavefunc )
        
        # gsq_wavefunc = [i**2 for i in gwavefunc]
        # fsq_wavefunc = [i**2 for i in fwavefunc]
    
        gnorm_wavefunc = [norm_const*i for i in gwavefunc]       
        fnorm_wavefunc = [norm_const*i for i in fwavefunc]
    
        #figure()
       
        subplot( 3,1,plot_incr  )
                             
        plot( arange(a,d+h,h)*(197.327/939.0), gnorm_wavefunc,'-b',label='g(x)' )
        plot( arange(a,d+h,h)*(197.327/939.0), fnorm_wavefunc,'-r',label='f(x)' )
       
        suptitle(' Finite Dirac Wavefunctions kappa = {} '.format(kappa) ) 
        legend( loc= 'upper right' )
        xlabel( 'r (fm)' )

        x = arange(a,d+h,h)
        
        gdata = open( '/home/jesse/Documents/Data_Files/gdata_scalar_{0}_{1}_{2}.dat'.format( str(m),str(S),str(plot_incr) ), 'w+' )
        fdata = open( '/home/jesse/Documents/Data_Files/fdata_scalar_{0}_{1}_{2}.dat'.format( str(m),str(S),str(plot_incr) ), 'w+' )
        
        for i in arange(0,len(x)):
            gdata.write( '{0:<9}{1:^1}\n'.format( str(x[i]),str(gnorm_wavefunc[i]) ) )
            fdata.write( '{0:<9}{1:^1}\n'.format( str(x[i]),str(fnorm_wavefunc[i]) ) )
        
        plot_incr += 1
        if plot_incr > 3: break # only want first 3 states
        
show()
