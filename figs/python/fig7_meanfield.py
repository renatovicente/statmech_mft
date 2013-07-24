import pylab
from numpy import arange
from math import *
from scipy import integrate
from scipy import optimize
import matplotlib.pyplot as plt

#***********************************************************************
def B(theta,kalpha,delta,m,r):
	a = (1. + delta)/2.
	b = (1. - delta)/2.
	x = a*m*cos(theta)-b*r*abs(cos(theta))
	return(exp(kalpha*x))

#***********************************************************************
def Z(kalpha,delta,m,r):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*sin(x)**3,0,pi)	
	return(res[0])
	
#***********************************************************************	
def m_int(m,r,kalpha,delta,Z):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*cos(x)*sin(x)**3,0,pi)	
	return(res[0]/Z)
	
#***********************************************************************	
def m2_int(m,r,kalpha,delta,Z):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*(cos(x)**2)*sin(x)**3,0,pi)	
	return(res[0]/Z)
	
#***********************************************************************	
def r_int(m,r,kalpha,delta,Z):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*abs(cos(x))*sin(x)**3,0,pi)	
	return(res[0]/Z)	
	
#***********************************************************************
def mf_eqs(p,kalpha,delta):
	m,r=p
	cte=Z(kalpha,delta,m,r)
	return(m_int(m,r,kalpha,delta,cte),r_int(m,r,kalpha,delta,cte))

#***********************************************************************		
def mr(kalpha,delta):
	p=(0.5,0.5)
	eps=1.
	n=0
	while (eps>1e-6 and n<1000):
		pnew=mf_eqs(p,kalpha,delta)
		eps=max(fabs(pnew[0]-p[0]),fabs(pnew[1]-p[1]))
		p=pnew
		n=n+1
	return(p,eps,n)

#***********************************************************************
def Cph(h,kalpha,delta,m,r):
	a = (1.+delta)/2.
	b = (1.-delta)/2.
	return((1.-h**2)*exp(kalpha*(a*m*h-b*r*abs(h))))

def C(kalpha,delta,m,r):
	res=integrate.quad(lambda x:Cph(x,kalpha,delta,m,r),-1.,1.)	
	return(res[0])

#***********************************************************************
def ph(h,kalpha,delta,m,r,C):
	return(Cph(h,kalpha,delta,m,r)/C)	
	
def ph_aprox(h,kalpha,delta,m):
	dkam=delta*kalpha*m	
	return(0.5*(dkam**2)*(1.-h**2)*exp(-dkam*(1.-h)))	
	
def avg_h_aprox(kalpha,delta):
	res=mr(kalpha,delta)
	m=res[0][0]
	dkam=delta*kalpha*m
	h=(exp(-dkam)/(dkam)**2)*(-6.+(dkam)**2 + 2*exp(dkam)*(3+dkam(dkam-3.)))
	return(h)

def vm(m,r,kalpha,delta):
    nZ=Z(kalpha,delta,m,r)
    m2=m2_int(m,r,kalpha,delta,nZ)
    return(m2-m**2)
    
def E(m,r,kalpha,delta,h):
    return(-0.5*(1+delta)*h*m+0.5*(1-delta)*abs(h)*r)
    
def spec_heat(m,r,kalpha,delta):
    nC=C(kalpha,delta,m,r)
    avgE=integrate.quad(lambda x:E(m,r,kalpha,delta,x)*ph(x,kalpha,delta,m,r,nC),-1.,1.)
    avgE2=integrate.quad(lambda x:(E(m,r,kalpha,delta,x)**2)*ph(x,kalpha,delta,m,r,nC),-1.,1.)
    return((kalpha*kalpha)*(avgE2[0]-avgE[0]*avgE[0]))
    
#***********************************************************************           
params = {'backend': 'ps',
          'axes.labelsize': 24,
          'text.fontsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 22,
          'ytick.labelsize': 22,
          'text.usetex': True}
pylab.rcParams.update(params)

fig = plt.figure(figsize=(9,6))
fig.subplots_adjust(hspace=0.32,wspace=0.34,top=0.93, right=0.94, bottom=0.10)

#***********************************************************************
fa = fig.add_subplot(2,2,1)
fb = fig.add_subplot(2,2,2)

delta=arange(0.,1.2,0.2)
kalpha=arange(0.,30.,0.1)
for d in delta:
    m=[]
    varm=[]
    for k in kalpha:
        res=mr(k,d)
        m.append(res[0][0])
        varm.append(vm(res[0][0],res[0][1],k,d))
		
    fa.plot(kalpha,m)
    fa.hold(True)
    fb.plot(kalpha,varm)
    fb.hold(True)
    
fa.set_xlabel(r'$k\alpha$',fontsize=20)
fa.set_ylabel(r'$m$',fontsize=20)
fa.text(2.5,0.1,'a',fontsize=24,backgroundcolor='w')

fb.set_xlabel(r'$k\alpha$',fontsize=20)
fb.set_ylabel(r'$v_m$',fontsize=20)
fb.text(2.5,0.025,'b',fontsize=24,backgroundcolor='w')

#***********************************************************************
fc = fig.add_subplot(2,2,3)
delta=arange(0.,1.01,0.01)
kalpha=[]
k0=20.
for d in delta:
       kn=optimize.bisect((lambda x:mr(x,d)[0][0]-0.0001),0.,k0)
       kalpha.append(kn)
       k0=kn 
       #print d,kn 
fc.plot(delta,kalpha)
fc.set_ylabel(r'$k\alpha$',fontsize=20)
fc.set_xlabel(r'$\delta$',fontsize=20)
fc.set_xlim((0.,1.))
fc.set_ylim((0.,20.))
fc.grid(True)
fc.text(.1,2.,'c',fontsize=24,backgroundcolor='w')
fc.text(.3,2.5,r'$m=0$',fontsize=24,backgroundcolor='w')
fc.text(.5,12,r'$m>0$',fontsize=24,backgroundcolor='w')

#***********************************************************************
fd = fig.add_subplot(2,2,4)
delta=arange(0.,1.01,0.01)
kalpha=[15.,20.,25.,30.,35.,40.,45.,50.]
for k in kalpha:
    varm=[]
    for d in delta:
        res=mr(k,d)
        varm.append(vm(res[0][0],res[0][1],k,d))
    fd.plot(delta,varm)
    fd.hold(True)
fd.set_xlabel(r'$\delta$',fontsize=20)
fd.set_ylabel(r'$v_m$',fontsize=20)
fd.text(.1,0.01,'d',fontsize=24,backgroundcolor='w')

#***********************************************************************
plt.savefig('fig7_meanfield.png')
plt.show()
