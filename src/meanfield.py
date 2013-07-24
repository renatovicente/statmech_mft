from math import *
from scipy import integrate
from scipy import optimize

def B(theta,kalpha,delta,m,r):
	a = (1. + delta)/2.
	b = (1. - delta)/2.
	x = a*m*cos(theta)-b*r*abs(cos(theta))
	return(exp(kalpha*x))
	
def Z(kalpha,delta,m,r):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*sin(x)**3,0,pi)	
	return(res[0])
	
def m_int(m,r,kalpha,delta,Z):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*cos(x)*sin(x)**3,0,pi)	
	return(res[0]/Z)

def r_int(m,r,kalpha,delta,Z):
	res=integrate.quad(lambda x:B(x,kalpha,delta,m,r)*abs(cos(x))*sin(x)**3,0,pi)	
	return(res[0]/Z)	

def mf_eqs(p,kalpha,delta):
	m,r=p
	cte=Z(kalpha,delta,m,r)
	return(m_int(m,r,kalpha,delta,cte),r_int(m,r,kalpha,delta,cte))
		
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

def Cph(h,kalpha,delta,m,r):
	a = (1.+delta)/2.
	b = (1.-delta)/2.
	return((1.-h**2)*exp(kalpha*(a*m*h-b*r*abs(h))))	

def C(kalpha,delta,m,r):
	res=integrate.quad(lambda x:Cph(x,kalpha,delta,m,r),-1.,1.)	
	return(res[0])

def ph(h,kalpha,delta,m,r,C):
	return(Cph(h,kalpha,delta,m,r)/C)	

def avg_h_aprox(kalpha,delta):
	res=mr(kalpha,delta)
	m=res[0][0]
	dkam=delta*kalpha*m
	h=(exp(-dkam)/(dkam)**2)*(-6.+(dkam)**2 + 2*exp(dkam)*(3+dkam(dkam-3.)))
	return(h)
	
def ph_aprox(h,kalpha,delta,m):
	dkam=delta*kalpha*m	
	return(0.5*(dkam**2)*(1.-h**2)*exp(-dkam*(1.-h)))

	
