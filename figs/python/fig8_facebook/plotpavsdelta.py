#!/usr/bin/env python
# -*- coding: utf-8 -*-


from pylab import *
import matplotlib.font_manager as fm
prop=fm.FontProperties(size=12)

fig = figure()

ax = fig.add_subplot(1,1,1) 
name='pavsdelta'
data = loadtxt(name+'.dat')
pa = data[:,0]
deltaBA = data[:,1]
deltaFB = data[:,2]
ax.plot(pa,deltaBA,'-o', label='BA', linewidth=1.5, color='red')
ax.plot(pa,deltaFB,'--^', label='FB', linewidth=2, color='green', markeredgewidth=0)
ax.set_xlim(0,7)
ax.set_ylabel(r'$\delta$', size = '20')
ax.set_xlabel(r'p.a.', size = '20')

leg=ax.legend(loc='upper left', prop=prop)
leg.draw_frame(0)


savefig('pavsdelta.eps')






