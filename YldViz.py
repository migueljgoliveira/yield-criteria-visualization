"""
	YldViz - Yield Criteria Visualization

	Created by: M.G. Oliveira

	Last Updated: 02/2021
"""

# IMPORT PACKAGES
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from math import sin, cos, pi, degrees
from matplotlib import cm, gridspec

# IMPORT F2PY MODULES
from ummdp.ummdp_yfunc_mises import jancae_mises
from ummdp.ummdp_yfunc_hill48 import jancae_hill_1948
from ummdp.ummdp_yfunc_yld2000_2d import jancae_yld2000_2d
from ummdp.ummdp_yfunc_yld2004_18p import jancae_yld2004_18p
from ummdp.ummdp_yfunc_cpb2006 import jancae_cpb2006
from ummdp.ummdp_yfunc_bbc2005 import jancae_bbc2005

# ----------------------------------------------------------------- YIELD PLOT
def yld_plot(shear,locus0,locus1,locus2,aniso):

	plt.figure(figsize=[12+12*0.15, 8])
	gs = gridspec.GridSpec(2,6,wspace=1.5,hspace=0.5)

	ax0 = plt.subplot(gs[0,0:2])
	ax0.set_title('yield locus\nxx-yy plane')
	ax0.plot([-1.5,1.5],[0,0],'-k',lw=1.0)
	ax0.plot([0,0],[-1.5,1.5],'-k',lw=1.0)
	for k in range(len(shear)):
		ax0.plot(locus0[:,0,k],locus0[:,1,k],'b-',lw=1.0)
	ax0.set_xlim(-1.5,1.5)
	ax0.set_ylim(-1.5,1.5)
	ax0.set_xlabel('longitudinal stress, $\sigma_{xx}$')
	ax0.set_ylabel('transversal stress, $\sigma_{yy}$')
	ax0.grid()
	
	ax1 = plt.subplot(gs[0,2:4])
	ax1.set_title('yield locus\nxx-xy plane')
	ax1.plot([-1.5,1.5],[0,0],'-k',lw=1.0)
	ax1.plot([0,0],[-1.5,1.5],'-k',lw=1.0)
	ax1.plot(locus1[:,0],locus1[:,2],'b-',lw=1.0)
	ax1.set_xlim(-1.5,1.5)
	ax1.set_ylim(-1.5,1.5)
	ax1.set_xlabel('longitudinal stress, $\sigma_{xx}$')
	ax1.set_ylabel('shear stress, $\sigma_{xy}$')
	ax1.grid()

	ax2 = plt.subplot(gs[0,4:6])
	ax2.set_title('yield locus\nyy-xy plane')
	ax2.plot([-1.5,1.5],[0,0],'-k',lw=1.0)
	ax2.plot([0,0],[-1.5,1.5],'-k',lw=1.0)
	ax2.plot(locus2[:,1],locus2[:,2],'b-',lw=1.0)
	ax2.set_xlim(-1.5,1.5)
	ax2.set_ylim(-1.5,1.5)
	ax2.set_xlabel('transversal stress, $\sigma_{yy}$')
	ax2.set_ylabel('shear stress, $\sigma_{xy}$')
	ax2.grid()

	ang = np.linspace(0,90,len(aniso[:,0]))

	ax3 = plt.subplot(gs[1,1:3])
	ax3.set_title('yield stress anisotropy')
	ax3.plot([0,90],[1.0,1.0],'-k',lw=1.0)
	ax3.plot(ang,aniso[:,0],'b-')
	ax3.set_ylim(0,3)
	ax3.set_xlim(0,90)
	ax3.set_xticks(np.linspace(0,90,7))
	ax3.set_xlabel('tensile orientation, $\\theta$')
	ax3.set_ylabel('yield stress, $\sigma_{y}$')
	ax3.grid()

	ax4 = plt.subplot(gs[1,3:5])
	ax4.set_title('Lankford coefficient anisotropy')
	ax4.plot([0,90],[1.0,1.0],'-k',lw=1.0)
	ax4.plot(ang,aniso[:,1],'b-')
	ax4.set_ylim(0,3)
	ax4.set_xlim(0,90)
	ax4.set_xticks(np.linspace(0,90,7))
	ax4.set_xlabel('tensile orientation, $\\theta$')
	ax4.set_ylabel('Lankford coefficient, $r$')
	ax4.grid()

	plt.show()

# -------------------------------------------------------------- YIELD WARNING
def yld_warning(yld,np,mp):
	
	text = 'Error: {:s} takes {:d} parameter, {:d} given.'.format(yld,np,mp)
	print(text)
	exit()

# ------------------------------------------------------------------ YIELD AUX
def yld_aux(x,s,ndyld,pryld,yld,nreq,flag):

	# STRESS TENSOR COMPONENTS
	stress = np.zeros(6)
	if flag in [0,-1]:
		stress[0] = s[0]*x
		stress[1] = s[1]*x
		stress[3] = s[3]
	elif flag == 1:
		stress[0] = s[0]*x
		stress[1] = s[1]
		stress[3] = s[3]*x
	elif flag == 2:
		stress[0] = s[0]
		stress[1] = s[1]*x
		stress[3] = s[3]*x

	# CALL YIELD FUNCTION FORTRAN SUBROUTINE
	se,dseds = yld(stress,ndyld,pryld,nreq)

	# IF YIELD LOCUS
	if flag >= 0:
		return se - 1.0

	# IF ANISOTROPY
	elif flag == -1:
		return se,dseds
		
# ------------------------------------------------------------- STRESS TENSOR 1
def stress_tensor(ang,s12,flag):

	# YIELD LOCUS XX-YY
	if flag == 0:
		xx = cos(ang)
		yy = sin(ang)
		xy = s12

	# YIELD LOCUS XX-XY
	elif flag == 1:
		if 3*pi/2 >= ang > pi/2:
			xx = -cos(ang)**2
		else:
			xx = cos(ang)**2
		yy = 0.0
		xy = sin(ang)*cos(ang)
		if ang == pi/2:
			xy = 1.0
		elif ang == 3*pi/2:
			xy = -1.0

	# YIELD LOCUS YY-XY
	elif flag == 2:
		xx = 0.0
		if 3*pi/2 >= ang > pi/2:
			yy = -sin(ang)**2
		else:
			yy = sin(ang)**2
		xy = sin(ang)*cos(ang)

	# ANISOTROPY
	elif flag == -1:
		xx = cos(ang)**2
		yy = sin(ang)**2
		xy = sin(ang)*cos(ang)

	s = np.array([xx, yy, 0.0, xy, 0.0, 0.0])

	return s

# ---------------------------------------------------------------- YIELD LOCUS
def yld_locus(shear,ndyld,pryld,yld):
	nreq = 0
	npts = 10000

	# YIELD LOCUS XX-YY
	flag = 0
	ang = np.linspace(0.0,2.0*pi,npts)
	locus0 = np.zeros((npts,3,len(shear)))

	for k in range(len(shear)):
		for i in range(npts):
			s = stress_tensor(ang[i],shear[k],flag)

			alpha = fsolve(yld_aux,1.0,args=(s,ndyld,pryld,yld,nreq,flag))

			locus0[i,0,k] = s[0]*alpha
			locus0[i,1,k] = s[1]*alpha
			locus0[i,2,k] = s[3]

	# YIELD LOCUS XX-XY
	flag = 1

	ang1 = np.linspace(0,pi/2,int(round(npts/4)))
	ang2 = np.flip(np.linspace(pi,3*pi/2,int(round(npts/4)),endpoint=False))
	ang3 = np.linspace(pi,pi/2,int(round(npts/4)),endpoint=False)
	ang4 = np.linspace(3*pi/2,2*pi,int(round(npts/4)))
	ang = np.concatenate((ang1,ang2,ang3,ang4))

	maxidp = np.where(ang == pi/2)
	maxidn = np.where(ang == 3*pi/2)

	locus1 = np.zeros((len(ang),3))

	for i in range(len(ang)):
		s = stress_tensor(ang[i],0,flag)

		alpha = fsolve(yld_aux,1.0,args=(s,ndyld,pryld,yld,nreq,flag))
		locus1[i,0] = s[0]*alpha
		locus1[i,1] = s[1]
		locus1[i,2] = s[3]*alpha

	# YIELD LOCUS YY-XY
	flag = 2

	ang1 = np.flip(np.linspace(0,pi/2,int(round(npts/4))))
	ang2 = np.flip(np.linspace(3*pi/2,pi,int(round(npts/4)),endpoint=False))
	ang3 = np.flip(np.linspace(pi,pi/2,int(round(npts/4)),endpoint=False))
	ang4 = np.flip(np.linspace(3*pi/2,2*pi,int(round(npts/4)),endpoint=False))
	ang = np.concatenate((ang1,ang2,ang3,ang4))

	locus2 = np.zeros((len(ang),3))

	for i in range(len(ang)):
		s = stress_tensor(ang[i],0,flag)

		alpha = fsolve(yld_aux,1.0,args=(s,ndyld,pryld,yld,nreq,flag))

		locus2[i,0] = s[0]
		locus2[i,1] = s[1]*alpha
		locus2[i,2] = s[3]*alpha

	maxid1 = int(round(npts/4))-1
	maxid2 = int(round(npts/4))*3-1

	locus2[maxid1,1] = 0.0
	locus2[maxid1,2] = locus1[maxidp,2]
	locus2[maxid2,1] = 0.0
	locus2[maxid2,2] = locus1[maxidn,2]
	locus2[-1,1] = locus2[0,1]
	locus2[-1,2] = locus2[0,2]

	return locus0,locus1,locus2

# ------------------------------------------------------------ YIELD ANISOTROPY
def yld_aniso(ndyld,pryld,yld):
	npts = 91
	nreq = 1
	flag = -1

	ang = np.linspace(0.0,pi/2,npts)
	aniso = np.zeros((npts,2))

	for i in range(npts):
		s = stress_tensor(ang[i],0,flag)

		se,dseds = yld_aux(1.0,s,ndyld,pryld,yld,nreq,flag)

		sy = 1.0/se
		aniso[i,0] = round(sy,5)

		j = int(len(dseds)/3)+1
		num1 = dseds[0]*sin(ang[i])**2
		num2 = dseds[1]*cos(ang[i])**2
		num3 = dseds[j]*sin(ang[i])*cos(ang[i])
		num = num3 - num1 - num2

		den = dseds[0] + dseds[1]

		#r = num/den
		r = (se/den) - 1.0
		aniso[i,1] = round(r,5)

	return aniso

# -------------------------------------------------------------- YIELD WRAPPER
def yld_wrapper(yldid,ndyld):

	# von Mises
	if yldid == 0:
		yld = jancae_mises
		n = 1
		if ndyld != n:
			yld_warning('von Mises',n,ndyld)

	# Hill 1948
	elif yldid == 1:
		yld = jancae_hill_1948
		n = 6
		if ndyld != n:
			yld_warning('Hill 1948',n,ndyld)

	# Yld2004-18p
	elif yldid == 2:
		yld = jancae_yld2004_18p
		n = 19
		if ndyld != n:
			yld_warning('Yld2004-18p',n,ndyld)
	
	# CPB 2006
	elif yldid == 3:
		yld = jancae_cpb2006
		n = 14
		if ndyld != n:
			yld_warning('CPB 2006',n,ndyld)

	# Yld2000-2d
	elif yldid == -2:
		yld = jancae_yld2000_2d
		n = 9
		if ndyld != n:
			yld_warning('Yld2000-2d',n,ndyld)

	# BBC 2005
	elif yldid == -4:
		yld = jancae_bbc2005
		n = 9
		if ndyld != n:
			yld_warning('BBC 2005',n,ndyld)

	return yld

# ----------------------------------------------------------------------- MAIN
def main():

	# IMPORT INPUT - SHEAR LEVELS AND PARAMETERS
	with open('YldPlot.inp','r') as f:
		shear = f.readline().strip().split(',')
		params = f.readline().strip().split(',')

	shear = [float(s) for s in shear]	
	params = [float(p) for p in params]

	yldid = params[0]
	pryld = params[1:]
	ndyld = len(pryld)

	# YIELD WRAPPER
	yld = yld_wrapper(yldid,ndyld)
	ndyld = len(pryld)

	# YIELD LOCUS
	locus0,locus1,locus2 = yld_locus(shear,ndyld,pryld,yld)

	# ANISOTROPY
	if yldid != 0:
		aniso = yld_aniso(ndyld,pryld,yld)
	else:
		aniso = np.ones((2,2))

	yld_plot(shear,locus0,locus1,locus2,aniso)

if __name__ == "__main__":
	main()