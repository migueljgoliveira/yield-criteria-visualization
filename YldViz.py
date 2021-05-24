"""
    YldViz - Yield Criteria Visualization

    Created by: M.G. Oliveira

    Last Updated: 05/2021
"""

# IMPORT PACKAGES
import warnings
import numpy as np
from math import sin, cos, pi
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import fsolve
from mpl_toolkits.mplot3d import axes3d


# PRE-LOADED COMMANDS
plt.style.use('seaborn')
warnings.filterwarnings("ignore")

# IMPORT F2PY FILES
import ummdp

def yld_3d(ndyld,pryld,yld):
    nreq = 0
    npts = 100

    # EXTRACT MAXIMUM SHEAR STRESS
    flag = 1

    ang = pi/2
    s = stress_tensor(ang,0,flag)
    alpha = fsolve(yld_aux,1.0,args=(s,ndyld,pryld,yld,nreq,flag))
    maxshear = s[3]*alpha

    nsh = 25
    shear = np.linspace(0.0,1.0,nsh)**0.85*maxshear

    # YIELD LOCUS XX-YY
    flag = 0
    ang = np.linspace(0.0,2.0*pi,npts)
    locus3d = np.zeros((1+npts*(len(shear)-1),3))

    for k in range(nsh-1):
        for i in range(npts):
            s = stress_tensor(ang[i],shear[k],flag)

            alpha = fsolve(yld_aux,1.0,args=(s,ndyld,pryld,yld,nreq,flag))

            idx = i+npts*k
            locus3d[idx,0] = s[0]*alpha
            locus3d[idx,1] = s[1]*alpha
            locus3d[idx,2] = s[3]

    locus3d[-1,2] = shear[-1]

    return locus3d

# ------------------------------------------------------------------- SAVE DATA
def yld_save(shear,locus0,locus1,locus2,aniso):

    nsh = 2*len(shear)
    ncol = nsh + 3

    if not isinstance(aniso,int):
        ncol = ncol + 3

    data = np.zeros((locus0.shape[0],ncol))

    # EXTRACT DATA FROM XX-YY LOCUS
    for i in range(len(shear)):
        data[:,i*2:i*2+2] = locus0[:,0:2,i]
    
    header = []
    for i in shear:
        tmp = [f'sigXX-{i}',f'sigYY-{i}'] 
        header = header + tmp
    
    # EXTRACT DATA FROM LOCUS-XY
    tmp = ['tau_sigXX','tau_sigYY','tau_sigXY']
    header = header + tmp
    data[:,nsh] = locus1[:,0]
    data[:,nsh+1] = locus2[:,1]
    data[:,nsh+2] = locus2[:,2]

    # EXTRACT DATA FROM ANISOTROPY
    if len(aniso) > 2:
        tmp = ['aniso_theta','aniso_sig0','aniso_r']
        header = header + tmp
        data[:aniso.shape[0],nsh+3:] = aniso
        data[aniso.shape[0]:,nsh+3:] = np.nan

    header = ';'.join(header)
    np.savetxt('YldViz-Data.csv',data,delimiter=';',fmt='%.6f',header=header,comments='')

    return

# ------------------------------------------------------------------ YIELD PLOT
def yld_plot(shear,locus0,locus1,locus2,aniso,locus3d):

    plt.figure(figsize=[12+12*0.15, 8])
    gs = gridspec.GridSpec(2,6,wspace=1.5,hspace=0.5)

    # YIELD LOCUS - XX-YY
    ax0 = plt.subplot(gs[0,0:2])
    ax0.set_title('yield locus\nxx-yy')
    ax0.axvline(0,c='k',lw=0.75)
    ax0.axhline(0,c='k',lw=0.75)
    for k in range(len(shear)):
        ax0.plot(locus0[:,0,k],locus0[:,1,k],'-',lw=1.0,color='tab:blue')
    ax0.set_xlim(-1.5,1.5)
    ax0.set_ylim(-1.5,1.5)
    ax0.set_xticks(np.arange(-1.5,2.0,0.5))
    ax0.set_yticks(np.arange(-1.5,2.0,0.5))
    ax0.set_xlabel('longitudinal stress, $\sigma_{xx}$')
    ax0.set_ylabel('transversal stress, $\sigma_{yy}$')
    ax0.set_aspect('equal')
    
    # YIELD LOCUS - XX-XY
    ax1 = plt.subplot(gs[0,2:4])
    ax1.set_title('yield locus\nxx-xy')
    ax1.axvline(0,c='k',lw=0.75)
    ax1.axhline(0,c='k',lw=0.75)
    ax1.plot(locus1[:,0],locus1[:,2],'-',lw=1.0,color='tab:blue')
    ax1.set_xlim(-1.5,1.5)
    ax1.set_ylim(-1.5,1.5)
    ax1.set_xticks(np.arange(-1.5,2.0,0.5))
    ax1.set_yticks(np.arange(-1.5,2.0,0.5))
    ax1.set_xlabel('longitudinal stress, $\sigma_{xx}$')
    ax1.set_ylabel('shear stress, $\sigma_{xy}$')
    ax1.set_aspect('equal')

    # YIELD LOCUS - YY-XY
    ax2 = plt.subplot(gs[0,4:6])
    ax2.set_title('yield locus\nyy-xy')
    ax2.axvline(0,c='k',lw=0.75)
    ax2.axhline(0,c='k',lw=0.75)
    ax2.plot(locus2[:,1],locus2[:,2],'-',lw=1.0,color='tab:blue')
    ax2.set_xlim(-1.5,1.5)
    ax2.set_ylim(-1.5,1.5)
    ax2.set_xticks(np.arange(-1.5,2.0,0.5))
    ax2.set_yticks(np.arange(-1.5,2.0,0.5))
    ax2.set_xlabel('transversal stress, $\sigma_{yy}$')
    ax2.set_ylabel('shear stress, $\sigma_{xy}$')
    ax2.set_aspect('equal')

    # ANISOTROPY - YIELD STRESS
    ax3 = plt.subplot(gs[1,2:4])
    ax3.set_title('anisotropy\nyield stress')
    ax3.axhline(1,c='k',lw=0.75)
    ax3.plot(aniso[:,0],aniso[:,1],'-',lw=1.0,color='tab:blue')
    ax3.set_ylim(0,3)
    ax3.set_xlim(0,90)
    ax3.set_xticks(np.linspace(0,90,7))
    ax3.set_xlabel('tensile orientation, $\\theta$')
    ax3.set_ylabel('yield stress, $\sigma_{y}$')

    # ANISOTROPY - LANKFORD COEFFICIENT
    ax4 = plt.subplot(gs[1,4:6])
    ax4.set_title('anisotropy\nLankford coefficient')
    ax4.axhline(1,c='k',lw=0.75)
    ax4.plot(aniso[:,0],aniso[:,2],'-',lw=1.0,color='tab:blue')
    ax4.set_ylim(0,3)
    ax4.set_xlim(0,90)
    ax4.set_xticks(np.linspace(0,90,7))
    ax4.set_xlabel('tensile orientation, $\\theta$')
    ax4.set_ylabel('Lankford coefficient, $r$')

    # YIELD LOCUS - 3-DIMENSIONAL
    ax5 = plt.subplot(gs[1,0:2],projection='3d')
    ax5.scatter(locus3d[:,0],locus3d[:,1],locus3d[:,2],s=1,color='tab:blue')
    ax5.set_xlim(-1.5,1.5)
    ax5.set_ylim(-1.5,1.5)
    ax5.set_zlim(0.0,1.0)
    ax5.set_xlabel('longitudinal stress, $\sigma_{xx}$')
    ax5.set_ylabel('transversal stress, $\sigma_{yy}$')
    ax5.set_zlabel('shear stress, $\sigma_{xy}$')
    ax5.set_facecolor('w')
    
    plt.show()

# --------------------------------------------------------------- YIELD WARNING
def yld_warning(yld,np,mp):
    
    text = 'Error: {:s} takes {:d} parameter, {:d} given.'.format(yld,np,mp)
    print(text)
    exit()

# ------------------------------------------------------------------- YIELD AUX
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

# ----------------------------------------------------------------- YIELD LOCUS
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

# ------------------------------------------------------------------ ANISOTROPY
def yld_aniso(ndyld,pryld,yld):
    npts = 91
    nreq = 1
    flag = -1

    ang = np.linspace(0.0,pi/2,npts)
    aniso = np.zeros((npts,3))

    for i in range(npts):
        s = stress_tensor(ang[i],0,flag)

        se,dseds = yld_aux(1.0,s,ndyld,pryld,yld,nreq,flag)

        sy = 1.0/se
        aniso[i,1] = round(sy,5)

        j = int(len(dseds)/3)+1
        num1 = dseds[0]*sin(ang[i])**2
        num2 = dseds[1]*cos(ang[i])**2
        num3 = dseds[j]*sin(ang[i])*cos(ang[i])
        num = num3 - num1 - num2

        den = dseds[0] + dseds[1]

        #r = num/den
        r = (se/den) - 1.0
        aniso[i,2] = round(r,5)

    aniso[:,0] = np.rad2deg(ang)

    return aniso

# --------------------------------------------------------------- YIELD WRAPPER
def yld_wrapper(yldid,ndyld):

    # von Mises
    if yldid == 0:
        yld = ummdp.jancae_mises
        n = 1
        if ndyld != n:
            yld_warning('von Mises',n,ndyld)

    # Hill 1948
    elif yldid == 1:
        yld = ummdp.jancae_hill_1948
        n = 6
        if ndyld != n:
            yld_warning('Hill 1948',n,ndyld)

    # Yld2004-18p
    elif yldid == 2:
        yld = ummdp.jancae_yld2004_18p
        n = 19
        if ndyld != n:
            yld_warning('Yld2004-18p',n,ndyld)
    
    # CPB 2006
    elif yldid == 3:
        yld = ummdp.jancae_cpb2006
        n = 14
        if ndyld != n:
            yld_warning('CPB 2006',n,ndyld)

    # Yld2000-2d
    elif yldid == -2:
        yld = ummdp.jancae_yld2000_2d
        n = 9
        if ndyld != n:
            yld_warning('Yld2000-2d',n,ndyld)

    # BBC 2005
    elif yldid == -4:
        yld = ummdp.jancae_bbc2005
        n = 9
        if ndyld != n:
            yld_warning('BBC 2005',n,ndyld)

    return yld

# ------------------------------------------------------------------------ MAIN
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

    # 3D LOCUS
    locus3d = yld_3d(ndyld,pryld,yld)

    # PLOT
    yld_plot(shear,locus0,locus1,locus2,aniso,locus3d)

    # SAVE
    yld_save(shear,locus0,locus1,locus2,aniso)

if __name__ == "__main__":
    main()