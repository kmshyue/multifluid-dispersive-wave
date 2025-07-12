from pylab import *
import scipy.io
import matplotlib
import linear_waves as LW

g = 9.81

RC = 300.0
lip = False

xlower = 0.
xupper = 5e3
h0 = 4000.

#
t0    = 0
tstop = 60 
nstop = 2
tout  = linspace(t0,tstop,nstop)
#
# Initial eta:

L     = xupper-xlower
mx_H1 = 5000
dx    = L/mx_H1
r     = linspace(dx/2,xupper-dx/2,mx_H1)
k     = linspace(1e-6,0.06,4000)

def eta0_crater(r,RC):
    """
    Initial conditions for a semi-circle crater with inner radius RC
    r is an array of radial distance values.  RC,r are all in meters.
    """

    eta = where(r<RC, -sqrt(RC**2 - r**2), 0.)
    return eta
    
print('Computing eta0 and Hankel transform eta0hat...')
eta0 = eta0_crater(r,RC)
eta0hat = LW.Htransform(r,eta0,k)

# set dispersion relation:
equation = 'Airy'

if equation == 'Airy':
    omega = lambda k: LW.omega_airy(k,h0)
elif equation == 'SGN':
    omega = lambda k: LW.omega_sgn(k,h0,alpha=1.153)
elif equation == 'MS':
    omega = lambda k: LW.omega_madsen(k,h0,B=1/15.)
else:
    raise InputError('unrecongized equation: %s' % equation)

# make plots of eta at various times:
#for j,t in enumerate([0., 10., 20., 30., 40., 50., 125., 250.]):
for j,t in enumerate(tout):
    eta,u = LW.eta_u_radial(t,r,k,eta0hat,omega,h0,direction='both')

    figure(j, figsize=(10,4))
    clf()
    plot(r/1e3, eta, 'b')
    plot(-r/1e3, eta, 'b')
    grid(True)
    xlim(-xupper/1e3, xupper/1e3)
    #ylim(-1.5*RC, 1.5*RC)
    xlabel('Radial distance (km)')
    ylabel('water surface eta (m)')
    title_str = '%s with RC=%s at t = %s sec' % (equation, int(RC), int(t))
    title(title_str)
    fname = '%s_RC%s_t%s.png' % (equation, int(RC), int(t))
    savefig(fname, bbox_inches='tight')
    print('Created ',fname)
#
    # save data             
    mdict = {'r':r,'eta':eta,'u':u}
    pname = '%s_RC%s_t%s.mat' % (equation, int(RC), int(t))
    scipy.io.savemat(pname,mdict)
