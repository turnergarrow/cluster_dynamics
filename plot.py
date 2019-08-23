import pickle
import numpy as np
from matplotlib import pyplot as plt
from astropy import units as u
from galpy.util import bovy_conversion
import matplotlib as mpl
from amuse.lab import *
from scipy import stats
from matplotlib.animation import FuncAnimation
import galpy as gp
from galpy.orbit import Orbit
from galpy.potential import *
from amuse.lab import *
from amuse.couple import bridge
from amuse.support.literature import LiteratureReferencesMixIn
from amuse.datamodel import Particles
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class galpy_profile(LiteratureReferencesMixIn):
    """
    User-defined potential from GALPY
    
    .. [#] Bovy, J; ApJSS, Volume 216, Issue 2, article id. 29, 27 pp. (2015)
    
    """
    def __init__(self,pot, t = 0. | units.Gyr, tgalpy = 0. | units.Gyr, ro=8, vo=220.):
        LiteratureReferencesMixIn.__init__(self)

        #Galpy scaling parameters
        self.ro=ro
        self.vo=vo

        #Galpy potential
        self.pot = pot

        #Initialize model time
        self.model_time=t

        #Initialize galpy time
        #Depending on how the galpy potential is setup, one may want to define a different initial galpy time
        #from which to evolve the potential (e.g. negative times) 
        self.tgalpy=tgalpy.value_in(units.Gyr)/bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)

    def evolve_model(self,t_end):
        dt=t_end-self.model_time
        self .model_time=t_end
        self.tgalpy+=dt.value_in(units.Gyr)/bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)

    def get_potential_at_point(self,eps,x,y,z):

        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        pot=gp.potential.evaluatePotentials(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy,ro=self.ro,vo=self.vo) | units.km**2*units.s**-2
                    
        return pot

    def get_gravity_at_point(self,eps,x,y,z):
        
        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        Rforce=gp.potential.evaluateRforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)
        phiforce=gp.potential.evaluatephiforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)/(R/self.ro)
        zforce=gp.potential.evaluatezforces(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy)

        ax=(Rforce*np.cos(phi)-phiforce*np.sin(phi))*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        ay=(Rforce*np.sin(phi)+phiforce*np.cos(phi))*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        az=zforce*bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo) | units.kms * units.myr**-1
        
        return ax,ay,az

    def mass_density(self,x,y,z):
        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        dens=gp.potential.evaluateDensities(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy,ro=self.ro,vo=self.vo) | units.MSun/(units.parsec**3.)
                        
        return dens

    def circular_velocity(self,r):
        vcirc=gp.potential.vcirc(self.pot,r.value_in(units.kpc)/self.ro,phi=0,ro=self.ro,vo=self.vo) | units.kms
        return vcirc
        
    def enclosed_mass(self,r):
        grav=4.302e-6 #kpc (km/s)^2/Msun
        vc2=gp.potential.vcirc(self.pot,r.value_in(units.kpc)/self.ro,phi=0,t=self.tgalpy,ro=self.ro,vo=self.vo)**2.
        menc= vc2*r.value_in(units.kpc)/grav | units.MSun
        return menc

    def stop(self):
        pass

with open('sim_ESO28006_0.pickle', 'rb') as handle:
    particles_MW = pickle.load(handle)
with open('sim_ESO28006_1.pickle', 'rb') as handle:
    particles_DG = pickle.load(handle)

ctrl = True
number_of_frames = 1200
parts = particles_MW if ctrl else particles_DG

def update_hist(i):
    plt.cla()
    print(i)
    
    plt.hist2d([i[0].value_in(units.kpc) for i in parts[i][0]], 
                         [i[1].value_in(units.kpc) for i in parts[i][0]],
                         bins=100,
                         range=((-20, 20),(-20, 20)))

fig = plt.figure()
hist = plt.hist2d([i[0].value_in(units.kpc) for i in parts[0][0]], 
                    [i[1].value_in(units.kpc) for i in parts[0][0]],
                    bins=100,
                    range=((-20, 20),(-20, 20)))#, norm = mpl.colors.LogNorm())

anim = animation.FuncAnimation(fig, update_hist, number_of_frames)
anim.save('control.gif', writer='imagemagick', fps=30)
