import numpy
import sys
from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse.ext.bridge import bridge
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
from matplotlib import pyplot
from amuse.ic.kingmodel import new_king_model
from amuse.support.literature import LiteratureReferencesMixIn
from datetime import datetime
import pickle
import matplotlib as mpl

from galpy.potential import *
from galpy.util import bovy_conversion
import galpy as gp
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from galpy.orbit import Orbit
from operator import itemgetter
from amuse.lab import *

conv = 28.12458453897696

args = sys.argv

cluster = args[1]
N = int(args[2])
dwarves = int(args[3])

print("Started")
print("Including dwarves: ", dwarves==True)

ts = np.linspace(0., -12*conv, 1201)
dt = (ts[1]-ts[0])/250

# SAGITTARIUS
sgr_M = 1.4E10*u.Msun
sgr_a = 7*u.kpc # 6.5kpc from https://arxiv.org/pdf/1608.04743.pdf
sgr_pot = HernquistPotential(amp=sgr_M, a=sgr_a, ro=8., vo=220.)
sgr_orbit = Orbit(vxvv=[283.8313, -30.5453, 26.67, -2.692, -1.359, 140], radec=True, ro=8.,vo=220.)
sgr_orbit.integrate(ts, MWPotential2014, dt=dt)
sgr_potential = MovingObjectPotential(sgr_orbit, pot=sgr_pot)

# LMC
lmc_M = 1E11*u.Msun # https://arxiv.org/pdf/1608.04743.pdf
lmc_a = 10.2*u.kpc  # https://iopscience.iop.org/article/10.1088/2041-8205/721/2/L97/pdf
lmc_pot = HernquistPotential(amp=lmc_M, a=lmc_a, ro=8., vo=220.)
lmc_orbit = Orbit(vxvv=[78.77, -69.01, 50.1, 1.850, 0.234, 262.2], radec=True, ro=8.,vo=220.)
lmc_orbit.integrate(ts, MWPotential2014, dt=dt)
lmc_potential = MovingObjectPotential(lmc_orbit, pot=lmc_pot)

# SMC
smc_M = 2.6E10*u.Msun # https://iopscience.iop.org/article/10.1088/2041-8205/721/2/L97/pdf
smc_a = 3.6*u.kpc     # https://iopscience.iop.org/article/10.1088/2041-8205/721/2/L97/pdf
smc_pot = HernquistPotential(amp=smc_M, a=smc_a, ro=8., vo=220.)
smc_orbit = Orbit(vxvv=[16.26, -72.42, 62.8, 0.797, -1.220, 145.6], radec=True, ro=8.,vo=220.)
smc_orbit.integrate(ts, MWPotential2014, dt=dt)
smc_potential = MovingObjectPotential(smc_orbit, pot=smc_pot)

# FORNAX
fnx_M = 2E9*u.Msun # https://academic.oup.com/mnras/article/368/3/1073/1022509
fnx_a = 3.4*u.kpc  # https://academic.oup.com/mnras/article-pdf/368/3/1073/18665310/mnras0368-1073.pdf
fnx_pot = HernquistPotential(amp=fnx_M, a=fnx_a)
fnx_orbit = Orbit(vxvv=[39.962, -34.511, 139.3, 0.374, -0.401, 55.3], radec=True, ro=8.,vo=220.)
fnx_orbit.integrate(ts, MWPotential2014, dt=dt)
fnx_potential = MovingObjectPotential(fnx_orbit, pot=fnx_pot)

# DRACO
dra_M = 6.3E9*u.Msun # https://iopscience.iop.org/article/10.1086/521543/pdf
dra_a = 7*u.kpc # https://iopscience.iop.org/article/10.1086/521543/pdf
dra_pot = HernquistPotential(amp=dra_M, a=dra_a)
dra_orbit = Orbit(vxvv=[260.06, 57.965, 79.07, -0.012, -0.158, -291], radec=True, ro=8.,vo=220.)
dra_orbit.integrate(ts, MWPotential2014, dt=dt)
dra_potential = MovingObjectPotential(dra_orbit, pot=dra_pot)

# URSA MINOR
umi_M = 2.5E9*u.Msun # https://iopscience.iop.org/article/10.1086/521543/pdf
umi_a = 5.4*u.kpc    # https://iopscience.iop.org/article/10.1086/521543/pdf
umi_pot = HernquistPotential(amp=umi_M, a=umi_a)
umi_orbit = Orbit(vxvv=[227.242, 67.222, 75.86, -0.184, 0.082, -246.9], radec=True, ro=8.,vo=220.)
umi_orbit.integrate(ts, MWPotential2014, dt=dt)
umi_potential = MovingObjectPotential(umi_orbit, pot=umi_pot)

satellites_pot = [sgr_potential, lmc_potential, smc_potential, fnx_potential, dra_potential, umi_potential]
total_pot = [sgr_potential, lmc_potential, smc_potential, fnx_potential, dra_potential, umi_potential, MWPotential2014]

print("Satellites done")

class ClusterArray:
    def __init__(self, time=-12., steps=1000, include_Sgr=False, method='symplec4_c'):
        self.time = time
        self.steps = steps
        self.include_Sgr = include_Sgr
        self.method = method
        self.clusters = []
        self.vxvvs = []
        self.names = []
        if include_Sgr:
            self.potential = total_pot
        else:
            self.potential = MWPotential2014

    def load_from_names(self, names):
        self.names = names
        self.orbits = Orbit.from_name(names, solarmotion=[-11.1,24.,7.25])
        self.vxvvs = self.orbits.vxvv
        self.ctrl_orbits = Orbit.from_name(names, solarmotion=[-11.1,24.,7.25])
        self.integrate_orbits()

    def integrate_orbits(self):
        # Integrate orbit and save values
        ts=np.linspace(0, self.time/bovy_conversion.time_in_Gyr(ro=8.,vo=220.), self.steps)
        self.orbits.integrate(ts, self.potential, method=self.method, dt=dt)
        self.ctrl_orbits.integrate(ts, MWPotential2014, method=self.method, dt=dt)

    def get_rs(self, times=np.array([0])):
        o = self.orbits.getOrbit()

        index = abs(times/self.time*self.steps)
        index = index.astype(int)
        pars = o[:, index]
        new_orbit = Orbit(vxvv=pars, ro=8.,vo=220., solarmotion=[-11.1,24.,7.25])
        ts=np.linspace(0,-12.,self.steps)*u.Gyr # Gyr
        new_orbit.integrate(ts, self.potential, dt=dt)
        bovy_t = times/bovy_conversion.time_in_Gyr(ro=8, vo=220)
        times = np.linspace(0, -3./bovy_conversion.time_in_Gyr(ro=8.,vo=220.), self.steps)
        return np.array([new_orbit.R(t) for t in times])

    def get_orbits(self, times=np.array([0]), ctrl=False):
        o = self.ctrl_orbits.getOrbit() if ctrl else self.orbits.getOrbit()

        index = abs(times/self.time*self.steps)
        index = index.astype(int)
        pars = o[:, index]
        new_orbit = Orbit(vxvv=pars, ro=8.,vo=220., solarmotion=[-11.1,24.,7.25])
        ts=np.linspace(0,-6.,self.steps)*u.Gyr # Gyr
        new_orbit.integrate(ts, self.potential, dt=dt)
        return new_orbit

class galpy_profile(LiteratureReferencesMixIn):
    """
    User-defined potential from GALPY

    .. [#] Bovy, J; ApJSS, Volume 216, Issue 2, article id. 29, 27 pp. (2015)

    """
    def __init__(self,pot, t = 0. | units.Gyr, tgalpy = 0. | units.Gyr, ro=8, vo=220.):
        LiteratureReferencesMixIn.__init__(self)
        self.ro=ro
        self.vo=vo
        self.fconv = bovy_conversion.force_in_kmsMyr(ro=self.ro,vo=self.vo)
        self.pot = pot
        self.model_time=t
        self.tstart = tgalpy
        self.tgalpy = tgalpy.value_in(units.Gyr)*conv
        #self.tgalpy=tgalpy.value_in(units.Gyr)*conv

    def evolve_model(self,t_end):
        self.model_time=t_end
        self.tgalpy=(self.tstart+t_end).value_in(units.Gyr)*conv

    def get_potential_at_point(self,eps,x,y,z):

        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))

        pot=gp.potential.evaluatePotentials(self.pot,R/self.ro,zed/self.ro,phi=phi,t=self.tgalpy,ro=self.ro,vo=self.vo) | units.km**2*units.s**-2

        return pot

    def get_gravity_at_point(self,eps,x,y,z, t=0|units.Gyr):

        R=np.sqrt(x.value_in(units.kpc)**2.+y.value_in(units.kpc)**2.)
        zed=z.value_in(units.kpc)
        phi=np.arctan2(y.value_in(units.kpc),x.value_in(units.kpc))
        Rnorm = R/self.ro
        znorm = zed/self.ro

        Rforce=gp.potential.evaluateRforces(self.pot,Rnorm,znorm,phi=phi,t=self.tgalpy)
        phiforce=gp.potential.evaluatephiforces(self.pot,Rnorm,znorm,phi=phi,t=self.tgalpy)/(Rnorm)
        zforce=gp.potential.evaluatezforces(self.pot,Rnorm,znorm,phi=phi,t=self.tgalpy)

        cosphi = np.cos(phi)
        sinphi = np.sin(phi)

        ax=(Rforce*cosphi-phiforce*sinphi)*self.fconv | units.kms * units.myr**-1
        ay=(Rforce*sinphi+phiforce*cosphi)*self.fconv | units.kms * units.myr**-1
        az=zforce*self.fconv | units.kms * units.myr**-1

        return ax,ay,az

    def circular_velocity(self,r):
        vcirc=gp.potential.vcirc(self.pot,r.value_in(units.kpc)/self.ro,phi=0,ro=self.ro,vo=self.vo) | units.kms
        return vcirc

    def stop(self):
        pass

def setup_cluster(nbodycode,N,Mcluster,Rcluster,Rinit,Vinit, parameters = []):

    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    stars=new_plummer_sphere(N,converter)

    stars.x += Rinit[0]
    stars.y += Rinit[1]
    stars.z += Rinit[2]
    stars.vx += Vinit[0]
    stars.vy += Vinit[1]
    stars.vz += Vinit[2]

    code=nbodycode(converter, number_of_workers=1)
    for name,value in parameters:
        setattr(code.parameters, name, value)
    code.particles.add_particles(stars)
    return code

def evolve_cluster_in_galaxy(N, pot, Rinit, Vinit, tstart, tend, timestep, M, R):
    galaxy_code = galpy_profile(pot, tgalpy=tstart)

    cluster_code = setup_cluster(BHTree,
                                 N,
                                 Mcluster,
                                 Rcluster,
                                 Rinit,
                                 Vinit,
                                 parameters=[("epsilon_squared",(0.1 | units.parsec)**2),
                                             ("timestep", 0.1 | units.Myr),
                                             ("opening_angle", 0.6)])


    stars = cluster_code.particles.copy()
    particles = []

    channel_from_stars_to_cluster_code=stars.new_channel_to(cluster_code.particles, attributes=["x", "y", "z", "vx", "vy", "vz"])
    channel_from_cluster_code_to_stars=cluster_code.particles.new_channel_to(stars, attributes=["x", "y", "z", "vx", "vy", "vz"])

    system = bridge(verbose=False)
    system.add_system(cluster_code, (galaxy_code,))
    system.add_system(galaxy_code)

    times = quantities.arange(0|units.Myr, tend-tstart, timestep)
    for i,t in enumerate(times):
        if (t.value_in(units.Myr)%10 < timestep.value_in(units.Myr)):
            print(t)
            particles.append([cluster_code.particles.position, cluster_code.particles.velocity])
        system.evolve_model(t, timestep=timestep)

    particles.append([cluster_code.particles.position, cluster_code.particles.velocity])
    cluster_code.stop()

    return particles

t_max = -12.
times = np.array([0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -11.999])
n_steps = 1201

names = [cluster]
ca = ClusterArray(include_Sgr=True, time=t_max, steps=n_steps)
ca_ctrl = ClusterArray(include_Sgr=False, time=t_max, steps=n_steps)
ca.load_from_names(names)
ca_ctrl.load_from_names(names)

os = ca.orbits
os_ctrl = ca.ctrl_orbits

print("Cluster done")

o = os.getOrbit() if dwarves else os_ctrl.getOrbit()
pars = o[0, -1]
new_orbit = Orbit(vxvv=pars, ro=8.,vo=220.)
ts_plus=np.linspace(0,12,1201)*u.Gyr # Gyr
new_orbit.integrate(ts_plus, MWPotential2014, dt=abs(dt))

a = [new_orbit.r(t) for t in ts_plus]
ind = min(enumerate(a), key=itemgetter(1))[0]
loc = new_orbit.getOrbit()[ind]*[8, 220, 220, 8, 220, 1]
rt_peri = gp.potential.rtide(MWPotential2014, R=loc[0], z=loc[3], phi=loc[5], M=10**5.5*u.Msun)

o = os
o_ctrl = os_ctrl

t = -12.*u.Gyr

orb = o if dwarves else o_ctrl

Rinit=[orb.x(t), orb.y(t),orb.z(t)] | units.kpc
Vinit=[orb.vx(t),orb.vy(t),orb.vz(t)] | units.km/units.s

print('Tidal radius at pericentre:', rt_peri)

r_h = 0.145*rt_peri
r_p = r_h/1.3

timestep= 0.5 | units.Myr
starttime = -12000 | units.Myr
endtime = -1. | units.Myr
Mcluster = 2*10**5 | units.MSun
Rcluster = r_p | units.kpc
print("Integrating", names[0])

poten = total_pot if dwarves else MWPotential2014
#Set Galactic Potential - note that initial galpy time can be set to a different value than model_time
# galaxy_code = galpy_profile(poten, tgalpy = tstart)
start_time = datetime.now()
particles = evolve_cluster_in_galaxy(N, poten, Rinit, Vinit, starttime, endtime, timestep,
                                    Mcluster, Rcluster)
print("Done.")
print("Took:", datetime.now()-start_time)

name = "sim_"+cluster+"_"+str(dwarves)+".pickle"
with open(name, 'wb') as handle:
    pickle.dump(particles, handle, protocol=pickle.HIGHEST_PROTOCOL)

print("Pickling done, exiting")
