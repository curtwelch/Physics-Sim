#!/usr/bin/python

"""
    physics-sim.py - Atomic particle simulation module

    2016-05-18 Started this file.  Took it from xnet.py from my AI work.
    2021-02-05 Created PyCharm Project, Put on Github, converted
                to python3.9, and started major cleanup of this old code.
    2021-02-21 Major code refactor and cleanup happening the past weeks.
                Rename from physics.py to physics_sim.py
                Turn into module.
                Move experiment run code into separate files.

"""

import time
import math
import numpy as np
from numpy import ndarray
from typing import Tuple, Dict, Callable
# import timeit
import pygame

import crt

BLACK = 0, 0, 0         # Color of an Electron
WHITE = 255, 255, 255   # Screen Background color
RED = 255, 0, 0         # Color of a Proton
BLUE = 0, 0, 255        # Color of a Neutron

ANGSTROM = 1.0e-10          # m/Å - One ANGSTROM is 1e-10 meters
CONST_C = 299792458.0       # Speed of light m/s - defined constant
# CONST_KE = 8.9875517873681764e9   # Coulomb's constant (1/4 πe)  K(sub)e
# CONST_KE = CONST_C² * 1.0e-7      # Same as above -- old now?
CONST_KE = 8.9875517923e9   # Coulomb's constant; New? 8.9875517923(14)×10^9
# CONST_G = 6.67408e-11  # Gravitational Constant, 2014 CODATA
CONST_G = 6.67430e-11  # Gravitational Constant, 2018 CODATA 6.67430(15)×10−11

CONST_P_CHARGE = +1.60218e-19
CONST_P_MASS = 1.672621898e-27

CONST_E_CHARGE = -CONST_P_CHARGE
CONST_E_MASS = 9.10938356e-31

# RLimit = 0.1 * ANGSTROM       # Radius limit hack
# RLimit = 0.0001 * ANGSTROM    # Radius limit hack
RLimit = 0.0000001 * ANGSTROM   # Radius limit hack

Energy_fix = False      # fix based on total PE+KE at start
Energy_fix2 = False     # TODO is this useful or should be be deleted?


def magnitude(vec: ndarray) -> float:
    """ Compute length of 3D vector.
        Just a wrapper for np.linalg.norm(vec)
        sqrt(sum of squares)
    """
    return np.linalg.norm(vec)


class ParticleState:
    def __init__(self, p: 'Particle'):
        """ Simulation State of a Particle. """
        self.r = np.zeros(3)  # x, y, z, Position
        self.v = np.zeros(3)  # x, y, z, Velocity
        self.f = np.zeros(3)  # x, y, z, Force
        self.p = p  # particle we are a state for

    def copy(self):
        """ Return copy of self. """
        ps = ParticleState(self.p)
        ps.r[:] = self.r
        ps.v[:] = self.v
        ps.f[:] = self.f
        return ps

    def __str__(self):
        return f"state: r:{self.r}\n       v:{self.v}\n       f:{self.f}"

    def momentum(self) -> ndarray:
        """ Momentum Vector. """
        return self.v * self.p.mass

    def f_from_v(self):
        """ Calculate frequency from v assuming it's in a circular orbit. """
        # Just assume it's in a circular orbit with a matching particle.
        # Given it's velocity, what is its orbital frequency?
        # Centripetal force is mv²/r. r is to center of mass.
        # Coulomb force is kqq / d² ... d is between particles
        # This math is for an Electron (flipped for Proton)
        # d = r + r * me / mp
        # d = r(1 + me/mp)
        # j = (1+me/mp)
        # d = rj
        # c force = kq²/r²j²
        # When these two are equal, we have a circular orbit.
        # kq²/r²j² = mv²/r  solving for r
        # kq²/j² = mv²r
        # r = kq²/mv²j²
        # The length of the orbit is c = πd, or 2πr.
        # cycles per second is distance/time so, d is c.
        # freq is laps per second.  We have v in m per second.
        # So we convert m/s to laps per second by dividing...
        # v / (distance of one lap).
        # freq = v / (2πr)
        # freq = v / (2πkq² / mv²j²)
        # freq = mv³j² / 2πkq²
        # v = magnitude(p1.v())
        if isinstance(self.p, Electron):
            j = 1.0 + CONST_E_MASS / CONST_P_MASS
        else:
            j = 1.0 + CONST_P_MASS / CONST_E_MASS

        # r = (p1.ke * p1.charge**2) / (p1.mass * v**2 * j**2)
        # d = r*j
        f = self.p.mass * magnitude(self.v) ** 3 * j ** 2 / \
            (2.0 * math.pi * CONST_KE * self.p.charge ** 2)
        return f

    def zero_force(self):
        """ Zero force vector. """
        self.f[:] = 0.0

    def total_force(self, p_state: 'ParticleState'):
        """ Add force on self created by p2_state. """

        if self.p is p_state.p:
            # Particles cause no force on self
            return np.zeros(3)

        if isinstance(self.p, Neutron):
            # these fake Neutrons only have forces with Protons
            # Not with each other, or not with Electrons.
            if not isinstance(p_state.p, Proton):
                return np.zeros(3)

        if isinstance(p_state.p, Neutron):
            # these fake Neutrons only have forces with Protons
            # Not with each other, or not with Electrons.
            if not isinstance(self.p, Proton):
                return np.zeros(3)

        # Add Electrostatic force
        # TODO might need to be R limited? was before
        f: ndarray = self.es_force(p_state)

        # Add electro-drag force
        if self.p.sim.do_v_force:
            f += self.v_force(p_state)

        return f

    # def add_force(self, p: 'Particle', p2_state: ParticleState):  # and set end force as well
    #     if p is self:
    #         return
    #
    #     # dx = (self.cur_state_r0 - p.cur_state_r0)
    #     # dy = (self.cur_state_r1 - p.cur_state_r1)
    #     # dz = (self.cur_state_r2 - p.cur_state_r2)
    #     dr = p2_state.r - p2_state.r
    #
    #     r2, l2 = self.distance2(p)
    #
    #     if r2 == 0.0:
    #         return  # Bogus but prevents DBZ
    #
    #     force = self.ke * (self.charge * p.charge) / l2
    #
    #     r = math.sqrt(r2)
    #
    #     # self.fx += force * dx / r
    #     # self.fy += force * dy / r
    #     # self.fz += force * dz / r
    #     self.cur_state.f += force * dr / r
    #
    #     if do_v_force:
    #         f = self.v_force(p.cur_state)
    #         # self.fx += f[0]
    #         # self.fy += f[1]
    #         # self.fz += f[2]
    #         self.cur_state.f += f
    #
    #     # self.end_fx = self.fx
    #     # self.end_fy = self.fy
    #     # self.end_fz = self.fz
    #     self.end_state.f = np.copy(self.cur_state.f)

    def es_force(self, p_state: 'ParticleState') -> ndarray:
        """ Electrostatic force on self by p2_state per Coulomb's law. """
        # Force on self, caused by p.
        # real force, not R limit limited force
        # Returns 0,0,0 instead of infinity for two particles located
        # on same spot.

        if self.p is p_state.p:
            # Particles cause no force on self
            return np.zeros(3)

        dr: ndarray = self.r - p_state.r

        r = np.linalg.norm(dr)

        dr_unit = dr / r

        if r == 0.0:
            return np.zeros(3)  # Bogus but prevents DBZ -- should be infinity

        return dr_unit * CONST_KE * self.p.charge * p_state.p.charge / (r * r)

    # def v_force_test(self, p2_state: 'ParticleState'):
    #     n = timeit.timeit('nv = self.v_force_new(p2_state)', '', number=100, globals=locals())
    #     o = timeit.timeit('ov = self.v_force_old(p2_state)', '', number=100, globals=locals())
    #     nv = self.v_force_new(p2_state)
    #     if n < o:
    #         # print(end='.')
    #         print(f"New Faster - New: {n:.3f} Old: {o:.3f}")
    #     else:
    #         print()
    #         print(f"Old Faster - New: {n:.3f} Old: {o:.3f}")
    #         ov = self.v_force_old(p2_state)
    #         print("  nv:", nv)
    #         print("  ov:", ov)
    #
    #     return nv

    def v_force(self, p_state: 'ParticleState'):
        """
            Faster version of v_force_old()
            Runs about 2x faster most the time.
            Slightly different results at the lowest bits.
        """
        # return self.v_force_old(p2_state)
        dv: ndarray = self.v - p_state.v
        dr: ndarray = self.r - p_state.r
        r = np.linalg.norm(dr)
        dr_hat = dr / r
        f_vec = (dr_hat * dv.dot(dr_hat) * CONST_KE * -1.0 *
                 np.abs(self.p.charge) *
                 np.abs(p_state.p.charge) / (r * r * CONST_C))
        return f_vec

    def v_force_old(self, p_state: 'ParticleState'):
        """ 2021-02-13 New idea.
            At least I hope it's new.  It was years ago I did the others.
            Use the velocity which the two particles are approaching to define
            the magnetic force.  Make the magnetic force act in the same line
            as the column force, but make it slow down velocity. So as to limit
            V to be the speed of light.  If V == c, then the magnetic force
            is just the opposite of the Coulomb force amd cancels it out.
        """
        # relative_v = subtract(self.v(), p.v())
        dv = self.v - p_state.v
        # r = (self.cur_state_r0 - p.cur_state_r0, self.cur_state_r1 - p.cur_state_r1, self.cur_state_r2 - p.cur_state_r2)
        # r_hat = product(1.0 / magnitude(r), r)
        # r points from p (logically at origin) to self.
        # r_hat is the unit vector pointing the same way.
        dr: ndarray = self.r - p_state.r
        r_hat: ndarray = dr / magnitude(dr)

        # Magnitude of v in line with r
        # vr = dot(relative_v, r_hat)
        vr = np.dot(dv, r_hat)
        # vr is the magnitude (and sign) of the relative velocity from
        # p to self.
        # es_f_mag = magnitude(self.es_force(p))
        es_f_mag = magnitude(self.es_force(p_state))
        # f_vec = product(es_f_mag * (-vr) / CONST_C, r_hat)
        f_vec: ndarray = (es_f_mag * -vr / CONST_C) * r_hat
        # First try at coding it:
        # v_mag = magnitude(relative_v)
        # print(f"{vr=:.3f} {v_mag=:.3f}")
        # F = self.product(-v_mag / p.c, es_f) # reduces es_force per abs(speed)
        return f_vec

    def distance_sq(self, p_state: 'ParticleState'):
        """ Distance**2 between self and p2_state.
            :returns: distance**2, limited_distance**2
        """

        d_sqr = np.sum((self.r - p_state.r) ** 2)

        return self.limited_distance_sqr(d_sqr)

    @staticmethod
    def limited_distance_sqr(d_sqr: float):
        """ Limit distance to RLimit to solve computational problems.
            return (real, limited) tuple
        """
        return d_sqr, max(d_sqr, RLimit ** 2)


class Particle:
    """ The root class for Protons and Electrons.
        Units are standard SI: meters, seconds, Newtons, Kg, Coulombs
    """

    def __init__(self, sim: 'Simulation',
                 r: Tuple[float, float, float],
                 v=(0.0, 0.0, 0.0),
                 ):
        """ The Base class for Electrons and Protons.

        Args:
            sim: Back reference to the Simulation()
            r: 3D position vector in meters.
            v: 3D position vector in m/s
        """
        self.sim = sim
        self.cur_state = ParticleState(self)
        self.end_state = ParticleState(self)
        self.cur_state.r[:] = r
        self.cur_state.v[:] = v

        self.static_f = np.zeros(3)  # Static forces for hack

        self.lock_in_place = False  # Don't move if true.

        self.avgKE = 0.0    # Running average of KE
        self.avgPE = 0.0    # Running average of PE

        self.charge = 0.0   # Defined by subclass for e and p
        self.mass = 0.0     # Defined in subclass for e and p
        self.symbol = 'e'   # Defined in subclass for e and p

    # def r(self):
    #     """ Current 3D Position Vector. """
    #     return self.cur_state.r

    def add_static_force(self):
        """ Add static forces to beginning and ending forces. """
        # self.fx += self.static_fx
        # self.fy += self.static_fy
        # self.fz += self.static_fz
        # self.end_fx += self.static_fx
        # self.end_fy += self.static_fy
        # self.end_fz += self.static_fz
        self.cur_state.f += self.static_f
        self.end_state.f += self.static_f

    def es_force(self, p: 'Particle'):
        # Electrostatic force between self and p per coulomb's law.
        # Force on self, caused by p.
        # real force, not R limit limited force
        # Returns 0,0,0 instead of infinity for two particles located
        # on same spot.

        if p is self:
            return np.zeros(3)

        # dx = (self.cur_state_r0 - p.cur_state_r0)
        # dy = (self.cur_state_r1 - p.cur_state_r1)
        # dz = (self.cur_state_r2 - p.cur_state_r2)
        dr: ndarray = self.cur_state.r - p.cur_state.r

        r2, l2 = self.distance_squared(p)

        if r2 == 0.0:
            return np.zeros(3)  # Bogus but prevents DBZ -- should be infinity

        force = CONST_KE * self.charge * p.charge / r2

        r = math.sqrt(r2)

        # return force * dx / r, force * dy / r, force * dz / r
        return dr * (force / r)

    def gravity_force(self, p: 'Particle'):
        # Magnitude of gravity between self and other particle
        r2, l2 = self.distance_squared(p)
        f = CONST_G * self.mass * p.mass / r2
        return f

    def magnetic_field(self, p: 'Particle'):
        """
            Returns a 3D field vector B = (x,y,z)
            Calculate the magnetic field created at self, by p.
            B = (1e-7 q1 v1 x r_hat) / r^2
            r_hat is the unit vector pointing from p to self
        """

        r2, l2 = self.distance_squared(p)

        if r2 == 0.0:
            return np.zeros(3)

        # print " distance is", r2, math.sqrt(r2)

        # r_hat = (self.cur_state_r0 - p.cur_state_r0, self.cur_state_r1 - p.cur_state_r1, self.cur_state_r2 - p.cur_state_r2)
        # r_hat = product(1.0 / math.sqrt(r2), r_hat)
        r_hat = (self.cur_state.r - p.cur_state.r) / math.sqrt(r2)

        # The books say based on current flow, this should be correct:
        # This follows the right hand rule for positive current flow
        # b_vec = product(1e-7 * p.charge / r2, cross(p.v(), r_hat))
        b_vec: ndarray = np.cross(p.cur_state.v, r_hat) * (
                1e-7 * p.charge / r2)

        # This is inverted
        # b_vec = product(-1.0, b_vec)  # backwards
        b_vec = b_vec * -1.0
        # 2021 -- I have no clue if I added this inversion per the note
        # below or if this inversion is what is needed per the books?

        # The correct way, means that motion causes electrons to push HARDER
        # apart from each other, and ep pairs to pull harder together.  So when
        # the velocity == c, the magnetic force equals the electrostatic force.
        # If they add, as the books imply they do, it means that going the
        # speed of light means the force doubles.  But that is fucking
        # pointless. It doesn't make the speed of light special in any sense
        # related to orbital mechanics as far as I can deduce. To make it
        # special, the speed of light has to be the point where these forces
        # cancel each other out!  Which makes me logically deduce that This
        # sign has to be backwards to work correctly. So I'm going to play with
        # making it backwards from how the books all say it should be to see
        # what happens.

        return b_vec

    def magnetic_force(self, p: 'Particle'):
        # Returns a 3D force vector (x,y,z)
        # Calculate the force created on self, by the magnetic
        # field generated by p.
        # B = (1e-7 q1 v1 x rHat) / r^2
        # F = q2 V2 x B
        # F = q2 V2 X (1e-7 q1 v1 x rHat) / r^2
        # rHat is the unit vector pointing from p to self

        b_vec = self.magnetic_field(p)

        # print "mag force"
        # print " B is", B

        # f_vec = product(self.charge, cross(self.v(), b_vec))
        f_vec: ndarray = self.charge * np.cross(self.cur_state.v, b_vec)

        return f_vec

    def calculate_end_velocity(self, dt):
        """
            Calculate end_state.v from cur_state v and f and
            end_state f.
        """
        # Assume linear change in acceleration from start (fx to end end_fx)
        # (yes this is the correct integral of a linear change in acceleration)
        # I had to recalculate it 10-4-2016 to verify

        avg_f = (self.cur_state.f + self.end_state.f) / 2.0
        dv = avg_f * dt / self.mass
        self.end_state.v = self.cur_state.v + dv

    def calculate_end_position(self, dt):
        """
            Calculate end position (r) using cur_state r, v and f, and
            end_state f.
        """

        # Assume force (acceleration) changes but is linear from start to end
        # x = 1/2 at² + vt + x
        # a is (2*as + ae)/3  --- where as is a start, and ae is a end.
        # end_x = 1/2 ((2as + ae)/3)t² + vt + x

        self.end_state.r = (self.cur_state.r + self.cur_state.v * dt +
                            (self.cur_state.f + 0.5 * self.end_state.f) /
                            (3.0 * self.mass) * dt ** 2)

    def move(self):
        """ Update cur_state from end_state. """
        # Make end state the current state
        # Save current state as old state
        # Leaves ending force the same as the starting

        if self.lock_in_place:
            # Hack to prevent this particle from moving.
            self.cur_state.v = np.zeros(3)
            self.cur_state.f = np.zeros(3)
            self.end_state.f = np.zeros(3)
            return

        self.cur_state = self.end_state.copy()

    def reset_state(self):
        # Reset so we can recompute step with new DT
        # Reset force end to match force begin
        # Everything else will be recomputed again.

        # self.end_fx = self.fx
        # self.end_fy = self.fy
        # self.end_fz = self.fz

        self.end_state.f[:] = self.cur_state.f

    def kinetic_energy(self):
        """ Kinetic energy of cur_state. """
        return self.kinetic_energy_calc(self.cur_state.v)

    def kinetic_end_energy(self):
        """ Kinetic energy of end_state """
        return self.kinetic_energy_calc(self.end_state.v)

    def kinetic_energy_calc(self, v: ndarray) -> float:
        """ 1/2 mv² """
        return 0.5 * self.mass * v.dot(v)

    def set_kinetic_energy(self, ke: float):
        """
            Update magnitude of velocity to match using given ke
            Keep direction of vector the same.
        """
        new_v2 = ke / (0.5 * self.mass)
        # old_v2 = (self.vx ** 2.0 + self.vy ** 2.0 + self.vz ** 2.0)
        # old_v2 = np.sum(self.v() ** 2)
        old_v2 = np.sum(self.cur_state.v ** 2)
        # print "in set kinetic new_v2 is", new_v2
        # print "in set kinetic old_v2 is", old_v2
        new_v = math.sqrt(new_v2)
        old_v = math.sqrt(old_v2)
        # # self.vx *= new_v2 / old_v2
        # # self.vy *= new_v2 / old_v2
        # # self.vz *= new_v2 / old_v2
        # self.vx *= new_v / old_v
        # self.vy *= new_v / old_v
        # self.vz *= new_v / old_v
        self.cur_state.v *= new_v / old_v

    def distance_squared(self, p: 'Particle'):
        """ Distance**2 between self and p.
            :returns: distance**2, limited_distance**2
        """

        if p is self:
            return self.limited_distance2(0.0)

        d_squared: float = np.sum((self.cur_state.r - p.cur_state.r) ** 2)

        return self.limited_distance2(d_squared)

    def end_distance_squared(self, p: 'Particle'):  # distance squared
        """ Distance**2 between self and p end states.
            :returns: distance**2, limited_distance**2
        """
        if p is self:
            return self.limited_distance2(0.0)

        # dx = (self.end_x - p.end_x)
        # dy = (self.end_y - p.end_y)
        # dz = (self.end_z - p.end_z)

        # d2 = dx ** 2.0 + dy ** 2.0 + dz ** 2.0

        d_squared: float = np.sum((self.end_state.r - p.end_state.r) ** 2)

        return self.limited_distance2(d_squared)

    @staticmethod
    def limited_distance2(d2: float):
        """ Limit distance to RLimit to solve computational problems.
            return (real, limited) tuple
        """
        return d2, max(d2, RLimit ** 2)

    def distance(self, p):
        r, r_limited = self.distance_squared(p)
        return math.sqrt(r), math.sqrt(r_limited)

    def end_distance(self, p):
        r, r_limited = self.end_distance_squared(p)
        return math.sqrt(r), math.sqrt(r_limited)

    def potential_energy(self, p):
        # potential energy between self and particle P

        if p is self:
            return 0.0  # No potential energy for self

        r, r_limited = self.distance(p)

        return self.potential_energy_for_distance(p, r)

    def potential_end_energy(self, p):
        # potential energy between self and particle P

        if p is self:
            return 0.0  # Bogus should be +infinity

        r, r_limited = self.end_distance(p)

        return self.potential_energy_for_distance(p, r)

    def potential_energy_for_distance(self, p, d):

        # if d == 0:
        # return 0.0	# Bogus should be +infinity

        if d >= RLimit:
            return CONST_KE * self.charge * p.charge / d

        self.sim.inside_r_limit_count += 1

        x = CONST_KE * self.charge * p.charge / RLimit

        return x + (x / RLimit) * (RLimit - d)


class Electron(Particle):
    """ A Single Electron. """
    def __init__(self, sim: 'Simulation',
                 r: Tuple[float, float, float],
                 v=(0.0, 0.0, 0.0)):
        """ A Single Electron in the simulation.

        Args:
            r: 3D location vector (x,y,z) in meters.
            v: 3D velocity vector in m/s
        """
        super().__init__(sim, r, v)
        self.charge = CONST_E_CHARGE
        self.mass = CONST_E_MASS
        self.symbol = 'e'


class Proton(Particle):
    """ A Single Proton. """
    def __init__(self, sim: 'Simulation',
                 r: Tuple[float, float, float],
                 v=(0.0, 0.0, 0.0),
                 n=1.0):
        """ A single Proton

        Args:
            r: 3D position vector in meters
            v: 3D velocity vector in m / s
            n: Proton has n times the normal charge and mass
        """
        super().__init__(sim, r, v)
        self.charge = n * CONST_P_CHARGE  # in Coulombs
        self.mass = n * CONST_P_MASS
        self.symbol = 'P'


class Neutron(Particle):
    """ A Single Fake Neutron.
        Mass of a Proton.  Charge of an Electron.
        Only interacts with Protons.
        Testing the idea of how this will "glue" together
        protons to form an atomic nucleus.
    """
    def __init__(self, sim: 'Simulation',
                 r: Tuple[float, float, float],
                 v=(0.0, 0.0, 0.0),
                 n=1.0,
                 nc=1.0,
                 ):
        """ A single Neutron

        Args:
            r: 3D position vector in meters
            v: 3D velocity vector in m / s
            n: Neutron has n times the normal charge and mass
            nc: n times the normal charge
        """
        super().__init__(sim, r, v)
        self.charge = nc * n * CONST_E_CHARGE  # in Coulombs
        self.mass = n * CONST_P_MASS
        self.symbol = 'N'


class Simulation:
    """ A 3D simulation of electrons and protons interacting.
        Includes a graphic display animation.
    """

    def __init__(self,
                 title='Physics-Sim',
                 screen_width=600,
                 screen_height=400,
                 screen_depth=500,
                 pixels_per_angstrom=200,
                 dt_min=1e-30,
                 dt_max=1.0,
                 dt_adjust=True,
                 stop_at=0.0,
                 do_v_force=True,
                 total_force: Callable[[ParticleState, ParticleState],
                                       ndarray] = None,
                 ):
        """
        3D Simulation of Protons and Electrons.

        Args:
            title: title for windows
            screen_width: pixels - window width
            screen_height: pixels - window height
            screen_depth: pixels - z dimension of simulation box
            pixels_per_angstrom: conversion factor- simulated space to screen
            dt_min: seconds - Smallest Time Step of simulation
            dt_max: seconds - Largest Time Stop of seconds
            dt_adjust: boolean - Adjust time step automatically
            stop_at: seconds - simulation stops when it goes past this.
                     0.0 means run forever.
            do_v_force: boolean - include v_force calculations
            total_force: A function to calculate force between to particles
        """
        self.title = title      # Screen Title
        self.screen_width = screen_width
        self.screen_height = screen_height
        self.screen_depth = screen_depth
        self.pixels_per_angstrom = pixels_per_angstrom
        self.dt_min = dt_min
        self.dt_max = dt_max
        self.dt_adjust = dt_adjust
        self.stop_at = stop_at      # Simulations time to stop
        self.do_v_force = do_v_force
        self.total_force = total_force

        self.world: list[Particle] = []

        self.fps_limit = 500    # pygame thing - Frames Per Second
        self.fps_avg = 1.0      # Computed average speed of simulation
        self.run_me = True
        self.cycle_count = 0
        self.now = 0.0          # Simulation time in seconds
        self.inside_r_limit_count = 0
        self.e_bounce_count = 0
        self.p_bounce_count = 0

        self.p_pair_last_distance: Dict[
            Tuple[float, float], Tuple[float, float, float]] = {}

        self.total_ke = 0.0
        self.total_pe = 0.0

        self.last_ke = self.total_ke
        self.last_pe = self.total_pe

        self.starting_total_ke = self.total_ke
        self.starting_total_pe = self.total_pe

        self.energy_diff_last = 0.0
        self.energy_diff_max = 0.0
        self.max_vc = 0.0  # Max V in units of CONST_C

        self.dt = self.dt_min

        self.screen = pygame.display.set_mode(
            (self.screen_width, self.screen_height))
        self.clock = pygame.time.Clock()
        pygame.display.set_caption(self.title)

    def add_e(self,
              r: Tuple[float, float, float],
              v=(0.0, 0.0, 0.0),
              ):
        """ Add a single Electron to the simulation.

        Args:
            r: 3D position vector in meters
            v: 3D velocity vector in m / s
        Returns:
            Electron -  the electron added
        """
        e = Electron(self, r, v)
        self.world.append(e)
        return e

    def add_p(self,
              r: Tuple[float, float, float],
              v=(0.0, 0.0, 0.0),
              n=1.0):
        """ Add a single Proton to the simulation.

        Args:
            r: 3D position vector in meters
            v: 3D velocity vector in m / s
            n: Proton has n times the normal charge and mass
        Returns:
            Proton: the Proton added
        """
        p = Proton(self, r, v, n)
        self.world.append(p)
        return p

    def add_n(self,
              r: Tuple[float, float, float],
              v=(0.0, 0.0, 0.0),
              n=1.0,
              nc=1.0):
        """ Add a single Proton to the simulation.

        Args:
            r: 3D position vector in meters
            v: 3D velocity vector in m / s
            n: Proton has n times the normal charge and mass
            nc: n times the normal charge
        Returns:
            Neutron: the Proton added
        """
        n = Neutron(self, r, v=v, n=n, nc=nc)
        self.world.append(n)
        return n

    def add_e_a(self,
                r: Tuple[float, float, float],
                v=(0.0, 0.0, 0.0),
                ):
        """ Add a single Electron to the simulation in ANGSTROM units.

        Args:
            r: 3D position vector in Angstroms
            v: 3D velocity vector in m / s
        Returns:
            Electron: The electron added
        """
        e = self.add_e(np.array(r) * ANGSTROM, v)
        return e

    def add_p_a(self,
                r: Tuple[float, float, float],
                v=(0.0, 0.0, 0.0),
                n=1.0):
        """ Add a single Proton to the simulation in ANGSTROM units.

        Args:
            r: 3D position vector in Angstroms
            v: 3D velocity vector in m / s
            n: Proton has n times the normal charge and mass
        Returns:
            Proton: The Proton added
        """
        p = self.add_p(np.array(r) * ANGSTROM, v, n)
        return p

    def add_n_a(self,
                r: Tuple[float, float, float],
                v=(0.0, 0.0, 0.0),
                n=1.0,
                nc=1.0):
        """ Add a single Neutron to the simulation in ANGSTROM units.

        Args:
            r: 3D position vector in Angstroms
            v: 3D velocity vector in m / s
            n: Neutron has n times the normal charge and mass
            nc: n times the normal charge
        Returns:
            Neutron: The Neutron added
        """
        n = self.add_n(np.array(r) * ANGSTROM, v=v, n=n, nc=nc)
        return n

    def add_ep_a(self, center: Tuple[float, float, float],
                 radius: float = None, velocity: float = None):
        """ Add an orbiting e/p pair to the Simulation.  A Hydrogen Atom.

        The location of the center of the mass of the orbiting pair is
        specified as 'center'. The Orbit can be specified as the radius of the electron
        orbit about the center of mass or by the velocity of the electron.
        If no radius is defined, .5Å is used making the atom have a diameter
        of 1Å which is in the rough ballpark of a real Hydrogen atom.

        The electron is placed to the right (+x) of the proton creating an orbit
        on the x-y plane.

        Args:
            center: The 3D location of the center of mass of the pair.
            radius: radius of electron orbit in Angstroms.
            velocity: velocity of electron in m/s
        Returns:
            Tuple[Electron, Proton]: The ep pair created
        """
        # j = (1+me/mp)
        # d = rj
        # c force = kq²/r²j²
        # When these two are equal, we have a circular orbit.
        # kq²/r²j² = mv²/r
        # solving for r:
        # kq²/j² = mv²r
        # r = kq²/mv²j²
        # solving for v:
        # v² = kq²/mrj²
        # v = sqrt(kq²/mrj²)
        cm_r = np.array(center) * ANGSTROM     # Center of mass
        rad = radius
        v = velocity
        j = 1 + CONST_E_MASS/CONST_P_MASS
        if rad is None and v is None:
            rad = 0.2           # Angstroms - Default
        if rad is not None:
            rad *= ANGSTROM     # Convert Angstrom to meters

        if v is None and rad is not None:
            # Calculate v from rad
            # v = sqrt(kq²/mrj²)
            v = math.sqrt(CONST_KE * CONST_P_CHARGE ** 2 /
                          (CONST_E_MASS * rad * (j**2)))

        if rad is None and v is not None:
            # Calculate rad from v
            # r = kq²/mv²j²
            rad = CONST_KE * (CONST_P_CHARGE ** 2) / \
                  (CONST_E_MASS * (v**2) * (j**2))

        # print()
        # print("rad is", rad, "v is", v)
        # print("rad is", rad/ANGSTROM, "Å")

        p_r = cm_r - np.array((rad * CONST_E_MASS/CONST_P_MASS, 0.0, 0.0))
        p_v = np.array((0.0, -v * CONST_E_MASS/CONST_P_MASS, 0.0))
        p = self.add_p(p_r, v=p_v)

        e_r = cm_r + np.array((rad, 0.0, 0.0))
        e_v = np.array((0.0, v, 0.0))
        e = self.add_e(e_r, v=e_v)

        return e, p

    def init_world(self):
        self.center_mass()
        self.zero_momentum()
        self.init_forces()

        self.total_ke = self.total_kinetic_energy()
        self.total_pe = self.total_potential_energy()

        self.starting_total_ke = self.total_ke
        self.starting_total_pe = self.total_pe

        self.last_ke = self.total_ke
        self.last_pe = self.total_pe

    def center_mass(self):
        """ Move world center of mass to center of screen. """
        center_m = self.center_of_mass()

        # Center of screen (0,0,0 is upper left front corner)
        center = (
            (self.screen_width / self.pixels_per_angstrom * ANGSTROM) / 2.0,
            (self.screen_height / self.pixels_per_angstrom * ANGSTROM) / 2.0,
            (self.screen_depth / self.pixels_per_angstrom * ANGSTROM) / 2.0)

        dr = center-center_m    # Distance vector to move everything

        for p in self.world:
            p.cur_state.r += dr

    def center_of_mass(self) -> ndarray:
        """ return center of mass of the world as pos vector. """
        c = np.zeros(3)
        total_mass = 0.0
        for p in self.world:
            c += p.cur_state.r * p.mass
            total_mass += p.mass
        return c / total_mass

    def zero_momentum(self):
        """
            Make cur_state momentum zero.
            Adjust all velocities by a constant to make total momentum zero.
        """
        tm = self.total_momentum()
        total_mass = sum([p.mass for p in self.world])
        dv = -tm / total_mass
        for p in self.world:
            p.cur_state.v += dv

    def total_momentum(self):
        """ Compute Momentum vector for all particles in world[]. """
        s = np.zeros(3)
        for p in self.world:
            s += p.cur_state.momentum()
        return s

    def init_forces(self):
        """
            Set up cur_state and end_state to start simulation.
            Cur_state position and velocity are the starting state.
        """
        # Copy to cur_state to end_state for all p.
        for p1 in self.world:
            p1.end_state = p1.cur_state.copy()
        # Compute the end forces using the end_state.
        self.compute_end_forces()
        # end_state now has valid r, p, and v.
        # Copy it all back to cur_state
        self.move_all()
        # cur_state and end_state are the same now
        # and both include valid forces.

    def compute_end_forces(self):
        """ Update end_state forces using end_state position and velocity. """

        for p in self.world:
            p.end_state.zero_force()

        for i1 in range(len(self.world)):
            p1 = self.world[i1]
            for i2 in range(i1+1, len(self.world)):
                p2 = self.world[i2]
                if self.total_force is not None:
                    f = self.total_force(p1.end_state, p2.end_state)
                else:
                    f = p1.end_state.total_force(p2.end_state)
                p1.end_state.f += f
                p2.end_state.f -= f
            # Add constant static force (experimental hack)
            p1.end_state.f += p1.static_f

    def total_kinetic_energy(self):
        total_ke = 0.0

        for p in self.world:
            total_ke += p.kinetic_energy()

        return total_ke

    def total_potential_energy(self):
        total_pe = 0.0

        for i in range(len(self.world)):
            for j in range(i + 1, len(self.world)):
                total_pe += self.world[i].potential_energy(self.world[j])

        return total_pe

    def total_end_kinetic_energy(self):
        total_ke = 0.0

        for p in self.world:
            total_ke += p.kinetic_end_energy()

        return total_ke

    def total_end_potential_energy(self):
        total_pe = 0.0

        for i in range(len(self.world)):
            for j in range(i + 1, len(self.world)):
                total_pe += self.world[i].potential_end_energy(self.world[j])

        return total_pe

    def move_all(self):
        """ Copy end_state to cur_state for all particles. """
        for p in self.world:
            p.move()

    def run(self):
        """ Run the simulation. """
        self.init_world()

        self.run_me = True
        loop_cnt = 0
        last_time = time.time()
        paused = False

        while self.run_me:
            loop_cnt += 1
            self.clock.tick(self.fps_limit)
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    # The run_me flag is not needed but I use it
                    # to prevent false "code unreachable" warnings.
                    self.run_me = False
                    break

            if paused:
                time.sleep(0.1)     # Maybe not needed?
                continue

            now = time.time()
            loop_time = now - last_time
            fps = 1 / loop_time
            self.fps_avg += (fps - self.fps_avg) * .001
            last_time = now

            if loop_cnt % 10 == 0:
                self.bounce_particles()
                self.draw_world()

            self.total_ke = self.total_kinetic_energy()
            self.total_pe = self.total_potential_energy()

            if loop_cnt % 30 == 0:
                self.print_stats()

                if 0.0 < self.stop_at < self.now:
                    paused = True
                    continue

            #####################################################################
            # Calculate next position now
            #####################################################################

            energy_diff_start, total_ke2 = self.calculate_end_state()

            #####################################################################
            # Make the move -- move current ending state to be the current
            # Also leaves ending force set to same as starting force
            #####################################################################

            for p in self.world:
                p.move()

            self.now += self.dt

            if Energy_fix:  # not compatible with v_force() - breaks PE
                # Fix total energy error
                # Adjust all velocities to remove energy error of energy_diff_start
                if energy_diff_start > total_ke2:
                    print("energy error greater than total KE error:",
                          energy_diff_start,
                          "total", total_ke2)
                    energy_diff_start = total_ke2 * 0.5
                    # Only try to take half of total KE out of the system
                    # time.sleep(1)
                    # sys.exit(1)

                for p in self.world:
                    # print
                    # print "Adjust KE for particle", i, "energy_diff_start error is", energy_diff_start
                    ke = p.kinetic_energy()
                    # print "Current ke is", ke
                    # print "percent of of total is", ke/total_ke2*100.0, "%"
                    # print "amount added of total energy_diff_start is", ke/total_ke2*energy_diff_start
                    new_ke = ke - ke / total_ke2 * energy_diff_start
                    # print "new should be", new_ke
                    p.set_kinetic_energy(new_ke)
                    # print "new KE is now", p1.kinetic_energy()

        pygame.quit()

    def draw_world(self):
        """ Draw all the particles on the screen. """
        self.screen.fill(WHITE)

        # Sort world by z dimension so we draw far to near.
        # When z is equal, draw Proton first, so the larger
        # Protons don't cover up the smaller electrons.

        def sort_key(p: Particle):
            z = p.cur_state.r[2]
            sort_order = 0
            if isinstance(p, Electron):
                sort_order = 1
            return z, sort_order
        s_list = sorted(self.world, key=sort_key)
        for p2 in s_list:
            self.draw_particle(p2)
        pygame.display.flip()
        return

    def draw_particle(self, p: Particle):
        """ Draw a single particle on the screen. """
        x = self.space_to_pixels(p.cur_state.r[0])
        y = self.space_to_pixels(p.cur_state.r[1])
        z = self.space_to_pixels(p.cur_state.r[2])

        if isinstance(p, Electron):
            color = BLACK
            size = 2
        elif isinstance(p, Proton):
            color = RED
            size = 4
        else:   # A Neutron
            color = BLUE
            size = 4

        min_scale = 1.0
        max_scale = 3.0
        size *= (z / float(self.screen_depth)) * (max_scale - min_scale) + min_scale
        size = abs(int(size))

        pygame.draw.circle(self.screen, color, (x, y), size, 0)

    def print_stats(self):
        """ Print statistics for this loop. """

        crt.clear_and_home()

        a_width = self.screen_width / self.pixels_per_angstrom
        a_height = self.screen_height / self.pixels_per_angstrom
        a_depth = self.screen_depth / self.pixels_per_angstrom
        print(f"PixelsPerAngstrom", self.pixels_per_angstrom, end='')
        print(f"  World ({a_width}A x {a_height:}A x {a_depth:}A)", end='')
        print(" ", self.title)

        print(f"Sim Time:{self.now * 1000_000_000:.20f} ns", end='')
        print("  DT:%4.1e" % self.dt, end='')
        print(f"  FPS:{self.fps_avg:.1f}")

        if self.stop_at > 0.0:
            print(f"    Stop at {self.stop_at * 1000_000_000:.20f} ns")

        print("do_v_force:", self.do_v_force)
        print()

        re = ((self.total_ke - self.starting_total_ke) + (
                self.total_pe - self.starting_total_pe))
        print("Total ke     %8.1e" % self.total_ke,
              "change from last %8.1e" % (self.total_ke - self.last_ke))
        print("Total pe     %8.1e" % self.total_pe,
              "change from last %8.1e" % (self.total_pe - self.last_pe))
        print("Total    e:  %8.1e" % (self.total_ke + self.total_pe))
        print("Relative e:  %8.1e" % (re,), end=' ')
        if self.total_ke:
            print("  %% of Total ke:%8.4f %% (change from start)" % (
                abs(re * 100.0 / self.total_ke)))
        else:
            print()

        print(f"Move e err:  {self.energy_diff_last:8.1e}", end='')
        print(f" max move error   {self.energy_diff_max:8.1e}")

        p_vec = self.total_momentum()
        print(f"Momentum err:{p_vec[0]:8.1e} {p_vec[1]:8.1e} "
              f"{p_vec[2]:8.1e}", end='')
        print(f"   mag: {magnitude(p_vec):.1e}")

        print()

        print("Inside R limit:", self.inside_r_limit_count,
              "  P Bounce:", self.p_bounce_count,
              "  E Bounce:", self.e_bounce_count)

        self.last_ke = self.total_ke
        self.last_pe = self.total_pe

        self.print_proton_distance()

        print()

        self.print_particle_stats()

        print(f"                   Max Velocity: {self.max_vc:7.5f}c")

    def print_proton_distance(self):
        """
            Print distance between Protons.
            Only first 10 shown.
        """
        cnt = 0
        for i in range(len(self.world)):
            p1 = self.world[i]
            if not isinstance(p1, Proton):
                continue
            for j in range(i + 1, len(self.world)):
                p2 = self.world[j]
                if not isinstance(p2, Proton):
                    continue
                d = p1.distance(p2)[0]
                pair = (i, j)
                last_d, last_time, avg_v = self.p_pair_last_distance.get(pair, (
                    d, self.now - self.dt, 0.0))
                v = (d - last_d) / (
                        self.now - last_time)  # positive v if moving away
                avg_v += (v - avg_v) * .01
                self.p_pair_last_distance[pair] = (
                    d, self.now, avg_v)  # save for next time
                print("Distance between", end=' ')
                print(f"P{i:02d} P{j:02d} {d / ANGSTROM:11.6f} A", end='')
                print(f"   v:{v:+.2e}", end='')
                print(f"   v:{avg_v:+.2e} m/s", end='')
                print()
                cnt += 1
                if cnt > 10:
                    break
            if cnt > 10:
                break

    def print_particle_stats(self):
        total_avg_ke = 0.0
        total_momentum_mag = 0.0

        for p in self.world:
            ke = p.kinetic_energy()
            p.avgKE += (ke - p.avgKE) * 0.0001
            total_avg_ke += p.avgKE
            total_momentum_mag += magnitude(p.cur_state.momentum())

        for i in range(len(self.world)):
            p = self.world[i]
            print(f"{p.symbol}{i:02}", end='')
            pv = p.cur_state.v
            print(f" vx:{pv[0]:10.2e}  vy:{pv[1]:10.2e}", end='')
            vc = magnitude(pv) / CONST_C
            self.max_vc = max(self.max_vc, vc)
            print(f" {vc:0.5f}c", end='')
            print(" x:%6.2f A" % (p.cur_state.r[0] / ANGSTROM), end=' ')
            # print " KE:%10.2e" % p1.avgKE
            # if total_avg_ke:
            #     print(" KE:%6.2f%%" % (p1.avgKE * 100 / total_avg_ke))
            # else:
            #     print(" KE:?????")
            if self.total_ke:
                print(" KE:%5.1f%%" % (
                        p.kinetic_energy() * 100 / self.total_ke), end='')
            else:
                print(" KE:?", end='')
            momentum = magnitude(p.cur_state.momentum())
            print(f"  p:{momentum:.2e}", end='')
            if total_momentum_mag:
                print(f" {momentum * 100 / total_momentum_mag:5.1f}%", end='')

            # Frequency of electron in orbit.
            if isinstance(p, Electron) or isinstance(p, Proton):
                f = p.cur_state.f_from_v()
                if f / 1e15 > 1000.0:
                    print(f" f:{f / 1e18:.1e} EHz", end='')
                else:
                    print(f" f:{f / 1e15:.1e} PHz", end='')
                print(f"  wl:{CONST_C * 1e9 / f:.3e} nm",
                      end='')  # wave length

            print()

    def calculate_end_state(self):
        """
            Calculate end_state from cur_state.
            Assumes end_state forces are a copy of cur_state forces.
            Adjusts self.dt as needed based on error estimates.
        """
        while True:  # DT restart loop

            # At this point, all particles have a known position, velocity, and
            # force and end_state force estimates ending force.  At the start
            # end_state force is just the same as the cur_state force.

            for it in range(3):
                # print "begin iteration", it
                for p in self.world:
                    p.calculate_end_velocity(self.dt)
                    p.calculate_end_position(self.dt)

                # Update ending force based on ending position calculated above.
                # Ending position is a function of ending force so when we update
                # the force, we can then loop back and calculate a new ending
                # Position

                # Could calculate PE in same loop we use for updating forces
                # because it duplicates a lot of the same work, but need to
                # restructure code to do that.

                self.compute_end_forces()

                # total_ke2 = total_end_kinetic_energy(self.world)
                # total_pe2 = total_end_potential_energy(self.world)
                # curEnergy = total_ke2 + total_pe2
                # error = (self.total_ke - total_ke2) + (self.total_pe - total_pe2)
                # print "Energy error after iteration     %18.10e" % error

            if Energy_fix2:
                # Ok, one last time, fudge velocity based on current position and
                # total force for this position This is just always good I think.
                # Probably shouldn't be an option.

                for p in self.world:
                    p.calculate_end_velocity(self.dt)

            # Now calculate total Energy after this move.

            total_ke2 = self.total_end_kinetic_energy()
            total_pe2 = self.total_end_potential_energy()

            # energy_diff = (total_ke2 + total_pe2) - (self.total_ke + self.total_pe) # last move
            # error only

            # Ah, a computational problem showed up.  When one of the two is many
            # orders of magnitude different from the other, the accuracy is all
            # lost in the difference if done the above way vs the way below!

            self.energy_diff_last = (total_ke2 - self.total_ke) + (
                    total_pe2 - self.total_pe)  # last move error only
            energy_diff_start = total_ke2 - self.starting_total_ke  # Error from start
            energy_diff_start += total_pe2 - self.starting_total_pe

            self.cycle_count += 1

            ######################################################################
            # Dynamically change dt to maintain error and maximise simulation speed
            ######################################################################

            if self.dt_adjust and total_ke2 != 0.0:
                # print "==DO DT ADJUST self.cycleCount", self.cycleCount, "abs(energy_diff)", abs(energy_diff),
                # print "total_ke2", total_ke2, "percent", abs(energy_diff) / total_ke2
                if self.dt < self.dt_max and self.cycle_count > 3 and abs(
                        self.energy_diff_last) / total_ke2 < 0.0001:
                    # print "SPEEDUP -- increase DT abs(diff)/total is", abs(energy_diff) / total_ke2
                    self.dt *= 2.0
                    self.dt = min(self.dt, self.dt_max)
                    # self.cycleCount = 0
                    # continue
                    # No need to restart -- take this move as fine but use a larger dt
                    # for the next move

                elif self.dt > self.dt_min and abs(
                        self.energy_diff_last) / total_ke2 > 0.001:
                    # print "SLOWDOWN -- reduce DT abs(diff)/total is", abs(energy_diff) / total_ke2
                    self.dt /= 2.0
                    self.dt = max(self.dt, self.dt_min)
                    for p1 in self.world:
                        p1.reset_state()
                    self.cycle_count = 0
                    continue

            self.energy_diff_max = max(self.energy_diff_max,
                                       abs(self.energy_diff_last))

            break  # End of DT adjust loop

        return energy_diff_start, total_ke2

    def bounce_particles(self):
        """
            Check for particles hitting the wall and make them bounce.
            Updates world momentum, and total_ke and total_pe if needed.
        """
        bounce_happened = False
        for p in self.world:
            bounce_happened |= self.bounce_particle(p)

        if bounce_happened:
            # Bounce happened one or more times for world.

            # We re-zero the momentum of the world on a bounce. The run starts
            # with a zero total momentum and a bounce messes that up. So this
            # sets back to to zero so we can still look at the total momentum
            # as an measure of error. This might be a messy problems if two or
            # more particles go off screen at same time and then we change the
            # velocity which could, under just the right conditions, cause a
            # particle that was heading back from a bounce to turn around and
            # head off screen. Could be a disaster for a large sim with lots of
            # fast moving particles bouncing all the time.

            self.zero_momentum()

            # We move particles back into the window on a bounce and that
            # changes the PE. And if we use ke_change_factor that changes the
            # KE. So we must reset energy totals and reset what we "started" at
            # so the bounce doesn't make it look like we have a large
            # simulation error (energy change from start of run).

            self.total_ke = self.total_kinetic_energy()
            self.total_pe = self.total_potential_energy()

            self.starting_total_ke = self.total_ke
            self.starting_total_pe = self.total_pe
            # TODO wait, do I need to fix end_state and forces???

        return

    def bounce_particle(self, p: Particle):
        """
            Check for particle off screen.
            Process bounce as needed.
            return True if there was a bounce.
        """
        x = self.space_to_pixels(p.cur_state.r[0])
        y = self.space_to_pixels(p.cur_state.r[1])
        z = self.space_to_pixels(p.cur_state.r[2])

        vx = p.cur_state.v[0]
        vy = p.cur_state.v[1]
        vz = p.cur_state.v[2]

        inset = 10  # Bounce at 10 pixels from edge of screen
        if isinstance(p, Proton):
            inset = 40

        # ke_change_factor = 0.25   # Remove energy on bounce (reduce to 25%)
        ke_change_factor = 1.0      # Don't remove energy on bounce

        bounce = False

        if x < inset and vx < 0:
            bounce = True
            vx *= -1
            p.cur_state.r[0] = self.pixels_to_space(inset)

        if x > self.screen_width - inset and vx > 0:
            bounce = True
            vx *= -1
            p.cur_state.r[0] = self.pixels_to_space(self.screen_width - inset)

        if y < inset and vy < 0:
            bounce = True
            vy *= -1
            p.cur_state.r[1] = self.pixels_to_space(inset)

        if y > self.screen_height - inset and vy > 0:
            bounce = True
            vy *= -1
            p.cur_state.r[1] = self.pixels_to_space(self.screen_height - inset)

        if z < inset and vz < 0:
            bounce = True
            vz *= -1
            p.cur_state.r[2] = self.pixels_to_space(inset)

        if z > self.screen_depth - inset and vz > 0:
            bounce = True
            vz *= -1
            p.cur_state.r[2] = self.pixels_to_space(self.screen_depth - inset)

        if bounce:
            # Process bounce for this particle
            # Update changes to velocity if needed
            p.cur_state.v[0] = vx
            p.cur_state.v[1] = vy
            p.cur_state.v[2] = vz

            if ke_change_factor != 1.0:
                # Adjust KE of the particle if we are using ke_change_factor to
                # reduce KE on bounce (slow particle down)
                p.set_kinetic_energy(p.kinetic_energy() * ke_change_factor)

            # Update statistics about bounces.
            if isinstance(p, Proton):
                self.p_bounce_count += 1
            else:
                self.e_bounce_count += 1

        return bounce

    def space_to_pixels(self, space):
        # 0,0 is the same in both and is the top left corner of the screen
        return int(self.pixels_per_angstrom * space / ANGSTROM)

    def pixels_to_space(self, pixels):
        return pixels * ANGSTROM / self.pixels_per_angstrom


def unit_test():
    """ Nothing here yet. """
    pass


if __name__ == '__main__':
    unit_test()
