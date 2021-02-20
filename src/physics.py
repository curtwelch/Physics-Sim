#!/usr/bin/python

"""
    physics.py - Atomic particle simulation experiments

    2016-05-18 Started this file.  Took it from xnet.py from my AI work.
    2021-02-05 Created PyCharm Project, Put on Github, converted
                to python3.9, and started major cleanup of this old code.

"""

# import random
import time
import sys
import math
# import heapq
import numpy as np
from numpy import ndarray
# import matplotlib.pyplot as plt
# from operator import attrgetter
from typing import List, Tuple, Dict
# import timeit

# import os
import pygame

BLACK = 0, 0, 0
WHITE = 255, 255, 255
RED = 255, 0, 0

Angstrom = 1.0e-10  # One Angstrom 10e-10 meters
CONST_C = 299792458.0  # Speed of light m/s - defined constant
# CONST_KE = 8.9875517873681764e9  # Coulomb's constant (1/4 pi e)  K(sub)e
# CONST_KE = self.c * self.c * 1.0e-7  # OLD?
CONST_KE = 8.9875517923e9  # Coulomb's constant; New? 8.9875517923(14)×10^9

# RLimit = 0.1 * Angstrom			# Radius limit hack
# RLimit = 0.0001 * Angstrom		# Radius limit hack
RLimit = 0.0000001 * Angstrom  # Radius limit hack

InsideRLimitCount = 0
eBounceCount = 0
pBounceCount = 0

doMagnetic = True
doMagneticInverse = False

dtMin = 1e-30
# dtMin = 1e-30

# dtMax = 1e-10
# dtMax = 1e-1
dtMax = 1e-18

dtAdjust = True  # Auto adjust dt time step based on error

energyFix = False  # fix based on total PE+KE at start
energyFix2 = False

# screen_size = screen_width, screen_height = 1000, 800
screen_size = screen_width, screen_height = 600, 400
screen_depth = 1000  # z dimension
pixelsPerAngstrom = 200.0

Stop_at = 0.000001 * 1e-9   # Pause sim after this clock time
# Stop_at = 0.0               # Run forever


def main():
    sol_test()


def mag_circle_test():
    # Copy forceCircleTest to crete mag_circle_test()
    # What is the mag field at different points?
    # Using Biot-Savant law the mag field is
    # dB = u0 I Dl x rHat / 4Pi R^2
    # I must integrate this all around the circle
    # for different points.
    # I'll just sum it up over the small steps
    # for different points.

    # Circle of radius 1. Using fake force as 1/d^2.
    # center at origin.
    # Field always points in same z direction.
    # So I just sum all the parts.

    y = 0.0
    for i in range(11):
        x = i * 0.1
        print("For x,y", x, y, end=' ')
        b_total = 0.0
        steps = 4000
        # dl = 2.0 * math.pi / steps
        last_px = 0.0  # First point
        last_py = 1.0
        # 2021 fixed the double use of 'i' but don't know if I fixed
        # it correctly.  Could have just broken the code.
        for j in range(1, steps + 1):  # UGH i used twice - changed to j
            angle_i = 2.0 * math.pi * j / steps
            px = math.sin(angle_i)
            py = math.cos(angle_i)
            dl = math.sqrt((last_px - px) ** 2.0 + (last_py - py) ** 2.0)
            dx = px - x
            dy = py - y
            r2 = dx ** 2 + dy ** 2
            # r = math.sqrt(r2)
            # For 1 amp current
            # b = 1e-7 * dl / r  # fake for experimenting
            # b = 1e-7 * dl / 1.0  # Fake again
            b = 1e-7 * dl / r2  # Real
            b_total += b
            last_px = px
            last_py = py
        print("total mag field is", b_total)
    sys.exit(1)

    # Conclusion ...
    # Mag field is weakest in middle, stronger near edges 
    # But is somewhat the same till it gets real near the edge?
    # Switching from 1/r^2 to 1/r doesn't make it the same!
    # But switching to 1/1 does!  But then we are just calculating 2 pi r I!  So of course!
    # For x,y 0.0 0.0 total mag field is 6.28318466122e-07
    # For x,y 0.1 0.0 total mag field is 6.34665117294e-07
    # For x,y 0.2 0.0 total mag field is 6.5449840221e-07
    # For x,y 0.3 0.0 total mag field is 6.90459852881e-07
    # For x,y 0.4 0.0 total mag field is 7.47998173954e-07
    # For x,y 0.5 0.0 total mag field is 8.37757954829e-07
    # For x,y 0.6 0.0 total mag field is 9.81747603315e-07
    # For x,y 0.7 0.0 total mag field is 1.2319969924e-06
    # For x,y 0.8 0.0 total mag field is 1.74532907256e-06
    # For x,y 0.9 0.0 total mag field is 3.30693929538e-06
    # For x,y 1.0 0.0 total mag field is 4.18946069487e+22  infinity? 0/0?


def force_circle_test():
    # What is the total force inside a circle of particles?
    # Is it the same at all places inside the circle or is it
    # different based on location.  Just do some quick math
    # test to find out.  In 2 dimensions.

    # Circle of radius 1. Using fake force as 1/d^2.
    # center at origin.

    y = 0.0
    z = 0.0
    for ii in range(9):
        x = ii * 0.1
        print("For x,y", x, y, end=' ')
        fx = 0.0
        fy = 0.0
        fz = 0.0
        steps = 400
        for i in range(steps):
            for j in range(steps):
                angle_i = 2.0 * math.pi * i / steps
                angle_j = 2.0 * math.pi * j / steps
                px = math.sin(angle_i)
                py = math.cos(angle_i)
                pz = math.cos(angle_j)
                dx = px - x
                dy = py - y
                dz = pz - z
                d2 = dx ** 2 + dy ** 2 + dz ** 2
                r = math.sqrt(d2)
                # force = 1.0 / r
                force = 1.0 / d2
                # Had to add /px below to adjust for the ds factor
                # that our step pattern was wrong
                fx += force * dx / r * px
                fy += force * dy / r * px
                fz += force * dz / r * px
        print("total force is", fx, fy, fz)
    sys.exit(1)

    # Conclusion ...
    # Force at point inside circle only sums to a constant if the
    # force equation is 1/x  not 1/x^2.
    # Force sums to constant when we sum over entire sphere for 1/x^22 force.
    # This makes sense because the size of the circle is 2pi r, or 4pi r^2 for surface
    # The use of 1/x^2 in physics is consistent with conservation in 3D space.


# def fast_test():
#     # p1 = Proton(0.0, 00.0, 0.0)
#     p1 = Electron(0.0, 00.0, 0.0)
#     for i in range(60):
#         x = (i - 30) * Angstrom / 10
#         e1 = Electron(x, 0.0, 0.0)
#         print("i", i, "x", x, x / Angstrom, "A", end=' ')
#         print(p1.potential_energy(e1))
#
#     sys.exit(1)


# def neutron_gravity_test():
#     # Some experiment on the idea that a neutron is
#     # is really a e p pair in orbit.  What would
#     # be the attraction between two such systems?
#     # Does it match gravity in relative force?
#     # The answer was the attraction is not even
#     # 1/x^2 in force!  oops.  That theory got blown
#     # out of the water!  The overlapping fields does
#     # create a 1/x^2 field!
#     # UGH -- 5-11-2018 thoughts: but do two systems in
#     # orbit actually end up moving closer together or
#     # further apart due to distorted orbits and orbital
#     # interactions -- this test did not answer that question
#     # and if the effect approximates 1/x^2 then it could be
#     # mistaken as 1/^x^2.  More testing required
#
#     close_space = 0.0001 * Angstrom
#     p1 = Proton(0.0, 0.0, 0.0)
#     e1 = Electron(0.0 + close_space, 0.0, 0.0)
#
#     p2 = Proton(100.0, 0.0, 0.0)
#     e2 = Electron(100.0 + close_space, 0.0, 0.0)
#
#     gravity_force = 0.0
#     gravity_force += p2.gravity_force(p1)
#     gravity_force += p2.gravity_force(e1)
#     gravity_force += e2.gravity_force(p1)
#     gravity_force += e2.gravity_force(e1)
#
#     print("Gravity force between two is:", gravity_force)
#
#     p2.zero_force()
#     p2.add_force(p1)
#     print("em force p2 to p1", p2.fx)
#     em_force = p2.fx
#
#     p2.zero_force()
#     p2.add_force(e1)
#     print("em force p2 to e1", p2.fx)
#     em_force += p2.fx
#
#     print("em force p2 to e1 and p1", em_force)
#
#     sys.exit(1)


def magnetic_test():
    world = []

    # p1 = Electron(0.0, 0.0, 0.0)
    p1 = Proton(0.0, 0.0, 0.0)
    # p1.vx = c/2.0
    # p1.vy = c/2.0
    # p1.vy = c/2.0
    # p1.vy = 2.0e8
    # p1.vz = 12.0e8
    world.append(p1)

    # p2 = Electron(1.0 * Angstrom, 0.0, 0.0)
    p2 = Electron(1.0 * Angstrom, 1.0 * Angstrom, 0.0)
    p2.vy = CONST_C / 2.0
    # p2.vy = -c/4.0
    # p2.vy = 1.0e8
    # p2.vz = 6.0e8
    world.append(p2)

    # p3 = Electron(1.0 * Angstrom, 0.0, 0.0)
    p3 = Electron(1.0 * Angstrom, 1.0 * Angstrom, 0.0)
    p3.vx = -CONST_C / 2.0
    p3.vy = -CONST_C / 2.0
    # world.append(p3)

    # maxV = CONST_C / math.sqrt(3.0)
    # print "max v is", maxV, "max mag is", math.sqrt(3.0 * maxV**2.0), "c is", c
    # zz

    # if False:
    #     speed = c / 2.0
    #     speed = maxV * random.random()
    #     speed = c
    #     p1.vx = (random.random() * 2.0 - 1.0)
    #     p1.vy = (random.random() * 2.0 - 1.0)
    #     p1.vz = (random.random() * 2.0 - 1.0)
    #     p1.vx, p1.vy, p1.vz = p1.product(speed / magnitude(p1.V()), p1.V())
    #
    #     speed = c / 2.0
    #     speed = maxV * random.random()
    #     speed = c
    #     p2.vx = (random.random() * 2.0 - 1.0)
    #     p2.vy = (random.random() * 2.0 - 1.0)
    #     p2.vz = (random.random() * 2.0 - 1.0)
    #     p2.vx, p2.vy, p2.vz = p2.product(speed / magnitude(p2.V()), p2.V())

    # if False:
    #     # Make the magnitude of the difference equal to c
    #     rV = p1.subtract(p1.V(), p2.V())
    #     s = magnitude(rV)
    #     p1.vx, p1.vy, p1.vz = p1.product(c / s, p1.V())
    #     p2.vx, p2.vy, p2.vz = p2.product(c / s, p2.V())

    # if False:
    #     p1.vx = 1.0
    #     p1.vy = 2.0
    #     p1.vz = 3.4
    #     p2.vx = 100.0
    #     p2.vy = 200.0
    #     p2.vz = 300.4

    # if False:
    #     # Time step test -- move particles forward in time
    #     # Based on velocity, see what dr and ds is!
    #     p0 = world[0]
    #     rp0 = world[0].R()
    #     rp1 = world[1].R()
    #     rBefore = p0.subtract(rp1, rp0)
    #     esBefore = world[0].es_force(world[1])
    #
    #     dt = 1e-20
    #
    #     magnetic_test2(world, dt=dt)
    #
    #     for it in range(3):
    #         for i in range(len(world)):
    #             p1 = world[i]
    #             p1.calculate_end_velocity(dt)
    #             p1.calculate_end_position(dt)
    #
    #         for p1 in world:
    #             p1.zero_end_force()
    #             for p2 in world:
    #                 p1.add_end_force(p2)
    #
    #     for p1 in world:
    #         p1.move()
    #
    #     print()
    #     print("Move Done!")
    #     print()
    #
    #     rp0 = world[0].R()
    #     rp1 = world[1].R()
    #     rAfter = p0.subtract(rp1, rp0)
    #     rDiff = p0.subtract(rAfter, rBefore)
    #     esAfter = world[0].es_force(world[1])
    #     esDiff = p0.subtract(esAfter, esBefore)
    #
    #     print(" dt is   ", dt)
    #     print(" r before", rBefore)
    #     print(" r after", rAfter)
    #     print(" dr =  after-before = ", rDiff)
    #     print(" dr diff magnitude", magnitude(rDiff))
    #     print(" es before  ", esBefore)
    #     print(" es after   ", esAfter)
    #     print(" es diff    ", esDiff)
    #     print(" es diff magnitude", magnitude(esDiff))
    #     print()
    #
    #     magnetic_test2(world, dt=dt)
    #     sys.exit(1)

    # magnetic_test2(world)
    # exit(1)
    #
    # # zz
    #
    # print("ZERO Momentum!")
    # print()
    # zero_momentum(world)
    # magnetic_test2(world)
    #
    # sys.exit(1)
    # if False:
    #     p1.vx = - p1.vx

    # Change velocities in different ways and print results again...

    # if False:
    #     # Change frame of reference
    #     print("Change inertial Frame of reference randomly")
    #
    #     # Because my test particles are on the same x axis, this
    #     # vx component is totally ignored in the equations so
    #     # I can assign random values without effecting the total force!
    #     # FALSE -- not true.  Changing dx changes fy and fz but not
    #     # fx.  So it is important, unless the vy and vz are zero.
    #     dx = maxV * (random.random() * 2.0 - 1.0)
    #     p1.vx += dx
    #     # dx = maxV * (random.random() * 2.0 - 1.0)
    #     p2.vx += dx
    #
    #     dy = maxV * (random.random() * 2.0 - 1.0)
    #
    #     p1.vy += dy
    #     p2.vy += dy
    #
    #     dz = maxV * (random.random() * 2.0 - 1.0)
    #
    #     p1.vz += dz
    #     p2.vz += dz

    # if False:
    #     print("Change X velocities randomly")
    #     # turns out, that doing this, has no
    #     # effect on the X force when the particles
    #     # are on the x axis!
    #     dx = maxV * (random.random() * 2.0 - 1.0)
    #     p1.vx += dx
    #     dx = maxV * (random.random() * 2.0 - 1.0)
    #     p2.vx += dx

    # if False:
    #     print("Change Y velocities randomly")
    #     # turns out, that doing this, has no
    #     # effect on the Z force when the particles
    #     # are on the x axis!
    #     dy = maxV * (random.random() * 2.0 - 1.0)
    #     p1.vy += dy
    #     dy = maxV * (random.random() * 2.0 - 1.0)
    #     p2.vy += dy

    # if False:
    #     print("Change Z velocities randomly")
    #     # turns out, that doing this, has no
    #     # effect on the Y force when the particles
    #     # are on the x axis!
    #     dz = maxV * (random.random() * 2.0 - 1.0)
    #     p1.vz += dz
    #     dz = maxV * (random.random() * 2.0 - 1.0)
    #     p2.vz += dz

    # if False:
    #     print("Randomly Rotate relative vx,vy vector!")
    #     # Coded this by mistake.  Meant to code vy,vz rotation which
    #     # I did below.  When we rotate vx,vy  vector as I'm doing
    #     # here, all forces change and the total magnitude of the force
    #     # changes as well!
    #     rV = p1.subtract(p1.V(), p2.V())
    #     old_vx = rV[0]
    #     old_vy = rV[1]
    #     old_r = math.sqrt(rV[0] ** 2.0 + rV[1] ** 2.0)  # x^2 + y^2
    #     print("starting relative velocity is:", rV)
    #
    #     # New random x,y direction:
    #     dx = maxV * (random.random() * 2.0 - 1.0)
    #     dy = maxV * (random.random() * 2.0 - 1.0)
    #
    #     # Now adjust length to match old:
    #     new_r = math.sqrt(dx ** 2.0 + dy ** 2)
    #     dx = dx * old_r / new_r
    #     dy = dy * old_r / new_r
    #
    #     print("old r2 was", old_r, end=' ')
    #     print("new r2 is", math.sqrt(dx ** 2.0 + dy ** 2.0))
    #
    #     # Now fudge one of the vx vy to change relative
    #     # vx,vy into dx,dy
    #
    #     p1.vx -= old_vx - dx
    #     p1.vy -= old_vy - dy
    #
    #     # Check if this worked, is the new relative V the same magnitude?
    #     rV = p1.subtract(p1.V(), p2.V())
    #     old_r = math.sqrt(rV[0] ** 2.0 + rV[1] ** 2.0)  # x^2 + y^2
    #     print("after rotation, magnitude of vx,vy is", old_r)
    #     print("ending   relative velocity is:", rV)
    #     print()

    # if False:
    #     print("Randomly Rotate relative vy,vz vector!")
    #     rV = p1.subtract(p1.V(), p2.V())
    #     old_vy = rV[1]
    #     old_vz = rV[2]
    #     old_r = math.sqrt(rV[1] ** 2.0 + rV[2] ** 2.0)  # y^2 + z^2
    #     print("starting relative velocity is:", rV)
    #
    #     # New random y,z direction:
    #     dy = maxV * (random.random() * 2.0 - 1.0)
    #     dz = maxV * (random.random() * 2.0 - 1.0)
    #
    #     # Now adjust length to match old:
    #     new_r = math.sqrt(dy ** 2.0 + dz ** 2)
    #     old_r2 = old_r * 2.0  # double magnitude to see what happens
    #     old_r2 = old_r
    #     dy = dy * old_r2 / new_r
    #     dz = dz * old_r2 / new_r
    #
    #     print("old r2 was", old_r, end=' ')
    #     print("new r2 is", math.sqrt(dy ** 2.0 + dz ** 2.0))
    #
    #     # Now fudge one of the vy vz to change relative
    #     # vy,vz into dy,dz
    #
    #     p1.vy -= old_vy - dy
    #     p1.vz -= old_vz - dz
    #
    #     # Check if this worked, is the new relative V the same magnitude?
    #     rV = p1.subtract(p1.V(), p2.V())
    #     old_r = math.sqrt(rV[1] ** 2.0 + rV[2] ** 2.0)  # y^2 + z^2
    #     print("after rotation, magnate of vy,vz is", old_r)
    #     print("ending   relative velocity is:", rV)
    #     print()

    # magnetic_test2(world)

    # sys.exit(1)


def magnitude(vec: ndarray) -> float:
    """ Compute length of 3D vector. """
    return np.linalg.norm(vec)


# def magnetic_test2(world):
#     for p1 in world:
#         p1.zero_force()
#         for p2 in world:
#             p1.add_force(p2)  # ES Only
#
#     for i in range(len(world)):
#         p1 = world[i]
#         print("%s%d X,Y %.1f, %.1f" % (
#             p1.symbol, i, p1.cur_state_r0 / Angstrom, p1.cur_state_r1 / Angstrom))
#
#     for i in range(len(world)):
#         p1 = world[i]
#         print("V of %s%i" % (p1.symbol, i), p1.v(),
#               "%5.3fc %s" % (magnitude(p1.v()) / CONST_C, p1.symbol))
#
#     # p1p22 = p1.magneticForce3(p2)
#     # p2p12 = p2.magneticForce3(p1)
#
#     for i in range(len(world)):
#         p1 = world[i]
#         print("ES  Total Force on %s%d" % (p1.symbol, i), p1.dt())
#
#     total_mag_force = (0.0, 0.0, 0.0)
#
#     for i in range(len(world)):
#         p1 = world[i]
#         for j in range(len(world)):
#             if i == j:
#                 continue
#             p2 = world[j]
#             f = p1.magnetic_force_total(p2)
#             print(
#                 "%s%d.magneticForceTotal(%s%d)" % (p1.symbol, i, p2.symbol, j),
#                 end=' ')
#             print(f, end=' ')
#             print(magnitude(f))
#
#             total_mag_force = vector_sum(total_mag_force, f)
#
#     print()
#     print("Total Mag force:", total_mag_force)
#     print()


def crt_clear():
    sys.stdout.write("\033[2J")


def crt_clear_and_home(row=0, col=0):
    crt_clear()
    crt_goto(row, col)


def crt_goto(row, col):
    sys.stdout.write("\033[%d;%dH" % (row, col))


def crt_mode(m=None):
    if m is None:
        m = 0  # normal mode
    elif m == "normal":
        m = 0
    elif m == "bold":
        m = 1
    elif m == "underline":
        m = 2
    elif m == "blinking":
        m = 5
    elif m == "reverse video":
        m = 7

    if m not in [0, 1, 2, 5, 7]:
        raise Exception("Invalid mode")

    sys.stdout.write(f"\033[{m}m")


def cross(v1, v2):
    # Cross product of two 3d vectors
    # returns a 3D vector
    x = v1[1] * v2[2] - v1[2] * v2[1]
    y = v1[2] * v2[0] - v1[0] * v2[2]
    z = v1[0] * v2[1] - v1[1] * v2[0]
    return x, y, z


def dot(v1, v2):
    """
    Dot product of two 3d vectors
    returns a magnitude value
    """
    sum_product = 0.0
    if len(v1) != len(v2):
        print("Vectors of different length in dot()", v1, v2)
        sys.exit(1)
    for i in range(len(v1)):
        sum_product += v1[i] * v2[i]
    return sum_product


def product(s, v):
    """ Multiply a scalar times a 3D vector and return a vector """
    return s * v[0], s * v[1], s * v[2]


def add(v1, v2):
    return v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2]


def subtract(v1, v2):
    return v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]


class ParticleState:
    def __init__(self, p: 'Particle'):
        """ Simulation State of Particle p. """
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
            j = 1.0 + Electron().mass / Proton().mass
        else:
            j = 1.0 + Proton().mass / Electron().mass

        # r = (p1.ke * p1.charge**2) / (p1.mass * v**2 * j**2)
        # d = r*j
        f = self.p.mass * magnitude(self.v) ** 3 * j ** 2 / \
            (2.0 * math.pi * CONST_KE * self.p.charge ** 2)
        return f

    def zero_force(self):
        """ Zero force vector. """
        self.f[:] = 0.0

    def add_force(self, p_state: 'ParticleState'):
        """ Add force on self created by p_state. """

        if self.p is p_state.p:
            # Particles cause no force on self
            return

        # Add Electrostatic force
        # TODO might need to be R limited?
        self.f += self.es_force(p_state)

        # Add electro-drag force
        if doMagnetic:
            self.f += self.v_force(p_state)

    def add_static_force(self):
        """ Add particle static forces to state. """
        self.f += self.p.static_f

    def es_force(self, p_state: 'ParticleState'):
        """ Electrostatic force on self by p_state per Coulomb's law. """
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

    # def v_force_test(self, p_state: 'ParticleState'):
    #     n = timeit.timeit('nv = self.v_force_new(p_state)', '', number=100, globals=locals())
    #     o = timeit.timeit('ov = self.v_force_old(p_state)', '', number=100, globals=locals())
    #     nv = self.v_force_new(p_state)
    #     if n < o:
    #         # print(end='.')
    #         print(f"New Faster - New: {n:.3f} Old: {o:.3f}")
    #     else:
    #         print()
    #         print(f"Old Faster - New: {n:.3f} Old: {o:.3f}")
    #         ov = self.v_force_old(p_state)
    #         print("  nv:", nv)
    #         print("  ov:", ov)
    #
    #     return nv

    def v_force(self, p_state: 'ParticleState'):
        """
            Faster version of v_force_old()
            Runs about 2x faster most the time.
        """
        # return self.v_force_old(p_state)
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
        # TODO verify this is working as expected. Testing needed.
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
        """ Distance**2 between self and p_state.
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
    """ The root class for protons and electrons. """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        # Units all standard SI: meters, seconds, Newtons, Kg, Coulombs
        self.cur_state = ParticleState(self)
        self.end_state = ParticleState(self)
        self.cur_state.r[:] = (x, y, z)

        self.static_f = np.zeros(3)  # Static forces for hack

        self.lock_in_place = False  # Don't move if true.

        self.avgKE = 0.0  # Running average of KE
        self.avgPE = 0.0  # Running average of PE

        self.charge = 0.0  # Defined by subclass for e and p
        self.mass = 0.0  # Defined in subclass for e and p
        self.symbol = 'e'  # Defined in subclass for e and p

    def r(self):
        """ Current 3D Position Vector. """
        return self.cur_state.r

    def v(self):
        """ Current 3D Velocity vector. """
        # return self.vx, self.vy, self.vz
        return self.cur_state.v

    def f(self):
        """ Current 3D force vector. """
        # return self.fx, self.fy, self.fz
        return self.cur_state.f

    def end_r(self):
        """ End state 3D Position Vector. """
        return self.end_state.r

    def end_v(self):
        """ End state Velocity vector. """
        # return self.vx, self.vy, self.vz
        return self.end_state.v

    def end_f(self):
        """ End state Force vector. """
        return self.end_state.f

    def zero_end_force(self):
        # self.end_fx = 0.0
        # self.end_fy = 0.0
        # self.end_fz = 0.0
        self.end_state.f[:] = 0.0

    # def add_force(self, p: 'Particle', p_state: ParticleState):  # and set end force as well
    #     if p is self:
    #         return
    #
    #     # dx = (self.cur_state_r0 - p.cur_state_r0)
    #     # dy = (self.cur_state_r1 - p.cur_state_r1)
    #     # dz = (self.cur_state_r2 - p.cur_state_r2)
    #     dr = p_state.r - p_state.r
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
    #     if doMagnetic:
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

        r2, l2 = self.distance2(p)

        if r2 == 0.0:
            return np.zeros(3)  # Bogus but prevents DBZ -- should be infinity

        force = CONST_KE * self.charge * p.charge / r2

        r = math.sqrt(r2)

        # return force * dx / r, force * dy / r, force * dz / r
        return dr * (force / r)

    def gravity_force(self, p: 'Particle'):
        # Magnitude of gravity between self and other particle
        g = 6.67408e-11  # 2014 CODATA recommended value
        r2, l2 = self.distance2(p)
        f = g * self.mass * p.mass / r2
        return f

    def magnetic_field(self, p: 'Particle'):
        """
            Returns a 3D field vector B = (x,y,z)
            Calculate the magnetic field created at self, by p.
            B = (1e-7 q1 v1 x r_hat) / r^2
            r_hat is the unit vector pointing from p to self
        """

        r2, l2 = self.distance2(p)

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
        f_vec: ndarray = self.charge * np.cross(self.v(), b_vec)

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
        """ 1/2 m v**2 """
        # ke = 0.5 * self.mass * (vx ** 2.0 + vy ** 2.0 + vz ** 2.0)
        # # print "KE CALC vx,vy,vz:", vx, vy, vz
        # # print "KE CALC answer =", ke
        ke = 0.5 * self.mass * v.dot(v)
        return ke

    def set_kinetic_energy(self, ke):
        # Back calculate velocity using given ke -- keep direction the same
        new_v2 = ke / (0.5 * self.mass)
        # old_v2 = (self.vx ** 2.0 + self.vy ** 2.0 + self.vz ** 2.0)
        old_v2 = np.sum(self.v() ** 2)
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

    def distance2(self, p: 'Particle'):
        """ Distance**2 between self and p.
            :returns: distance**2, limited_distance**2
        """

        if p is self:
            return self.limited_distance2(0.0)

        # dx = (self.cur_state_r0 - p.cur_state_r0)
        # dy = (self.cur_state_r1 - p.cur_state_r1)
        # dz = (self.cur_state_r2 - p.cur_state_r2)
        #
        # d2 = dx ** 2.0 + dy ** 2.0 + dz ** 2.0

        d2 = np.sum((self.r() - p.r()) ** 2.0)

        return self.limited_distance2(d2)

    def end_distance2(self, p: 'Particle'):  # distance squared
        """ Distance**2 between self and p end states.
            :returns: distance**2, limited_distance**2
        """
        if p is self:
            return self.limited_distance2(0.0)

        # dx = (self.end_x - p.end_x)
        # dy = (self.end_y - p.end_y)
        # dz = (self.end_z - p.end_z)

        # d2 = dx ** 2.0 + dy ** 2.0 + dz ** 2.0

        d2: float = np.sum((self.end_r() - p.end_r()) ** 2.0)

        return self.limited_distance2(d2)

    @staticmethod
    def limited_distance2(d2: float):
        """ Limit distance to RLimit to solve computational problems.
            return (real, limited) tuple
        """
        return d2, max(d2, RLimit ** 2)

    def distance(self, p):
        r, r_limited = self.distance2(p)
        return math.sqrt(r), math.sqrt(r_limited)

    def end_distance(self, p):
        r, r_limited = self.end_distance2(p)
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

        global InsideRLimitCount
        InsideRLimitCount += 1

        x = CONST_KE * self.charge * p.charge / RLimit

        return x + (x / RLimit) * (RLimit - d)


class Electron(Particle):
    """ A Single Electron. """

    def __init__(self, x=0.0, y=0.0, z=0.0):
        Particle.__init__(self, x, y, z)
        self.charge = -1.60218e-19  # in Coulombs
        self.mass = 9.10938356e-31
        self.symbol = 'e'


class Proton(Particle):
    def __init__(self, x=0.0, y=0.0, z=0.0, n=1.0):
        Particle.__init__(self, x, y, z)
        self.charge = n * +1.60218e-19  # in Coulombs
        self.mass = n * 1.672621898e-27
        self.symbol = 'p'


class ParticleImage:
    """" Particle that can draw itself on screen. """

    # TODO -- bad way to factor this, do it better
    def __init__(self, p: Particle):
        self.p: Particle = p

    def draw_particle(self, screen, world):
        """ Draw Particle on Screen.
            Implement bounce logic.
            returns: True if energy reset is needed
        """
        # TODO bounce logic shouldn't be in draw routine.
        x = self.space_to_pixels(self.p.cur_state.r[0])
        y = self.space_to_pixels(self.p.cur_state.r[1])
        z = self.space_to_pixels(self.p.cur_state.r[2])

        if isinstance(self.p, Electron):
            color = BLACK
            size = 2
        else:
            color = RED
            size = 4

        min_scale = 1.0
        max_scale = 3.0
        size *= (z / float(screen_depth)) * (max_scale - min_scale) + min_scale
        size = abs(int(size))

        inset = 10  # Bounce at 10 pixels from edge of screen
        if isinstance(self.p, Proton):
            inset = 40

        # e_change = 0.25   # Remove energy on bounce
        e_change = 1.00  # Don't remove energy on bounce

        reset_energy_needed = False
        bounce = False

        if x < inset and self.p.v()[0] < 0:
            self.p.v()[0] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r0 = self.pixels_to_space(inset)
            bounce = True

        if x > screen_width - inset and self.p.v()[0] > 0:
            self.p.v()[0] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r0 = self.pixels_to_space(screen_width - inset)
            bounce = True

        if y < inset and self.p.v()[1] < 0:
            self.p.v()[1] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r1 = self.pixels_to_space(inset)
            bounce = True

        if y > screen_height - inset and self.p.v()[1] > 0:
            self.p.v()[1] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r1 = self.pixels_to_space(screen_height - inset)
            bounce = True

        if z < inset and self.p.v()[2] < 0:
            self.p.v()[2] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r2 = self.pixels_to_space(inset)
            bounce = True

        if z > screen_depth - inset and self.p.v()[2] > 0:
            self.p.v()[2] *= -1
            self.p.set_kinetic_energy(self.p.kinetic_energy() * e_change)
            reset_energy_needed = True
            self.p.cur_state_r2 = self.pixels_to_space(screen_depth - inset)
            bounce = True

        global eBounceCount
        global pBounceCount

        if bounce:
            zero_momentum(world)  # fix momentum
            if isinstance(self.p, Proton):
                pBounceCount += 1
            else:
                eBounceCount += 1

        # Notice -- the way we change position below
        # to move particle back inside window frame without
        # adjusting velocity changes the total energy in the system.
        # So setting e_change = 1.0 does not conserve energy correctly.
        # We adjust the systems total target for fixing energy
        # on a bounce which forces us to accept this broken
        # energy total as correct.
        # But if we turn off the Bounce check, that allows the energyFix
        # system to keep fixing this error for us.  Hence, the test
        # below to turn Bounce flag off if e_change is 1.00

        if e_change == 1.0:
            reset_energy_needed = False  # Don't reset energy total -- fix it instead

        x = self.space_to_pixels(self.p.cur_state.r[0])
        y = self.space_to_pixels(self.p.cur_state.r[1])

        # print "x y is", x, y
        pygame.draw.circle(screen, color, (x, y), size, 0)

        return reset_energy_needed

    @staticmethod
    def space_to_pixels(space):
        # 0,0 is the same in both and is the top left corner of the screen
        return int(pixelsPerAngstrom * space / Angstrom)

    @staticmethod
    def pixels_to_space(pixels):
        return pixels * Angstrom / pixelsPerAngstrom


def total_momentum(world: list[Particle]):
    """ Total Momentum vector for all particles in world[]. """
    s = np.zeros(3)
    for p in world:
        s += p.cur_state.momentum()
    return s


def total_kinetic_energy(world: list[Particle]):
    total_ke = 0.0

    for p in world:
        total_ke += p.kinetic_energy()

    return total_ke


def total_potential_energy(world: list[Particle]):
    total_pe = 0.0

    for i in range(len(world)):
        for j in range(i + 1, len(world)):
            total_pe += world[i].potential_energy(world[j])

    return total_pe


def total_end_kinetic_energy(world: list[Particle]):
    total_ke = 0.0

    for p in world:
        total_ke += p.kinetic_end_energy()

    return total_ke


def total_end_potential_energy(world: list[Particle]):
    total_pe = 0.0

    for i in range(len(world)):
        for j in range(i + 1, len(world)):
            total_pe += world[i].potential_end_energy(world[j])

    return total_pe


#####################################################################
# Normalize momentum
#####################################################################

def zero_momentum(world):
    """
        Adjust all velocities by a constant to make total momentum zero.
    """
    tm = total_momentum(world)
    total_mass = sum([p.mass for p in world])
    dv = -tm / total_mass
    for p in world:
        p.cur_state.v += dv

    # print("starting m is", tm)
    # print("ending m is", total_momentum(world))
    # tm = total_momentum(world)
    # total_mass = sum([p.mass for p in world])
    # dv = -tm / total_mass
    # for p in world:
    #     p.cur_state.v += dv
    # print("starting m is", tm)
    # print("ending m is", total_momentum(world))
    # exit()


#####################################################################
# Move objects to place center of mass in the middle of the screen
#####################################################################

def center_mass(world: list[Particle]):
    """ Move all objects to put the center of mass in the center of
        the Screen.
    """
    center_m = center_of_mass(world)

    center_x = (screen_width / pixelsPerAngstrom * Angstrom) / 2.0
    center_y = (screen_height / pixelsPerAngstrom * Angstrom) / 2.0
    center_z = (screen_depth / pixelsPerAngstrom * Angstrom) / 2.0

    for p in world:
        # TODO improve this
        p.cur_state.r[0] += center_x - center_m[0]
        p.cur_state.r[1] += center_y - center_m[1]
        p.cur_state.r[2] += center_z - center_m[2]


def center_of_mass(world: List[Particle]) -> ndarray:
    """ return center of mass of the world as pos vector. """
    # cx = cy = cz = 0.0
    # tm = 0.0
    # for p in world:
    #     cx += p.cur_state_r0 * p.mass
    #     cy += p.cur_state_r1 * p.mass
    #     cz += p.cur_state_r2 * p.mass
    #     tm += p.mass
    # x = cx / tm
    # y = cy / tm
    # z = cz / tm
    # return x, y, z
    c = np.zeros(3)
    total_mass = 0.0
    for p in world:
        c += p.cur_state.r * p.mass
        total_mass += p.mass
    return c / total_mass


# fastTest()
# neutron_gravity_test()
# magneticTest()
# forceCircleTest()
# mag_circle_test()

""" Speed of light experiment.
    2021-02-11.
    Create a 1D string of electrons and make them stable
    by adding and extra artificial force to hold it in place.
    Then wiggle the first one, and watch how the rest respond.
    I'm hoping it will it show a wave moving through the string
    even though the simulation uses no time delay for the coulomb
    force moving them.
    
    First calculation.
    The 10 electrons, passed the energy of the first e to the 9th
    in this much time.  This is where the KE of the 9th maxed.
    spacing is 2/10 A so 9 times that is
    18/10 A in this much time:
    Time now is 0.0000000000000000243099999999994800664042
    
    Output is:
    Speed: speed=7404360.349, c=299792458.000 speed/c=0.025

"""


def sol_test():
    # speed = ((18/10) * Angstrom) / 0.00000000000000002430999999
    # print(f"Speed: {speed=:.3f}, {c=:.3f} {speed/c=:.3f}")
    # exit()
    # num = 8
    # spacing = (2 / num) * Angstrom
    # for cnt in range(num):
    #     pi_world.append(ParticleImage(Electron(cnt * spacing, 0.0, 0.0)))
    # # pi_world[0].p.lock_in_place = True
    # pi_world[0].p.vx = .5*c
    # # TODO - just a marker to find this code
    world: List[Particle] = [Proton(0.2 * Angstrom, 0.0, 0.0),
                             Proton(0.5 * Angstrom, 0.0, 0.0),
                             Proton(1.0 * Angstrom, 0.1 * Angstrom, 0.0),
                             Proton(1.2 * Angstrom, 0.1 * Angstrom,
                                    1.0 * Angstrom),
                             ]
    e1 = Electron(0.0 * Angstrom, 0.0, 0.0)
    world.append(e1)
    e1.cur_state.v[1] = .011 * CONST_C
    e2 = Electron(0.8 * Angstrom, 0.0, 0.0)
    world.append(e2)
    e2.cur_state.v[1] = .011 * CONST_C
    world.append(Electron(0.0 * Angstrom, 0.4 * Angstrom, 0.2 * Angstrom))
    world.append(Electron(0.0 * Angstrom, 0.1 * Angstrom, 0.2 * Angstrom))

    sim = Simulation(world)
    sim.run()


# p1 = Proton(Angstrom * 0.0, 0.0, 0.0 * Angstrom, n=3.0)
# e1 = Electron(Angstrom * 0.25, 0.0, 0.0 * Angstrom)
# e2 = Electron(Angstrom * -0.25, 0.0 * Angstrom, 0.0 * Angstrom)

# e1.vy = 200000.0
# e1.vy = 100000.0
# e1.vy = 0.0
# e1.vy = 1600000.0
# e1.vy = 16000000.0
# # e1.vx = 800000.0
# e1.vy = 3000000.0

# # This is what the below calculates for the perfect circular orbit
# # Velocity for the ep pair spaced at 0.25A
# e1.vy = 3181993.44316 + 100
# # p1.vy = -3181993.44316 * e1.mass / p1.mass
# e2.vy = -3181993.44316
# # e2.vz = -3181993.44316 / 10
# # e1.vy = p1.vy = 0.0
# # e1.vz = 800000.0

# if False:  # I have no clue what this code does 5-11-2018 CW
#     # 2021-02-06 it's calculating the velocity of a circular orbit.
#     r = Angstrom * 0.25
#     print("r is", r)
#     r = e1.x - p1.x
#     print("r is", r)
#     rCenter = r * p1.mass / (p1.mass + e1.mass)
#     e1.vz = e1.c * math.sqrt(
#         abs(1e-7 * e1.charge * p1.charge * rCenter / (r * r * e1.mass)))
#     rCenter = r * e1.mass / (p1.mass + e1.mass)
#     p1.vz = -p1.c * math.sqrt(
#         abs(1e-7 * e1.charge * p1.charge * rCenter / (r * r * p1.mass)))
#     print("p1.vz is", p1.vz, "momentum is", p1.momentum())
#     print("e1.vz is", e1.vz, "momentum is", e1.momentum())
#     print("p1.vz - e1.vz", p1.vz - e1.vz)
#
#     p1.vy = p1.vz
#     e1.vy = e1.vz
#     p1.vz = 0.0
#     e1.vz = 0.0
#     e1.vy /= 2.0

# sys.exit(1)
# e1.vz = 1000000.0

# Does adding a third electron away from an orbiting pair get sucked in?
# p2 = Proton(100.0 * Angstrom, 0.0 * Angstrom, 0.0 * Angstrom)
# e2 = Electron(100.25 * Angstrom, 0.0 * Angstrom, 0.0 * Angstrom)
# e2.vy = e1.vy # copy orbital velocity for second pair
# p2.vy = p1.vy


# p2 =   Proton(Angstrom * 0.5, Angstrom * 0.0,  0.0 * Angstrom)
# e2 = Electron(Angstrom * 1.5, Angstrom * 0.0, -0.5 * Angstrom)
# e2.vy = 1600000.0
# e2.vy = 800000.0

# p3 = Proton(Angstrom * 6.0, Angstrom * 0.0, 0.0)
# e3 = Electron(Angstrom * 7.0, Angstrom * 0.0, 0.0)
# e3.vy = 800000.0
# e3.vy = 1600000.0
#
# pi_world = []

# if True:
#     pi1 = ParticleImage(p1)
#     pi_world.append(pi1)
#     pi2 = ParticleImage(e1)
#     pi_world.append(pi2)
#     # pi3 = ParticleImage(p2)
#     # pi_world.append(pi3)
#     pi4 = ParticleImage(e2)
#     pi_world.append(pi4)

# if False:
#     pi3 = ParticleImage(p2)
#     pi_world.append(pi3)
#     if False:
#         pi4 = ParticleImage(e2)
#         pi_world.append(pi4)

# if False:
#     pi5 = ParticleImage(p3)
#     pi_world.append(pi5)
#     pi6 = ParticleImage(e3)
#     pi_world.append(pi6)

# if False:  # Two electrons
#     e1 = Electron(Angstrom * 0.1, Angstrom * 1.0, 0.0)
#     pi_world.append(ParticleImage(e1))
#     e2 = Electron(Angstrom * 5.1, Angstrom * 1.0, 0.0)
#     pi_world.append(ParticleImage(e2))

# if False:
#     # Random particles
#
#     p = 2
#     e = 2
#     percentOfScreen = 0.10
#     e1 = Electron()  # Just need any particle
#     momentum = e1.mass * e1.c / 1000.0
#     world = []
#     x_width = (screen_width * percentOfScreen / pixelsPerAngstrom) * Angstrom
#     y_width = (screen_height * percentOfScreen / pixelsPerAngstrom) * Angstrom
#
#     for i in range(p):
#         x = random.random() * x_width
#         y = random.random() * y_width
#         z = random.random() * min(x_width, y_width)
#         # z = 0.0
#         p = Proton(x, y, z)
#         world.append(p)
#         p.vx = random.random() * momentum / p.mass
#         p.vy = random.random() * momentum / p.mass
#         p.vz = random.random() * momentum / p.mass
#
#     for i in range(e):
#         x = random.random() * x_width
#         y = random.random() * y_width
#         z = random.random() * min(x_width, y_width)
#         # z = 0.0
#         p = Electron(x, y, z)
#         world.append(p)
#         p.vx = random.random() * momentum / p.mass
#         p.vy = random.random() * momentum / p.mass
#         p.vz = random.random() * momentum / p.mass
#
#     for p in world:
#         pi_world.append(ParticleImage(p))


class Simulation:
    """ A 3D simulation of electrons and protons interacting.
        Includes a graphic display window.
    """

    def __init__(self, world: List[Particle]):
        self.world = world

        self.fps_limit = 500    # pygame thing - Frames Per Second
        self.fps_avg = 1.0      # Computed speed of simulation
        self.run_me = True
        self.cycle_count = 0
        self.now = 0.0

        self.lastD = 0.0
        self.last_now = self.now
        self.lastV = 0.0
        self.lastA = 0.0

        self.p_pair_last_distance: Dict[
            Tuple[float, float], Tuple[float, float]] = {}

        self.pi_world: List[ParticleImage] = [ParticleImage(p) for p in
                                              self.world]

        center_mass(self.world)
        zero_momentum(self.world)
        self.init_forces()

        self.total_ke = total_kinetic_energy(self.world)
        self.total_pe = total_potential_energy(self.world)

        self.last_ke = self.total_ke
        self.last_pe = self.total_pe

        self.starting_total_ke = self.total_ke
        self.starting_total_pe = self.total_pe

        self.energy_diff_last = 0.0
        self.energy_diff_max = 0.0
        self.max_vc = 0.0  # Max V in units of CONST_C

        self.dt = dtMin

        self.screen = pygame.display.set_mode(screen_size)
        self.clock = pygame.time.Clock()
        pygame.display.set_caption('Physics-Sim')

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

        # if False:
        # shouldn't be done here -- should be done in test code.
        #     # compute static starting force to keep particles in place
        #     p1.static_fx = - p1.fx
        #     p1.static_fy = - p1.fy
        #     p1.static_fz = - p1.fz
        #     p1.add_static_force()

    def compute_end_forces(self):
        """ Update end_state forces using end_state position and velocity. """
        for p1 in self.world:
            p1.end_state.zero_force()
            for p2 in self.world:
                p1.end_state.add_force(p2.end_state)
            p1.end_state.add_static_force()

    def move_all(self):
        """ Copy end_state to cur_state for all particles. """
        for p in self.world:
            p.move()

    def run(self):
        """ Run the simulation. """

        self.run_me = True
        loop_cnt = 0
        reset_energy_needed = False
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
            self.fps_avg += (fps - self.fps_avg) * .01
            last_time = now

            if loop_cnt % 5 == 0:
                reset_energy_needed = self.draw_world()

            self.total_ke = total_kinetic_energy(self.world)
            self.total_pe = total_potential_energy(self.world)

            if reset_energy_needed:
                # Bounce happened on display above.
                # If it's configured to reduce energy on a bounce
                # then we re-define the starting_totals so if we
                # are using the fix energy system, it will fix it
                # based on the energy as it is now, not as it was
                # when the simulation started.
                # Currently - we aren't removing energy on a bounce
                # and we aren't using the fix code, so this is not used.
                self.starting_total_ke = self.total_ke
                self.starting_total_pe = self.total_pe

            if loop_cnt % 20 == 0:
                self.print_stats()

                if 0.0 < Stop_at < self.now:
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

            ######################################################################
            # Fix total energy error
            # Adjust all velocities to remove energy error of energy_diff_start
            ######################################################################

            if energyFix:
                if energy_diff_start > total_ke2:
                    print("energy error greater than total KE error:",
                          energy_diff_start,
                          "total", total_ke2)
                    energy_diff_start = total_ke2 * 0.5
                    # Only try to take half of total KE out of the system
                    # time.sleep(1)
                    # sys.exit(1)

                for i in range(len(self.world)):
                    # print
                    # print "Adjust KE for particle", i, "energy_diff_start error is", energy_diff_start
                    p1 = self.world[i]
                    ke = p1.kinetic_energy()
                    # print "Current ke is", ke
                    # print "percent of of total is", ke/total_ke2*100.0, "%"
                    # print "amount added of total energy_diff_start is", ke/total_ke2*energy_diff_start
                    new_ke = ke - ke / total_ke2 * energy_diff_start
                    # print "new should be", new_ke
                    p1.set_kinetic_energy(new_ke)
                    # print "new KE is now", p1.kinetic_energy()

        pygame.quit()

    def draw_world(self):
        self.screen.fill(WHITE)
        reset_energy_needed = False
        s_list = sorted(self.pi_world, key=lambda arg: arg.p.cur_state.r[2])
        for p in s_list:
            reset_energy_needed |= p.draw_particle(self.screen, self.world)
        pygame.display.flip()
        return reset_energy_needed

    def print_stats(self):
        """ Print statistics for this loop. """

        crt_clear_and_home()

        print()
        print("PixelsPerAngstrom", pixelsPerAngstrom, end=' ')
        print("Screen (%.1fx%.1f) A" % (
            screen_width / pixelsPerAngstrom,
            screen_height / pixelsPerAngstrom))
        print(f"Time now is {self.now * 1000_000_000:.20f} ns", end='')
        print("  DT is: %4.1e" % self.dt, end='')
        print(f"  FPS:{self.fps_avg:.1f}")

        if Stop_at > 0.0:
            print(f"    Stop at {Stop_at * 1000_000_000:.20f} ns")

        p_vec = total_momentum(self.world)
        print(f"Total Momentum {p_vec[0]:8.1e} {p_vec[1]:8.1e} "
              f"{p_vec[2]:8.1e}", end='')
        print(f"   mag: {magnitude(p_vec):.1e}")
        print()
        print("doMagnetic:", doMagnetic, "  doMagneticInverse:",
              doMagneticInverse)
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
        print()

        print("Inside R limit ", InsideRLimitCount, " pBounce:", pBounceCount,
              " eBounce:", eBounceCount)

        self.last_ke = self.total_ke
        self.last_pe = self.total_pe

        self.print_proton_distance()

        print()

        self.print_particle_stats()

        print("Max Velocity: %5.3fc" % self.max_vc)

        print("Energy error for move is %12.4e" % self.energy_diff_last,
              end='')
        print(" max error %12.4e" % self.energy_diff_max, end='')
        print()

        ###########################################
        # Debug center of mass of pairs for
        # special e/p pair drift test
        ###########################################
        # lastCm = None
        # for i in range(len(self.world) - 1):
        #     p1 = self.world[i]
        #     p2 = self.world[i + 1]
        #     if isinstance(p2, Electron) and isinstance(p1, Proton):
        #         cm = center_of_mass((p1, p2))
        #         # print("Center of Mass of pair %20.15f %20.15f %20.15f" % (
        #         #     cm[0] / Angstrom, cm[1] / Angstrom, cm[2] / Angstrom))
        #         if lastCm is not None:
        #             d = math.sqrt((cm[0] - lastCm[0]) ** 2 +
        #                           (cm[1] - lastCm[1]) ** 2 +
        #                           (cm[2] - lastCm[2]) ** 2)
        #             # print("Distance is %20.15f" % (d / Angstrom))
        #             dd = d - last_d
        #             # print("Change   is %20.15f" % (dd / Angstrom))
        #             move_dt = self.now - self.last_now
        #             # print("DT is %20.30f" % move_dt)
        #             v = 0.0
        #             if move_dt != 0.0:
        #                 v = dd / move_dt  # Velocity
        #                 print(f"M/S is {dd / move_dt:.3e} <0 is moving together")
        #                 a = (v - self.lastV) / move_dt
        #                 print(f"Acceleration is {a:+.4e}")
        #                 self.lastA += (a - self.lastA) * .001
        #                 print(f"       avg A is {self.lastA:.4e}", end='')
        #             self.lastV = v
        #             self.last_now = self.now
        #             last_d = d
        #             gravity_force = 0.0
        #             gravity_force += p1.gravity_force(self.world[i - 1])
        #             gravity_force += p1.gravity_force(self.world[i - 2])
        #             gravity_force += p2.gravity_force(self.world[i - 2])
        #             gravity_force += p2.gravity_force(self.world[i - 1])
        #             # print("Gravity Force is", gravity_force)
        #             # es = p1.es_force(p2)
        #             # print("es from p1 to p2 is", es)
        #             # es = p1.es_force(self.world[i - 2])
        #             # print("es from p1 to e2 is", es)
        #             # print("es from 1 to p4 is", es)
        #             print(f"  Gravity A is {gravity_force / (p1.mass + p2.mass):.3e}")
        #         lastCm = cm
        # print()

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
                last_d, last_time = self.p_pair_last_distance.get(pair, (
                    d, self.now - self.dt))
                self.p_pair_last_distance[pair] = (
                    d, self.now)  # save for next time
                v = (d - last_d) / (
                        self.now - last_time)  # positive v if moving away
                print("Distance between", end=' ')
                print(f"P{i:02d} P{j:02d} {d / Angstrom:11.6f} A", end='')
                print(f"   v:{v:+.2e}")
                cnt += 1
                if cnt > 10:
                    break
            if cnt > 10:
                break

    def print_particle_stats(self):
        total_avg_ke = 0.0
        total_momentum_mag = 0.0

        for i in range(len(self.world)):
            p1 = self.world[i]
            ke = p1.kinetic_energy()
            p1.avgKE += (ke - p1.avgKE) * 0.0001
            total_avg_ke += p1.avgKE
            total_momentum_mag += magnitude(p1.cur_state.momentum())

        for i in range(len(self.world)):
            p1 = self.world[i]
            print(f"{p1.symbol}{i}", end='')
            print(f" vx:{p1.v()[0]:10.2e}  vy:{p1.v()[1]:10.2e}", end='')
            vc = magnitude(p1.v()) / CONST_C
            self.max_vc = max(self.max_vc, vc)
            print(f" {vc:0.5f}c", end='')
            print(" x:%6.2f A" % (p1.cur_state.r[0] / Angstrom), end=' ')
            # print " KE:%10.2e" % p1.avgKE
            # if total_avg_ke:
            #     print(" KE:%6.2f%%" % (p1.avgKE * 100 / total_avg_ke))
            # else:
            #     print(" KE:?????")
            if self.total_ke:
                print(" KE:%5.1f%%" % (
                        p1.kinetic_energy() * 100 / self.total_ke), end='')
            else:
                print(" KE:?", end='')
            momentum = magnitude(p1.cur_state.momentum())
            print(f"  p:{momentum:.2e}", end='')
            if total_momentum_mag:
                print(f" {momentum * 100 / total_momentum_mag:5.1f}%", end='')

            # Frequency of electron in orbit.
            if isinstance(p1, Electron) or isinstance(p1, Proton):
                f = p1.cur_state.f_from_v()
                print(f" f:{f / 1_000_000_000_000_000:.1e} PHz", end='')
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

            if energyFix2:
                # Ok, one last time, fudge velocity based on current position and
                # total force for this position This is just always good I think.
                # Probably shouldn't be an option.

                for p in self.world:
                    p.calculate_end_velocity(self.dt)

            # Now calculate total Energy after this move.

            total_ke2 = total_end_kinetic_energy(self.world)
            total_pe2 = total_end_potential_energy(self.world)

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

            if dtAdjust and total_ke2 != 0.0:
                # print "==DO DT ADJUST self.cycleCount", self.cycleCount, "abs(energy_diff)", abs(energy_diff),
                # print "total_ke2", total_ke2, "percent", abs(energy_diff) / total_ke2
                if self.dt < dtMax and self.cycle_count > 3 and abs(
                        self.energy_diff_last) / total_ke2 < 0.0001:
                    # print "SPEEDUP -- increase DT abs(diff)/total is", abs(energy_diff) / total_ke2
                    self.dt *= 2.0
                    self.dt = min(self.dt, dtMax)
                    # self.cycleCount = 0
                    # continue
                    # No need to restart -- take this move as fine but use a larger dt
                    # for the next move

                elif self.dt > dtMin and abs(
                        self.energy_diff_last) / total_ke2 > 0.001:
                    # print "SLOWDOWN -- reduce DT abs(diff)/total is", abs(energy_diff) / total_ke2
                    self.dt /= 2.0
                    self.dt = max(self.dt, dtMin)
                    for p1 in self.world:
                        p1.reset_state()
                    self.cycle_count = 0
                    continue

            self.energy_diff_max = max(self.energy_diff_max,
                                       abs(self.energy_diff_last))

            break  # End of DT adjust loop

        return energy_diff_start, total_ke2


if __name__ == '__main__':
    main()
