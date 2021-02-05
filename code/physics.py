#!/usr/bin/python

"""
    physics.py - Atomic particle simulation experments

    2016-05-18 Started this file.  Took it from xnet.py from my AI work.
    2021-02-05 Created PyCharm Project, Put on Github, converted
                to python3.9, and started major cleanup of this old code.

"""

import random
# import time
import sys
import math
# import heapq
# import numpy as np
# import matplotlib.pyplot as plt
from operator import attrgetter

# import os
import pygame

black = 0, 0, 0
white = 255, 255, 255
red = 255, 0, 0

Angstrom = 1.0e-10  # One Angstrom 10e-10 meters

# RLimit = 0.1 * Angstrom			# Radius limit hack
# RLimit = 0.0001 * Angstrom		# Radius limit hack
RLimit = 0.0000001 * Angstrom  # Radius limit hack

InsideRLimitCount = 0
eBounceCount = 0
pBounceCount = 0

FakeConstants = False

doMagnetic = False
doMagneticInverse = False

# dt = 2e-17
# dt = 2e-18
# dt = 2e-19

# dtMin = 1e-20
dtMin = 1e-30

# dtMax = 1e-10
dtMax = 1e-1

dt = dtMin
dtAdjust = True

energyFix = False  # fix based on total PE+KE at start -- doesn't work when magnatism added
energyFix2 = False

# screen_size = screen_width, screen_height = 1000, 800
screen_size = screen_width, screen_height = 600, 400
screen_depth = 1000  # z dimension

# pixelsPerAngstrom = 4000.0
# pixelsPerAngstrom = 0.00001
# pixelsPerAngstrom = 100000000.0
# pixelsPerAngstrom = 10000.0
# pixelsPerAngstrom = 5.0
pixelsPerAngstrom = 200.0


def mag_circle_test():
    # Copy forceCircleTest to crete magfCircleTest()
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
    z = 0.0
    for i in range(11):
        x = i * 0.1
        print("For x,y", x, y, end=' ')
        bTotal = 0.0
        steps = 4000
        dl = 2.0 * math.pi / steps
        lastpx = 0.0  # First point
        lastpy = 1.0
        for i in range(1, steps + 1):
            anglei = 2.0 * math.pi * i / steps
            px = math.sin(anglei)
            py = math.cos(anglei)
            dl = math.sqrt((lastpx - px) ** 2.0 + (lastpy - py) ** 2.0)
            dx = px - x
            dy = py - y
            r2 = dx ** 2 + dy ** 2
            r = math.sqrt(r2)
            # For 1 amp current
            b = 1e-7 * dl / r  # fake for expermenting
            b = 1e-7 * dl / 1.0  # Fake again
            b = 1e-7 * dl / r2  # Real
            bTotal += b
            lastpx = px
            lastpy = py
        print("total mag field is", bTotal)
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
    for i in range(9):
        x = i * 0.1
        print("For x,y", x, y, end=' ')
        fx = 0.0
        fy = 0.0
        fz = 0.0
        steps = 400
        for i in range(steps):
            for j in range(steps):
                anglei = 2.0 * math.pi * i / steps
                anglej = 2.0 * math.pi * j / steps
                px = math.sin(anglei)
                py = math.cos(anglei)
                pz = math.cos(anglej)
                dx = px - x
                dy = py - y
                dz = pz - z
                d2 = dx ** 2 + dy ** 2 + dz ** 2
                r = math.sqrt(d2)
                force = 1.0 / r
                force = 1.0 / d2
                # Had to add /px below to adjust for the ds factor
                # that our step pattern was wrong
                fx += force * dx / r * px
                fy += force * dy / r * px
                fz += force * dz / r * px
        print("total force is", fx, fy, fz)
    sys.exit(1)

    # Conclusion ...
    # Force at point inside cirucle only sums to a constant if the
    # force equation is 1/x  not 1/x^2.
    # Force sums to constant when we sum over entire spher for 1/x^22 force.
    # This makes sense because the size of the circle is 2pi r, or 4pi r^2 for surface
    # The use of 1/x^2 in physics is consistent with conservation in 3D space.


def fast_test():
    p1 = Proton(0.0, 00.0, 0.0)
    p1 = Electron(0.0, 00.0, 0.0)
    for i in range(60):
        x = (i - 30) * Angstrom / 10
        e1 = Electron(x, 0.0, 0.0)
        print("i", i, "x", x, x / Angstrom, "A", end=' ')
        print(p1.potentialEnergy(e1))

    sys.exit(1)


def neutronGravityTest():
    # Some experiment on the idea that a neutron is
    # is really a e p pair in orbit.  What would
    # be the attraction between two such systems?
    # Does it match gravity in relative force?
    # The answer was the attraction is not even
    # 1/x^2 in force!  Opps.  That theory got blown
    # out of the water!  The overlaping fields does
    # create a 1/x^2 field!
    # UGH -- 5-11-2018 thoughts: but do two systems in
    # orbit actually end up moving closer together or
    # further apart due to distroted orbits and orbital
    # intereactions -- this test did not answer that question
    # and if the effect appxoimates 1/x^2 then it could be
    # mistaken as 1/^x^2.  More testing required

    closeSpace = 0.0001 * Angstrom
    p1 = Proton(0.0, 0.0, 0.0)
    e1 = Electron(0.0 + closeSpace, 0.0, 0.0)

    p2 = Proton(100.0, 0.0, 0.0)
    e2 = Electron(100.0 + closeSpace, 0.0, 0.0)

    gravityForce = 0.0
    gravityForce += p2.gravityForce(p1)
    gravityForce += p2.gravityForce(e1)
    gravityForce += e2.gravityForce(p1)
    gravityForce += e2.gravityForce(e1)

    print("Gravity force between two is:", gravityForce)

    p2.zeroForce()
    p2.addForce(p1)
    print("em force p2 to p1", p2.fx)
    emForce = p2.fx

    p2.zeroForce()
    p2.addForce(e1)
    print("em force p2 to e1", p2.fx)
    emForce += p2.fx

    print("em force p2 to e1 and p1", emForce)

    sys.exit(1)


def magnetic_test():
    c = 299792458.0  # Speed of light m/s

    world = []

    p1 = Electron(0.0, 0.0, 0.0)
    p1 = Proton(0.0, 0.0, 0.0)
    # p1.vx = c/2.0
    # p1.vy = c/2.0
    # p1.vy = c/2.0
    # p1.vy = 2.0e8
    # p1.vz = 12.0e8
    world.append(p1)

    p2 = Electron(1.0 * Angstrom, 0.0, 0.0)
    p2 = Electron(1.0 * Angstrom, 1.0 * Angstrom, 0.0)
    p2.vy = c / 2.0
    # p2.vy = -c/4.0
    # p2.vy = 1.0e8
    # p2.vz = 6.0e8
    world.append(p2)

    p3 = Electron(1.0 * Angstrom, 0.0, 0.0)
    p3 = Electron(1.0 * Angstrom, 1.0 * Angstrom, 0.0)
    p3.vx = -c / 2.0
    p3.vy = -c / 2.0
    # world.append(p3)

    maxV = c / math.sqrt(3.0)
    # print "max v is", maxV, "max mag is", math.sqrt(3.0 * maxV**2.0), "c is", c
    # zz

    if False:
        speed = c / 2.0
        speed = maxV * random.random()
        speed = c
        p1.vx = (random.random() * 2.0 - 1.0)
        p1.vy = (random.random() * 2.0 - 1.0)
        p1.vz = (random.random() * 2.0 - 1.0)
        p1.vx, p1.vy, p1.vz = p1.product(speed / magnitude(p1.V()), p1.V())

        speed = c / 2.0
        speed = maxV * random.random()
        speed = c
        p2.vx = (random.random() * 2.0 - 1.0)
        p2.vy = (random.random() * 2.0 - 1.0)
        p2.vz = (random.random() * 2.0 - 1.0)
        p2.vx, p2.vy, p2.vz = p2.product(speed / magnitude(p2.V()), p2.V())

    if False:
        # Make the magnitude of the difference equal to c
        rV = p1.subtract(p1.V(), p2.V())
        s = magnitude(rV)
        p1.vx, p1.vy, p1.vz = p1.product(c / s, p1.V())
        p2.vx, p2.vy, p2.vz = p2.product(c / s, p2.V())

    if False:
        p1.vx = 1.0
        p1.vy = 2.0
        p1.vz = 3.4
        p2.vx = 100.0
        p2.vy = 200.0
        p2.vz = 300.4

    if False:
        # Time step test -- move particles forward in time
        # Based on velocity, see what dr and ds is!
        p0 = world[0]
        rp0 = world[0].R()
        rp1 = world[1].R()
        rBefore = p0.subtract(rp1, rp0)
        esBefore = world[0].esForce(world[1])

        dt = 1e-20

        magnetic_test2(world, dt=dt)

        for it in range(3):
            for i in range(len(world)):
                p1 = world[i]
                p1.calculateEndVelocity(dt)
                p1.calculateEndPosition(dt)

            for p1 in world:
                p1.zeroEndForce()
                for p2 in world:
                    p1.addEndForce(p2)

        for p1 in world:
            p1.move()

        print()
        print("Move Done!")
        print()

        rp0 = world[0].R()
        rp1 = world[1].R()
        rAfter = p0.subtract(rp1, rp0)
        rDiff = p0.subtract(rAfter, rBefore)
        esAfter = world[0].esForce(world[1])
        esDiff = p0.subtract(esAfter, esBefore)

        print(" dt is   ", dt)
        print(" r before", rBefore)
        print(" r after", rAfter)
        print(" dr =  after-before = ", rDiff)
        print(" dr diff magnatude", magnitude(rDiff))
        print(" es before  ", esBefore)
        print(" es after   ", esAfter)
        print(" es diff    ", esDiff)
        print(" es diff magnatude", magnitude(esDiff))
        print()

        magnetic_test2(world, dt=dt)
        sys.exit(1)

    magnetic_test2(world)
    sys.exit(1)

    # zz

    print("ZERO Monumtum!")
    print()
    zero_momentum(world)
    magnetic_test2(world)

    sys.exit(1)
    if False:
        p1.vx = - p1.vx

    # Change velocities in different ways and print results again...

    if False:
        # Change frame of reference
        print("Change inertial Frame of reference randomly")

        # Because my test particles are on the same x axis, this
        # vx component is totally ignored in the equations so
        # I can assign random values without effecting the total force!
        # FALSE -- not true.  Changing dx changes fy and fz but not
        # fx.  So it is important, unless the vy and vz are zero.
        dx = maxV * (random.random() * 2.0 - 1.0)
        p1.vx += dx
        # dx = maxV * (random.random() * 2.0 - 1.0)
        p2.vx += dx

        dy = maxV * (random.random() * 2.0 - 1.0)

        p1.vy += dy
        p2.vy += dy

        dz = maxV * (random.random() * 2.0 - 1.0)

        p1.vz += dz
        p2.vz += dz

    if False:
        print("Change X velocities randomly")
        # turns out, that doing this, has no
        # effect on the X force when the particles
        # are on the x axis!
        dx = maxV * (random.random() * 2.0 - 1.0)
        p1.vx += dx
        dx = maxV * (random.random() * 2.0 - 1.0)
        p2.vx += dx

    if False:
        print("Change Y velocities randomly")
        # turns out, that doing this, has no
        # effect on the Z force when the particles
        # are on the x axis!
        dy = maxV * (random.random() * 2.0 - 1.0)
        p1.vy += dy
        dy = maxV * (random.random() * 2.0 - 1.0)
        p2.vy += dy

    if False:
        print("Change Z velocities randomly")
        # turns out, that doing this, has no
        # effect on the Y force when the particles
        # are on the x axis!
        dz = maxV * (random.random() * 2.0 - 1.0)
        p1.vz += dz
        dz = maxV * (random.random() * 2.0 - 1.0)
        p2.vz += dz

    if False:
        print("Randomly Rotate relative vx,vy vector!")
        # Coded this by mistake.  Meant to code vy,vz rotation which
        # I did below.  When we rotate vx,vy  vextor as I'm doinging
        # here, all forces change and the total magnature of the force
        # changes as well!
        rV = p1.subtract(p1.V(), p2.V())
        oldvx = rV[0]
        oldvy = rV[1]
        oldr = math.sqrt(rV[0] ** 2.0 + rV[1] ** 2.0)  # x^2 + y^2
        print("starting relative velocity is:", rV)

        # New random x,y direction:
        dx = maxV * (random.random() * 2.0 - 1.0)
        dy = maxV * (random.random() * 2.0 - 1.0)

        # Now adjust length to match old:
        newr = math.sqrt(dx ** 2.0 + dy ** 2)
        dx = dx * oldr / newr
        dy = dy * oldr / newr

        print("old r2 was", oldr, end=' ')
        print("new r2 is", math.sqrt(dx ** 2.0 + dy ** 2.0))

        # Now fudge one of the vx vy to change relative
        # vx,vy into dx,dy

        p1.vx -= oldvx - dx
        p1.vy -= oldvy - dy

        # Checck if this worked, is the new relatie V the same magnatude?
        rV = p1.subtract(p1.V(), p2.V())
        oldr = math.sqrt(rV[0] ** 2.0 + rV[1] ** 2.0)  # x^2 + y^2
        print("after rotation, magnature of vx,vy is", oldr)
        print("ending   relative velocity is:", rV)
        print()

    if False:
        print("Randomly Rotate relative vy,vz vector!")
        rV = p1.subtract(p1.V(), p2.V())
        oldvy = rV[1]
        oldvz = rV[2]
        oldr = math.sqrt(rV[1] ** 2.0 + rV[2] ** 2.0)  # y^2 + z^2
        print("starting relative velocity is:", rV)

        # New random y,z direction:
        dy = maxV * (random.random() * 2.0 - 1.0)
        dz = maxV * (random.random() * 2.0 - 1.0)

        # Now adjust length to match old:
        newr = math.sqrt(dy ** 2.0 + dz ** 2)
        oldr2 = oldr * 2.0  # double magnatude to see what hapens
        oldr2 = oldr
        dy = dy * oldr2 / newr
        dz = dz * oldr2 / newr

        print("old r2 was", oldr, end=' ')
        print("new r2 is", math.sqrt(dy ** 2.0 + dz ** 2.0))

        # Now fudge one of the vy vz to change relative
        # vy,vz into dy,dz

        p1.vy -= oldvy - dy
        p1.vz -= oldvz - dz

        # Checck if this worked, is the new relatie V the same magnatude?
        rV = p1.subtract(p1.V(), p2.V())
        oldr = math.sqrt(rV[1] ** 2.0 + rV[2] ** 2.0)  # y^2 + z^2
        print("after rotation, magnature of vy,vz is", oldr)
        print("ending   relative velocity is:", rV)
        print()

    magnetic_test2(world)

    sys.exit(1)


def magnitude(v):
    return math.sqrt(v[0] ** 2.0 + v[1] ** 2.0 + v[2] ** 2.0)


def magnetic_test2(world, dt=None):
    doMagnetic = True
    doMagneticInverse = True

    for p1 in world:
        p1.zeroForce()
        for p2 in world:
            p1.addForce(p2)  # ES Only

    c = 299792458.0  # Speed of light m/s

    for i in range(len(world)):
        p1 = world[i]
        print("%s%d X,Y %.1f, %.1f" % (
        p1.symbol, i, p1.x / Angstrom, p1.y / Angstrom))

    for i in range(len(world)):
        p1 = world[i]
        print("V of %s%i" % (p1.symbol, i), p1.V(),
              "%5.3fc %s" % (magnitude(p1.V()) / c, p1.symbol))

    # p1p22 = p1.magneticForce3(p2)
    # p2p12 = p2.magneticForce3(p1)

    for i in range(len(world)):
        p1 = world[i]
        print("ES  Total Force on %s%d" % (p1.symbol, i), p1.F())

    totalMagForce = (0.0, 0.0, 0.0)

    for i in range(len(world)):
        p1 = world[i]
        for j in range(len(world)):
            if i == j:
                continue
            p2 = world[j]
            f = p1.magnetic_force_total(p2, dt=dt)
            print(
                "%s%d.magneticForceTotal(%s%d)" % (p1.symbol, i, p2.symbol, j),
                end=' ')
            print(f, end=' ')
            print(magnitude(f))

            totalMagForce = vectorSum(totalMagForce, f)

    print()
    print("Total Mag force:", totalMagForce)
    print()

    # zz


def crtClear():
    sys.stdout.write("\033[2J")


def crtClearAndHome(row=0, col=0):
    crtClear()
    crtGOTO(row, col)


def crtGOTO(row, col):
    sys.stdout.write("\033[%d;%dH" % (row, col))


def crtMode(m=None):
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
        raise Exception("Invlaid mode")

    sys.stdout.write("\033[%dm" % (m))


class Particle:
    # The basic particle code for protons and electons
    def __init__(self, x=0.0, y=0.0, z=0.0):
        # Units all standard SI meters, seconds, Newtons, Kg, Coulombs
        self.x = x
        self.y = y
        self.z = z
        self.vx = 0.0
        self.vy = 0.0
        self.vz = 0.0
        self.fx = 0.0
        self.fy = 0.0
        self.fz = 0.0
        self.c = 299792458.0  # Speed of light m/s -- defined exact constant
        self.ke = 8.9875517873681764e9  # Coulomb's constant (1/4 pi e) written as K(sub)e
        self.ke = self.c * self.c * 1.0e-7  # exactly the same as above

        if FakeConstants:
            self.ke = 1.0  # Coulomb's constant (1/4 pi e) written as K(sub)e
        # self.charge = -1.0

        self.avgKE = 0.0  # Running average of KE
        self.avgPE = 0.0  # Running average of PE

    def R(self):
        return (self.x, self.y, self.z)

    def V(self):
        return (self.vx, self.vy, self.vz)

    def F(self):
        return (self.fx, self.fy, self.fz)

    def zeroForce(self):  # and end force as well
        self.fx = 0.0
        self.fy = 0.0
        self.fz = 0.0
        self.zeroEndForce()

    def zeroEndForce(self):
        self.endFx = 0.0
        self.endFy = 0.0
        self.endFz = 0.0

    def addForce(self, p):  # and set end force as well
        if p is self:
            return

        dx = (self.x - p.x)
        dy = (self.y - p.y)
        dz = (self.z - p.z)

        r2, l2 = self.distance2(p)

        if r2 == 0.0:
            return  # Bogus but prevents DBZ

        force = self.ke * (self.charge * p.charge) / l2

        r = math.sqrt(r2)

        self.fx += force * dx / r
        self.fy += force * dy / r
        self.fz += force * dz / r

        if doMagnetic:
            f = self.magnetic_force_total(p)
            self.fx += f[0]
            self.fy += f[1]
            self.fz += f[2]

        self.endFx = self.fx
        self.endFy = self.fy
        self.endFz = self.fz

    def esForce(self, p):
        # Electrostic force between self and p per coulomb's law.
        # Force on self, caused by p.
        # real force, not Rlimit limited force
        # Returns 0,0,0 instead of infinity for two particls located
        # on same spot.

        if p is self:
            return (0.0, 0.0, 0.0)

        dx = (self.x - p.x)
        dy = (self.y - p.y)
        dz = (self.z - p.z)

        r2, l2 = self.distance2(p)

        if r2 == 0.0:
            return (
            0.0, 0.0, 0.0)  # Bogus but prevents DBZ -- should be infinity

        force = self.ke * self.charge * p.charge / r2

        r = math.sqrt(r2)

        return (force * dx / r, force * dy / r, force * dz / r)

    def gravityForce(self, p):
        G = 6.67408e-11  # 2014 CODATA recommended value
        r2, l2 = self.distance2(p)
        f = G * self.mass * p.mass / r2
        return f

    def cross(self, v1, v2):
        # Cross product of two 3d vectors
        # returns a 3D vector
        x = v1[1] * v2[2] - v1[2] * v2[1]
        y = v1[2] * v2[0] - v1[0] * v2[2]
        z = v1[0] * v2[1] - v1[1] * v2[0]
        return (x, y, z)

    def dot(self, v1, v2):
        # Dot product of two 3d vectors
        # returns a magnatude value
        sum = 0.0
        if len(v1) != len(v2):
            print("Vectors of different length in dot()", v1, v2)
            sys.exit(1)
        for i in range(len(v1)):
            sum += v1[i] * v2[i]
        return sum

    def product(self, s, v):
        # Multiply a scaler times a 3D vector and return a vector
        return (s * v[0], s * v[1], s * v[2])

    def add(self, v1, v2):
        return (v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2])

    def subtract(self, v1, v2):
        return (v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2])

    def magneticField(self, p):
        # Returns a 3D field vector B = (x,y,z)
        # Calculate the magnetic field created on self, by p.
        # B = (1e-7 q1 v1 x rHat) / r^2
        # rHat is the unit vector pointing from p to self

        (r2, l2) = self.distance2(p)

        if r2 == 0.0:
            return (0.0, 0, 0, 0.0)

        # print " distance is", r2, math.sqrt(r2)

        rHat = (self.x - p.x, self.y - p.y, self.z - p.z)
        rHat = self.product(1.0 / math.sqrt(r2), rHat)

        # The books say based on current flow, this should be correct:
        # This follows the right hand rule for postive current flow
        B = self.product(1e-7 * p.charge / r2, self.cross(p.V(), rHat))

        # This is inverted
        B = self.product(-1.0, B)  # backwards

        # The correct way, means that motion causes electrons to push HARDER
        # apart from each other, and ep pairs to pull harder together.  So
        # when the velociety == c, the magnetic force equals the electostatic force.
        # If they add, as the books imply they do, it means that going the speed
        # of light means the force doubles.  But that is fucking pointless.
        # It doesn't make the speed of light special in any sense related to
        # orbital mechanics as far as I can deduce.
        # To make it special, the speed of light has to be the point where these
        # forces cancel each other out!  Which makes me logically deduce that
        # This sign has to be backwards to work correctly.
        # So I'm going to play with making it backwards from how the books all
        # say it should be to see what happens.

        return B

    def magneticForce(self, p):
        # Returns a 3D force vector (x,y,z)
        # Calculate the force created on self, by the magnetic
        # field generated by p.
        # B = (1e-7 q1 v1 x rHat) / r^2
        # F = q2 V2 x B
        # F = q2 V2 X (1e-7 q1 v1 x rHat) / r^2
        # rHat is the unit vector pointing from p to self

        B = self.magneticField(p)

        # print "mag force"
        # print " B is", B

        F = self.product(self.charge, self.cross(self.V(), B))

        return F

    def magneticForce2(self, p):
        # Same as mag force but assume the magenetic field
        # is moving with p.  So use our relative speed
        # as our speed moving THROUGH the mag field.
        # Returns a 3D force vector (x,y,z)
        # Calculate the force created on self, by the magnetic
        # field generated by p.
        # B = (1e-7 q1 v1 x rHat) / r^2
        # F = q2 V2 x B
        # F = q2 V2 X (1e-7 q1 v1 x rHat) / r^2
        # rHat is the unit vector pointing from p to self

        B = self.magneticField(p)

        # print "mag force"
        # print " B is", B

        relativeV = self.subtract(self.V(), p.V())

        F = self.product(self.charge, self.cross(relativeV, B))

        return F

    def magnetic_force_total(self, p, dt=None):
        return self.magneticForceTotal5(p, dt)

        # Total force on self, caused by relative self and p velocity
        # Combine the calculation of total force into one simpler formula.
        # Force on p is the negative of this

        # Becomes:
        # F = u0/4pi q1 q2 Vx(Vxr) / r^3

        r2, l2 = self.distance2(p)

        if r2 == 0.0:
            # Should be infinity I guess?
            return 0.0, 0.0, 0.0

        r3 = r2 ** (3.0 / 2.0)

        relativeV = self.subtract(self.V(), p.V())  # order not important

        # r points from self to p
        r = (p.x - self.x, p.y - self.y, p.z - self.z)
        if doMagneticInverse:
            r = self.product(-1.0, r)
        f = self.cross(relativeV, self.cross(relativeV, r))
        F = self.product(1e-7 * self.charge * p.charge / r3, f)

        return F

    def magnetic_force_total2(self, p):
        # Second version for testing new ideas
        # Lets try new idea.  Mag force acts
        # to slow down relative motion parallel
        # to r vector. No what what the sign of
        # the charges, the force is always acting
        # against forward motion.
        # It's equal and opposite on the two particles.
        # so what we calculate for one is just the
        # inverse of what it is for the other.
        #  What force we add to make it slow down
        # We must also take away from ES force
        # to keep total force on the particle the same,
        # after we sum ES and MAG.  So this MAG effect
        # is really just twisting the ES force field.
        # AKA DISTORTING THE ES FORCE FIELD? :)
        # Maybe that's relativity right there???
        # By twisting instead of adding new forces,
        # we mantain potential vs ke energy truth.
        # We aren't creating a new potential energy
        # force field to store mag energy.
        # Idea, is that when two particles pass
        # Parallel to each other, at the speed of light
        # The force acting to slow it down will be the
        # full ES force and the old ES force will be
        # down to zero.  So we will have twisted the
        # ES field a full 90 deg at that point.

        # So to code, find magnatde of v perperndicular to r.
        # is that v dot rHat?
        # Adjust as precent relative ot the speed of light.
        # Divide by v/c.  or v^2/c^2?

        relativeV = self.subtract(self.V(), p.V())  # order not important
        r = (p.x - self.x, p.y - self.y, p.z - self.z)

        # Total force on self, due to other particle p

        # Becomes:
        # F = u0/4pi q1 q2 Vx(Vxr) / r^3

        (r2, l2) = self.distance2(p)

        if r2 == 0.0:
            # Should be infinity I guess?
            return (0.0, 0.0, 0.0)

        r3 = r2 ** (3.0 / 2.0)

        relativeV = self.subtract(self.V(), p.V())  # order not important

        if True:
            # new ES twist to slow down logic sort of.
            # r points from self to p

            r = (p.x - self.x, p.y - self.y, p.z - self.z)
            rHat = self.product(1.0 / magnitude(r), r)
            mag = self.cross(self.cross(relativeV, rHat), rHat)

            # Ok, fucking messy with abs() and magnatude(V) to get v^2.
            # But leave it for now.  Try it, and if it looks good, look for way
            # to simplify the math later.

            magF = self.product(
                1e-7 * abs(self.charge * p.charge) * magnitude(relativeV) / r2,
                mag)

            # So, magF if I coded it correctly, is perpendicular to r, in plane of V.
            # acting against velociety to slow it down in the direction perpendicular to r.
            # But also need to reduce ES by this same magnatude to keep total force the same.
            # OHHH NOO. That's not that simple.  That won't keep total the same!
            # It won't keep the the magnatude of the force vector the same.
            # OK -- decided to punt and try it anyway.  These are v^2 terms, so maybe just maybe
            # there is justification for doing it this way?
            # So, magF, is the force vector, 90 deg to r.  now I have to reduce r by this same amount.
            # So I have to add a vector of this same length, pointing in rHat direction.

            esAdjustF = self.product(magnitude(magF), rHat)
            if self.charge * p.charge < 0:
                esAdjustF = self.product(-1.0, esAdjustF)
            F = vectorSum(magF, esAdjustF)

            if False:
                print("TotalForce2")
                print(" self.V() is ", self.V(), magnitude(self.V()))
                print(" p.V() is    ", p.V(), magnitude(p.V()))
                print(" relativeV is", relativeV, magnitude(relativeV))
                print(" r is        ", r)
                print(" rhat is     ", rHat, "magnatude", magnitude(rHat))
                print(" magF is     ", magF, "magnatude", magnitude(magF))
                print(" esAdjustF is", esAdjustF, "magnatude",
                      magnitude(esAdjustF))
                print(" F is        ", F, "magnatude", magnitude(F))
                print()
        else:
            # Old first way
            # r points from self to p
            r = (p.x - self.x, p.y - self.y, p.z - self.z)

            if doMagneticInverse:
                r = self.product(-1.0, r)
            # f = self.cross(self.V(), self.cross(p.V(), r))
            f = self.cross(relativeV, self.cross(relativeV, r))
            F = self.product(1e-7 * self.charge * p.charge / r3, f)

        return F

    def magneticForceTotal3(self, p, dt=None):
        # Thrid verson for testing new ideas for magnetic force

        # -- THIS IS BORKEN. Tried to do what I'm doing in Total 4
        # then tried to turn it into DF/dt -- but failed to calcuate
        # that correctly and now I've given up.  Because the code
        # below is not getting the signs right becaues I'm not
        # calculating F using vectors and the derivivae of F is
        # not using vectors correctly!
        # So I'm giving up beuase I don't like this idea anyhow!
        # Leving it becuase I might come back to it some day.

        # dF/dt = (- c^2 1e7 q1 q2 / 2 r^3) v dot r
        # Error!  Should have been:
        # dF/dt = (-2 c^2 1e7 q1 q2 / r^3) v dot rHat

        (r2, l2) = self.distance2(p)

        if r2 == 0.0:
            # Should be infinity I guess?
            return (0.0, 0.0, 0.0)

        # relativeV = self.subtract(self.V(), p.V())
        relativeV = self.subtract(p.V(), self.V())

        r3 = r2 ** (3.0 / 2.0)

        r = (p.x - self.x, p.y - self.y, p.z - self.z)
        rHat = self.product(-1.0 / magnitude(r), r)
        # rHat points from p back to self
        # rHat * self.q * self.q must define force on self

        # dr = self.dot(relativeV, rHat)
        # df = -2.0 * self.ke * self.charge * p.charge * dr / r3
        # magFactor = 1.0/(self.c*self.c)
        # F = self.product(df*magFactor, rHat)

        dr = relativeV
        factor = -2.0 * self.ke * self.charge * p.charge / r3
        df = self.product(factor, dr)

        magFactor = 1.0 / (self.c * self.c)
        F = self.product(magFactor, df)

        if True:
            print("ForceTotal3")
            print(" self.V() is ", self.V(), magnitude(self.V()))
            print(" p.V() is    ", p.V(), magnitude(p.V()))
            print(" relativeV is", relativeV, magnitude(relativeV))
            print(" r is        ", r)
            print(" rhat is     ", rHat, magnitude(rHat))
            es = self.esForce(p)
            print(" esForce is  ", es, magnitude(es))
            if dt is not None:
                print(" dt is       ", dt)
            print(" dr/dt is    ", dr)
            if dt is not None:
                print(" dr/dt * dt  ", self.product(dt, dr))
                print(" dr/dt * dt m", magnitude(dr) * dt)
            print(" df/dt is    ", df)
            if dt is not None:
                print(" df/dt * dt  ", self.product(dt, df))
                print(" df/dt * dt m", magnitude(df) * dt)
            print(" F is        ", F, magnitude(F))
            print()

        return F

    def magneticForceTotal4(self, p, dt=None):
        # Forth verson for testing new ideas for magnetic force
        # This is what I tried to code with Total3 but then
        # relized was thinking about it all wrong.  I
        # want to calucate how the strenght of E is changing
        # over time, and then modify that strenth, but not
        # the direction of of the force vector, as a a function
        # of how fast the strength of E is changing! aka
        # d(magnatude(e))/dt.

        # Calculate D|e|/dt and make magenetic force
        # Nothing more than an increase in the ES force based
        # on D|e|/dt.
        # This is a very different approach, but yet seems like
        # it could produces
        # the same sort of results for current -- aka the
        # sum of many electons moving at once.  But it has
        # the intersting effect of having no magnetic force
        # at all when the two particles pass side by side
        # becuase de/dt = 0 at that point!  Which could be
        # useful in not messing up orbits!
        # And it might force eleptical high speed orbits into
        # circles.

        # d|F|/dt = (-2 c^2 1e7 q1 q2 / r^3) v dot rHat

        (r2, l2) = self.distance2(p)

        if r2 == 0.0:
            # Should be infinity I guess?
            return (0.0, 0.0, 0.0)

        # relativeV = self.subtract(p.V(), self.V())
        relativeV = self.subtract(self.V(), p.V())

        r3 = r2 ** (3.0 / 2.0)

        r = (p.x - self.x, p.y - self.y, p.z - self.z)
        rHat = self.product(-1.0 / magnitude(r), r)
        # rHat points from p back to self
        # rHat * self.q * self.q must define force on self

        # dr = self.dot(relativeV, rHat)
        # df = -2.0 * self.ke * self.charge * p.charge * dr / r3
        # magFactor = 1.0/(self.c*self.c)
        # F = self.product(df*magFactor, rHat)

        es = self.esForce(p)
        dr = self.dot(relativeV, rHat)
        beta2 = (dr / self.c) ** 2
        F = self.product(beta2, es)

        if False:
            print("ForceTotal4")
            print(" self.V() is ", self.V(), magnitude(self.V()))
            print(" p.V() is    ", p.V(), magnitude(p.V()))
            print(" relativeV is", relativeV, magnitude(relativeV))
            print(" r is        ", r)
            print(" rhat is     ", rHat, magnitude(rHat))
            es = self.esForce(p)
            print(" esForce is  ", es, magnitude(es))
            if dt is not None:
                print(" dt is       ", dt)
            print(" dr/dt is    ", dr)
            if dt is not None:
                print(" dr/dt * dt  ", dr * dt)
            print(" df/dt is    ", df)
            if dt is not None:
                print(" df/dt * dt  ", df * dt)
            print(" F is        ", F, magnitude(F))
            print()

        return F

    def magneticForceTotal5(self, p, dt=None):
        # Fith version for testing yet another magnatism
        # idea.  This one is based on the idea that magnatism
        # needs to translate the force field into kenetic
        # energy in a separate dimention from E.  But
        # when in an orbit, the idea is not to speed or slow
        # down the orbit -- that's the dimension E is using
        # store ke in the orbit.  And not to make the orbting
        # particle go into a higher or lower orbit -- again
        # that's the same dimension as E is using. But instead
        # Make it turn, say from an equator orbit, to a polar
        # orbit!  So the turn doesn't mess up the kenetic energy
        # of the E orbit at all!  But instead adds NEW kenetic
        # energy of a complex sprial or spinning orbit, or something
        # complex I can't even grasp at the moment!  But a separate
        # kenetic energy store for the E orbit!

        # So, to do this.  We apply the force of the B field
        # at 90 degree to R and 90 deg to (relaveV)! aka vxr.
        # AKA, what B flux lines are, but not how B force is
        # described!

        # But, we are going to guess that the mangnatude of
        # B should be v^2/c^2 relative to E strength.  Not
        # with V reduced by vxrHat. (aka using it's full
        # velociety not it's sin theta reduced velociety.
        # EDIT -- coded vxr instead. It's easier. Cleaner.

        # And it must be equal and opposit of course.
        # so p1.mag(p2) is the inverse of p2.mag(p1).

        # d|F|/dt = (-2 c^2 1e7 q1 q2 / r^3) v dot rHat

        (r2, l2) = self.distance2(p)

        if r2 == 0.0:
            # Should be infinity I guess?
            return (0.0, 0.0, 0.0)

        # relativeV = self.subtract(p.V(), self.V())
        # relativeV = self.subtract(self.V(), p.V())

        # r3 = r2 ** (3.0/2.0)

        # r = (p.x - self.x, p.y - self.y, p.z - self.z)
        # rHat = self.product(-1.0/magnatude(r), r)

        # rHat points from p back to self
        # rHat * self.q * self.q must define sign of force on self?

        if False:
            # Old attempt
            b = self.cross(relativeV, r)
            factor = 1e-7 * self.charge * p.charge * magnitude(relativeV) / r3
            if self.charge < 0:
                factor *= -1.0  # well, needed to make ep pair push opposit directions!
            F = self.product(factor, b)

        # So it's 1e-7 q * q * v * v / r^2 in total whcih gets all the unints
        # consitent with E. But in the direction of v x r.

        if False:
            # Try making B driven by dE/dt. The faster E is changing, the
            # strong the B field will be.
            # Keep B at 90 deg to E.  Seems to be required.  But
            # we can point it either at vxr, so it's perpendicular to
            # r AND v. Or we can be in the same plane with v and r, while
            # being perpendicular to R still.
            # Lets try vxr first.
            # And I really mean d|E|/dt, not dE/dt.
            #### NOTE: from above: d|F|/dt = (-2 c^2 1e7 q1 q2 / r^3) v dot rHat
            dr = self.dot(relativeV, rHat)
            # df = -2.0 * self.ke * self.charge * p.charge * dr / r3
            # Figured out how to convert Dr into Force!  Stupid me for not
            # seeing this before.  Dr is V!  And we know the formula is v^2!
            df = self.ke * self.charge * p.charge * dr * dr / r2
            #### [[[ error above, should be 1e-7 not self.ke?? ]]]]
            b = self.cross(relativeV, rHat)
            F = self.product(df, b)
            if True:
                print("ForceTotal5")
                print(" self.V() is ", self.V(), magnitude(self.V()))
                print(" p.V() is    ", p.V(), magnitude(p.V()))
                print(" relativeV is", relativeV, magnitude(relativeV))
                print(" r is        ", r)
                print(" rhat is     ", rHat, magnitude(rHat))
                es = self.esForce(p)
                print(" esForce is  ", es, magnitude(es))
                print(" dr is       ", dr)
                print(" df is       ", df)
                print(" b is        ", b, magnitude(b))
                print(" F is        ", F, magnitude(F))

        if True:
            # Total5(C) -- third attempt coded in total5.
            # Ok, figured out we can not use vxr as the direction
            # of B. It makes equal and opposit impossible for ++ and --
            # paris.  But making V in the plane with v and r, while
            # 90 deg to E is possible.  So lets do that!  We already
            # tried a version of this in Total2()  But we return.
            # The difference is that we use V dot r, instead of V x r
            # this time.  That means this is maximal B when the two
            # particles are closing on each other the fastest - which
            # makes them serve away BTW and not run into each other.
            # I would think.

            # Think of p is being at the orign standing still.
            # self particles is elsewhere and moving.

            relativeV = self.subtract(self.V(), p.V())

            # Sign of relativeV is of the velcoity of self if the
            # the velociaty of p is zero.

            r = (self.x - p.x, self.y - p.y, self.z - p.z)
            rHat = self.product(1.0 / magnitude(r), r)

            # r points from p (logically at orign) to self.

            # Magnatude of v in line with r
            vr = self.dot(relativeV, rHat)

            vHat = self.product(1.0 / magnitude(relativeV), relativeV)
            bHat = self.cross(self.cross(vHat, rHat), rHat)
            bMag = magnitude(bHat)
            if bMag != 0.0:
                bHat = self.product(1.0 / bMag, bHat)
                # Otherwise, leave it as 0,0,0 and let it
                # return F of zero.
                # This is a problem.  It happens because we don't
                # know which direction to point B!  R and V are
                # in the same, so we don't know which direction
                # is both 90 deg to R and in the same plane as V.
                # But yet, at this point, B should be maximal value!
                # This brings this whole idea under question as to
                # whether this is logically valid to begin with.
            factor = 1.0e-7 * self.charge * p.charge * vr * vr / r2
            F = self.product(factor, bHat)
            # F = self.product(-1.0, F) # backwards

            if False:
                print("ForceTotal5")
                print("self x y is  ", self.x / Angstrom, self.y / Angstrom)
                print("P   x y is   ", p.x / Angstrom, p.y / Angstrom)
                print(" self.V() is ", self.V(), magnitude(self.V()) / self.c)
                print(" p.V() is    ", p.V(), magnitude(p.V()) / self.c)
                print(" relativeV is", relativeV, magnitude(relativeV))
                print(" VHat is     ", vHat, magnitude(vHat))
                print(" r is        ", r)
                print(" rhat is     ", rHat, magnitude(rHat))
                es = self.esForce(p)
                print(" esForce is  ", es, magnitude(es))
                print(" vr is       ", vr)
                print(" bHat is     ", bHat, magnitude(bHat))
                print(" factor is   ", factor)
                print(" F is        ", F, magnitude(F))
                print(" F dot rHat  ", self.dot(F, rHat))

        if False:
            # Hard Code fake mag force at 90 deg to v and r for an electon
            # only, at 1/2 the force of E.
            F = (0.0, 0.0, 0.0)
            if isinstance(self, Electron):
                eForce = abs(self.ke * self.charge * p.charge / r2)
                vHat = self.product(1.0 / magnitude(relativeV), relativeV)
                b = self.cross(vHat, rHat)
                bHat = self.product(1.0 / magnitude(b), b)
                F = self.product(eForce, bHat)
                # print "v is", self.V()
                # print "rhat is", rHat
                # print "bHat is", bHat
                # print "eForce is", eForce
                # es = self.esForce(p)
                # print " esForce is  ", es, magnatude(es)
                # print "F    is", F
                # sys.exit(1)

        return F

    def addEndForce(self, p):
        if p is self:
            return

        dx = (self.endX - p.endX)
        dy = (self.endY - p.endY)
        dz = (self.endZ - p.endZ)

        r2, l2 = self.endDistance2(p)

        force = self.ke * (self.charge * p.charge) / l2

        r = math.sqrt(r2)

        self.endFx += force * dx / r
        self.endFy += force * dy / r
        self.endFz += force * dz / r

        if doMagnetic:
            f = self.magnetic_force_total(p)
            self.endFx += f[0]
            self.endFy += f[1]
            self.endFz += f[2]

    def calculateEndVelocity(self, dt):
        # Assume linear change in acceleration from start (fx to end endFx)
        # (yes this is the correct intergral of a linear chagne in acceleration)
        # I had to recalcuate it 10-4-2016 to verify

        self.endVx = self.vx + ((self.fx + self.endFx) / 2.0) * dt / self.mass
        self.endVy = self.vy + ((self.fy + self.endFy) / 2.0) * dt / self.mass
        self.endVz = self.vz + ((self.fz + self.endFz) / 2.0) * dt / self.mass

    def calculateEndPosition(self, dt):  # Calcluate end X Y Z
        # Calculate ending position using both start and end velocity

        # Assume force (acceleration) changes but is linear from start to end
        # x = 1/2 a t^2 + vt + x
        # a is (2*as + ae)/3  --- where as is a start, and ae is a end.
        # endx = 1/2 (2as + ae)/3 t^2 + v t + x

        # self.endX = self.x + self.vx * dt + 0.5 * (2.0*self.fx + self.endFx)/(3.0*self.mass) * dt ** 2.0

        self.endX = self.x + self.vx * dt + (self.fx + 0.5 * self.endFx) / (
                    3.0 * self.mass) * dt ** 2.0
        self.endY = self.y + self.vy * dt + (self.fy + 0.5 * self.endFy) / (
                    3.0 * self.mass) * dt ** 2.0
        self.endZ = self.z + self.vz * dt + (self.fz + 0.5 * self.endFz) / (
                    3.0 * self.mass) * dt ** 2.0

    def move(self):
        # Make end state the current state
        # Save current state as old state
        # Leaves ending force the same as the starting

        self.x = self.endX
        self.y = self.endY
        self.z = self.endZ

        self.vx = self.endVx
        self.vy = self.endVy
        self.vz = self.endVz

        self.fx = self.endFx
        self.fy = self.endFy
        self.fz = self.endFz

    def resetState(self):
        # Reset so we can recompute step with new DT
        # Reset force end to match force begin
        # Evertyhing else will be recomputed again.

        self.endFx = self.fx
        self.endFy = self.fy
        self.endFz = self.fz

    def keneticEnergy(self):
        # print "KE CALC BEGIN"
        return self.keneticEnergyCalc(self.vx, self.vy, self.vz)

    def keneticEndEnergy(self):
        # print "KE CALC END"
        return self.keneticEnergyCalc(self.endVx, self.endVy, self.endVz)

    def keneticEnergyCalc(self, vx, vy, vz):
        ## 1/2 m v**2
        ke = 0.5 * self.mass * (vx ** 2.0 + vy ** 2.0 + vz ** 2.0)
        # print "KE CALC vx,vy,vz:", vx, vy, vz
        # print "KE CALC answer =", ke
        return ke

    def setKeneticEnergy(self, ke):
        # Back calculate velocity using given ke -- keep direction the same
        newV2 = ke / (0.5 * self.mass)
        oldV2 = (self.vx ** 2.0 + self.vy ** 2.0 + self.vz ** 2.0)
        # print "in set kenetic newV2 is", newV2
        # print "in set kenetic oldV2 is", oldV2
        newV = math.sqrt(newV2)
        oldV = math.sqrt(oldV2)
        # self.vx *= newV2 / oldV2
        # self.vy *= newV2 / oldV2
        # self.vz *= newV2 / oldV2
        self.vx *= newV / oldV
        self.vy *= newV / oldV
        self.vz *= newV / oldV

    def distance2(self, p):  # distance squared

        if p is self:
            return 0.0

        dx = (self.x - p.x)
        dy = (self.y - p.y)
        dz = (self.z - p.z)

        d2 = dx ** 2.0 + dy ** 2.0 + dz ** 2.0

        return self.limitedDistance2(d2)

    def endDistance2(self, p):  # distance squared

        if p is self:
            return 0.0

        dx = (self.endX - p.endX)
        dy = (self.endY - p.endY)
        dz = (self.endZ - p.endZ)

        d2 = dx ** 2.0 + dy ** 2.0 + dz ** 2.0

        return self.limitedDistance2(d2)

    def limitedDistance2(self, d2):

        # Limit distance to RLimit to solve computational problems.
        # return (real, limited) tuple

        return (d2, max(d2, RLimit ** 2))

    def distance(self, p):
        r, l = self.distance2(p)
        return (math.sqrt(r), math.sqrt(l))

    def endDistance(self, p):
        r, l = self.endDistance2(p)
        return (math.sqrt(r), math.sqrt(l))

    def potentialEnergy(self, p):
        # potential energy between self and particle P

        if p is self:
            return 0.0  # Bogus should be +infinity

        r, l = self.distance(p)

        return self.potentialEnergyForDistance(p, r)

    def potentialEndEnergy(self, p):
        # potential energy between self and particle P

        if p is self:
            return 0.0  # Bogus should be +infinity

        r, l = self.endDistance(p)

        return self.potentialEnergyForDistance(p, r)

    def potentialEnergyForDistance(self, p, d):

        # if d == 0:
        # return 0.0	# Bogus should be +infinity

        if d >= RLimit:
            return self.ke * self.charge * p.charge / d

        global InsideRLimitCount
        InsideRLimitCount += 1

        x = self.ke * self.charge * p.charge / RLimit

        return x + (x / RLimit) * (RLimit - d)

        # self.charge = -1.60218e-19 # in Coulombs
        # self.ke = 8.9875517873681764e9  # Coulomb's constant (1/4 pi e) written as K(sub)e

    def momentum(self):
        # Return (mx,my,mz) tuple
        return self.mass * self.vx, self.mass * self.vy, self.mass * self.vz

    def addMomentum(self, m):
        self.vx += m[0] / self.mass
        self.vy += m[1] / self.mass
        self.vz += m[2] / self.mass

    def addVelocity(self, v):
        self.vx += v[0]
        self.vy += v[1]
        self.vz += v[2]


class Electron(Particle):
    def __init__(self, x=0.0, y=0.0, z=0.0):
        Particle.__init__(self, x, y, z)
        self.charge = -1.60218e-19  # in Coulombs
        self.mass = 9.10938356e-31
        self.symbol = 'e'
        if FakeConstants:
            self.charge = -1.0
            self.mass = 1.0


class Proton(Particle):
    def __init__(self, x=0.0, y=0.0, z=0.0, n=1.0):
        Particle.__init__(self, x, y, z)
        self.charge = n * +1.60218e-19  # in Coulombs
        self.mass = n * 1.672621898e-27
        self.symbol = 'p'
        if FakeConstants:
            self.charge = +1.0
            # Keep real ratio the same -- mass of electon is 1.0
            # mass of proton is 1836.15267376
            self.mass = 1.672621898e-27 / 9.10938356e-31  # units for mass e == 1.0


global ResetEnergy
ResetEnergy = False


class ParticleImage:
    def __init__(self, p):
        self.p = p  # Particle

    def display(self, zMin, zMax):

        x = self.spaceToPixels(self.p.x)
        y = self.spaceToPixels(self.p.y)
        z = self.spaceToPixels(self.p.z)

        if isinstance(self.p, Electron):
            color = black
            size = 2
        else:
            color = red
            size = 4

        minScale = 1.0
        maxScale = 3.0
        size *= (z / float(screen_depth)) * (maxScale - minScale) + minScale
        size = abs(int(size))

        inset = 10
        if isinstance(self.p, Proton):
            inset = 40

        # eChange = 0.25
        eChange = 1.00

        global ResetEnergy
        bounce = False
        if x < inset and self.p.vx < 0:
            self.p.vx *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.x = self.pixelsToSpace(inset)
            bounce = True
        if x > screen_width - inset and self.p.vx > 0:
            self.p.vx *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.x = self.pixelsToSpace(screen_width - inset)
            bounce = True
        if y < inset and self.p.vy < 0:
            self.p.vy *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.y = self.pixelsToSpace(inset)
            bounce = True
        if y > screen_height - inset and self.p.vy > 0:
            self.p.vy *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.y = self.pixelsToSpace(screen_height - inset)
            bounce = True
        if z < inset and self.p.vz < 0:
            self.p.vz *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.z = self.pixelsToSpace(inset)
            bounce = True
        if z > screen_depth - inset and self.p.vz > 0:
            self.p.vz *= -1
            self.p.setKeneticEnergy(self.p.keneticEnergy() * eChange)
            ResetEnergy = True
            self.p.z = self.pixelsToSpace(screen_depth - inset)
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
        # So setting eChange = 1.0 does not conserve energy correctly.
        # We adjust the systems total target for fixing energy
        # on a bounce which forces us to accept this broken
        # energy total as correct.
        # But if we turn off the Bounce check, that allows the energyFix
        # system to keep fixing this error for us.  Hince, the test
        # below to turn Bounce flag off if eChange is 1.00

        if eChange == 1.0:
            ResetEnergy = False  # Don't reset energy total -- fix it instead

        x = self.spaceToPixels(self.p.x)
        y = self.spaceToPixels(self.p.y)

        # print "x y is", x, y
        pygame.draw.circle(screen, color, (x, y), size, 0)

    def spaceToPixels(self, space):
        # 0,0 is the same in both and is the top left corner of the screen
        return int(pixelsPerAngstrom * space / Angstrom)

    def pixelsToSpace(self, pixels):
        return pixels * Angstrom / pixelsPerAngstrom


def vectorSum(a, b):
    # breaks if vectors not the same length
    l = max(len(a), len(b))
    s = [0.0] * l
    for i in range(l):
        s[i] = a[i] + b[i]
    return tuple(s)


def totalMomentum(world):
    s = (0.0, 0.0, 0.0)
    for i in range(len(world)):
        s = vectorSum(s, world[i].momentum())
    return s


def total_kinetic_energy(world):
    totalKE = 0.0

    for i in range(len(world)):
        totalKE += world[i].keneticEnergy()

    return totalKE


def total_potential_energy(world):
    total_pe = 0.0

    for i in range(len(world)):
        for j in range(i + 1, len(world)):
            total_pe += world[i].potentialEnergy(world[j])

    return total_pe


def total_end_kinetic_energy(world):
    total_ke = 0.0

    for i in range(len(world)):
        total_ke += world[i].keneticEndEnergy()

    return total_ke


def total_end_potential_energy(world):
    total_pe = 0.0

    for i in range(len(world)):
        for j in range(i + 1, len(world)):
            total_pe += world[i].potentialEndEnergy(world[j])

    return total_pe


#####################################################################
# Normalize momentum
#####################################################################

def zero_momentum(world):
    """ Normalize momentum to zero to freeze frame of reference. """
    tm = totalMomentum(world)
    n = len(world)
    total_mass = 0.0
    for i in range(n):
        total_mass += world[i].mass
    fudge = (-tm[0] / total_mass, -tm[1] / total_mass, -tm[2] / total_mass)
    for i in range(n):
        world[i].addVelocity(fudge)


#####################################################################
# Move objects to place center of mass in the middle of the screen
#####################################################################

def center_mass(world):
    center_x = (screen_width / pixelsPerAngstrom * Angstrom) / 2.0
    center_y = (screen_height / pixelsPerAngstrom * Angstrom) / 2.0
    center_z = (screen_depth / pixelsPerAngstrom * Angstrom) / 2.0
    cx = cy = cz = 0.0
    tm = 0.0
    for p in world:
        cx += p.x * p.mass
        cy += p.y * p.mass
        cz += p.z * p.mass
        tm += p.mass
    for p in world:
        p.x = p.x - cx / tm + center_x
        p.y = p.y - cy / tm + center_y
        p.z = p.z - cz / tm + center_z


def center_of_mass(world):
    cx = cy = cz = 0.0
    tm = 0.0
    for p in world:
        cx += p.x * p.mass
        cy += p.y * p.mass
        cz += p.z * p.mass
        tm += p.mass
    x = cx / tm
    y = cy / tm
    z = cz / tm
    return x, y, z


# fastTest()
# neutronGravityTest()
# magneticTest()
# forceCircleTest()
# mag_circle_test()

screen = pygame.display.set_mode(screen_size)
clock = pygame.time.Clock()
pygame.display.set_caption('Physics')

p1 = Proton(Angstrom * 0.0, 0.0, 0.0 * Angstrom, n=10.0)
e1 = Electron(Angstrom * 0.25, 0.0, 0.0 * Angstrom)
e2 = Electron(Angstrom * -0.25, 0.0 * Angstrom, 0.0 * Angstrom)

e1.vy = 200000.0
e1.vy = 100000.0
e1.vy = 0.0
e1.vy = 1600000.0
e1.vy = 16000000.0
# e1.vx = 800000.0
e1.vy = 3000000.0

# This is what the below calculates for the perfect circular orbit
# Velocity for the ep pair spaced at 0.25A
e1.vy = 3181993.44316
p1.vy = -3181993.44316 * e1.mass / p1.mass
e2.vy = -3181993.44316
e2.vz = -3181993.44316 / 10
# e1.vy = p1.vy = 0.0
# e1.vz = 800000.0

if False:  # I have no clue what this code does 5-11-2018 CW
    r = Angstrom * 0.25
    print("r is", r)
    r = e1.x - p1.x
    print("r is", r)
    rCenter = r * p1.mass / (p1.mass + e1.mass)
    e1.vz = e1.c * math.sqrt(
        abs(1e-7 * e1.charge * p1.charge * rCenter / (r * r * e1.mass)))
    rCenter = r * e1.mass / (p1.mass + e1.mass)
    p1.vz = -p1.c * math.sqrt(
        abs(1e-7 * e1.charge * p1.charge * rCenter / (r * r * p1.mass)))
    print("p1.vz is", p1.vz, "momentum is", p1.momentum())
    print("e1.vz is", e1.vz, "momentum is", e1.momentum())
    print("p1.vz - e1.vz", p1.vz - e1.vz)

    p1.vy = p1.vz
    e1.vy = e1.vz
    p1.vz = 0.0
    e1.vz = 0.0
    e1.vy /= 2.0

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

p3 = Proton(Angstrom * 6.0, Angstrom * 0.0, 0.0)
e3 = Electron(Angstrom * 7.0, Angstrom * 0.0, 0.0)
e3.vy = 800000.0
e3.vy = 1600000.0

piworld = []

if True:
    pi1 = ParticleImage(p1)
    piworld.append(pi1)
    pi2 = ParticleImage(e1)
    piworld.append(pi2)
    # pi3 = ParticleImage(p2)
    # piworld.append(pi3)
    pi4 = ParticleImage(e2)
    piworld.append(pi4)

if False:
    pi3 = ParticleImage(p2)
    piworld.append(pi3)
    if False:
        pi4 = ParticleImage(e2)
        piworld.append(pi4)

if False:
    pi5 = ParticleImage(p3)
    piworld.append(pi5)
    pi6 = ParticleImage(e3)
    piworld.append(pi6)

if False:  # Two electrons
    e1 = Electron(Angstrom * 0.1, Angstrom * 1.0, 0.0)
    piworld.append(ParticleImage(e1))
    e2 = Electron(Angstrom * 5.1, Angstrom * 1.0, 0.0)
    piworld.append(ParticleImage(e2))

if False:
    # Random particles

    p = 2
    e = 2
    percentOfScreen = 0.10
    e1 = Electron()  # Just need any particle
    monumtum = e1.mass * e1.c / 1000.0
    world = []
    xwidth = (screen_width * percentOfScreen / pixelsPerAngstrom) * Angstrom
    ywidth = (screen_height * percentOfScreen / pixelsPerAngstrom) * Angstrom

    for i in range(p):
        x = random.random() * xwidth
        y = random.random() * ywidth
        z = random.random() * min(xwidth, ywidth)
        # z = 0.0
        p = Proton(x, y, z)
        world.append(p)
        p.vx = random.random() * monumtum / p.mass
        p.vy = random.random() * monumtum / p.mass
        p.vz = random.random() * monumtum / p.mass

    for i in range(e):
        x = random.random() * xwidth
        y = random.random() * ywidth
        z = random.random() * min(xwidth, ywidth)
        # z = 0.0
        p = Electron(x, y, z)
        world.append(p)
        p.vx = random.random() * monumtum / p.mass
        p.vy = random.random() * monumtum / p.mass
        p.vz = random.random() * monumtum / p.mass

    for p in world:
        piworld.append(ParticleImage(p))

world = []
for p in piworld:
    world.append(p.p)

lastKE = 0.0
lastPE = 0.0

zero_momentum(world)
center_mass(world)
# e2.x = 40 * Angstrom

#####################################################################
# Calculate all initial forces
# Sets starting and ending forces at the same time
#####################################################################

for p1 in world:
    p1.zeroForce()
    for p2 in world:
        p1.addForce(p2)

startingTotalKE = total_kinetic_energy(world)
startingTotalPE = total_potential_energy(world)

#####################################################################
# Main draw and simulation loop
#####################################################################

fps_limit = 60
run_me = True
i = 0
cycleCount = 0.0
now = 0.0
energyDiffMax = None
maxVc = 0.0

lastD = 0.0
lastNow = now
lastV = 0.0
lastA = 0.0

while run_me:
    clock.tick(fps_limit)
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run_me = False

    screen.fill(white)

    zMin = None
    zMax = None
    for p in piworld:
        z = p.p.z
        if zMin is None or z < zMin:
            zMin = z
        if zMax is None or z > zMax:
            zMax = z

    slist = sorted(piworld, key=attrgetter('p.z'))
    for p in slist:
        p.display(zMin, zMax)

    if False:
        for p in piworld:
            if isinstance(p.p, Proton):
                p.display(zMin, zMax)
        for p in piworld:
            if isinstance(p.p, Electron):
                p.display(zMin, zMax)

    pygame.display.flip()

    ###
    #  Start of physics simulation step loop
    ###

    crtClearAndHome()

    print()
    print("Time now is %25.40f" % now)
    print("Total Momentum %8.1e %8.1e %8.1e" % totalMomentum(world))
    print()
    print("doMagnetic:", doMagnetic, "  doMagneticInverse:", doMagneticInverse)
    print()

    totalKE = total_kinetic_energy(world)
    totalPE = total_potential_energy(world)

    if ResetEnergy:
        # Bounce happened on display above
        # Bounce will remove energy, we need to cope...
        startingTotalKE = totalKE
        startingTotalPE = totalPE
        ResetEnergy = False

    re = ((totalKE - startingTotalKE) + (totalPE - startingTotalPE))

    print("Total ke     %8.1e" % totalKE,
          "change from last %8.1e" % (totalKE - lastKE))
    print("Total pe     %8.1e" % totalPE,
          "change from last %8.1e" % (totalPE - lastPE))
    print("Total    e:  %8.1e" % (totalKE + totalPE))
    print("Relative e:  %8.1e" % (re), end=' ')
    print("  %% of Total ke:%8.4f" % (abs(re * 100.0 / totalKE)))
    print()
    print("DT is: %4.1e" % dt)
    print()
    print("PixelsPerAngstrom", pixelsPerAngstrom, end=' ')
    print("Screen (%.1fx%.1f) A" % (
        screen_width / pixelsPerAngstrom, screen_height / pixelsPerAngstrom))
    print()

    print("Inside Rlimit ", InsideRLimitCount, " pBounce:", pBounceCount,
          " eBounce:", eBounceCount)
    print()

    lastKE = totalKE
    lastPE = totalPE

    cnt = 0
    for i in range(len(world)):
        p1 = world[i]
        if not isinstance(p1, Proton):
            continue
        for j in range(i + 1, len(world)):
            p2 = world[j]
            if not isinstance(p2, Proton):
                continue
            print("Distance between", end=' ')
            print("P%02d P%02d %11.6f" % (i, j, p1.distance(p2)[0] / Angstrom))
            cnt += 1
            if cnt > 10:
                break
        if cnt > 10:
            break
    print()

    ###########################################
    # Debug print out particle stats
    ###########################################

    totalAvgKE = 0.0
    for i in range(len(world)):
        p1 = world[i]
        p1.avgKE += (p1.keneticEnergy() - p1.avgKE) * 0.0001
        totalAvgKE += p1.avgKE

    for i in range(len(world)):
        p1 = world[i]
        sys.stdout.write(p1.symbol)
        vc = magnitude(p1.V()) / p1.c
        maxVc = max(maxVc, vc)
        print("%d vx:%10.2e  vy:%10.2e %5.3fc" % (i, p1.vx, p1.vy, vc),
              end=' ')
        print("x:%10.5f A" % (p1.x / Angstrom), end=' ')
        # print " KE:%10.2e" % p1.avgKE
        print(" KE:%6.2f%%" % (p1.avgKE * 100 / totalAvgKE))

    print()
    print("Max Velociety: %5.3fc" % maxVc)
    print()

    ###########################################
    # Debug center of mass of pairs for
    # special e/p pair drift test
    ###########################################

    lastCm = None
    for i in range(len(world) - 1):
        p1 = world[i]
        p2 = world[i + 1]
        if isinstance(p2, Electron) and isinstance(p1, Proton):
            cm = center_of_mass((p1, p2))
            print("Center of Mass of pair %20.15f %20.15f %20.15f" % (
                cm[0] / Angstrom, cm[1] / Angstrom, cm[2] / Angstrom))
            if lastCm is not None:
                d = math.sqrt((cm[0] - lastCm[0]) ** 2 +
                              (cm[1] - lastCm[1]) ** 2 +
                              (cm[2] - lastCm[2]) ** 2)
                print("Distance is %20.15f" % (d / Angstrom))
                dd = d - lastD
                print("Change   is %20.15f" % (dd / Angstrom))
                move_dt = now - lastNow
                print("DT is %20.30f" % move_dt)
                v = 0.0
                if move_dt != 0.0:
                    v = dd / move_dt  # Velocity
                    print("M/S is ", dd / move_dt)
                    a = (v - lastV) / move_dt
                    print("A is ", a)
                    lastA += (a - lastA) * .001
                    print("avg A is", lastA)
                lastV = v
                lastNow = now
                lastD = d
                gravityForce = 0.0
                gravityForce += p1.gravityForce(world[i - 1])
                gravityForce += p1.gravityForce(world[i - 2])
                gravityForce += p2.gravityForce(world[i - 2])
                gravityForce += p2.gravityForce(world[i - 1])
                print("Gravity Force is", gravityForce)
                es = p1.esForce(p2)
                print("es from p1 to p2 is", es)
                es = p1.esForce(world[i - 2])
                print("es from p1 to e2 is", es)
                print("es from 1 to p4 is", es)
                print("Gravity A is", gravityForce / (p1.mass + p2.mass))
            lastCm = cm
    print()

    #####################################################################
    # Calculate next position now
    #####################################################################

    while True:  # DT restart loop

        # At this point, all particles have a known position and velocity
        # and the force is calculated to match the position.  Ending force
        # is also set to match the current force but will be updated.

        for it in range(3):
            # print "begin iteration", it
            for i in range(len(world)):
                p1 = world[i]
                p1.calculateEndVelocity(dt)
                p1.calculateEndPosition(dt)

            # Update ending force based on ending position calculated above.
            # Ending position is a function of ending force so when we update
            # the force, we can then loop back and calculate a new ending
            # Position

            # Could calculate PE in same loop we use for updating forces
            # because it duplicates a lot of the same work, but need to
            # restructure code to do that.

            for p1 in world:
                p1.zeroEndForce()
                for p2 in world:
                    p1.addEndForce(p2)

            # totalKE2 = total_end_kinetic_energy(world)
            # totalPE2 = total_end_potential_energy(world)
            # curEnergy = totalKE2 + totalPE2
            # error = (totalKE - totalKE2) + (totalPE - totalPE2)
            # print "Energy error after iteration     %18.10e" % error

        if energyFix2:
            # Ok, one last time, fudge velocity based on current position
            # and total force for this position
            # This is just always good I think.  Probably shouldn't be an option.
            for i in range(len(world)):
                p1 = world[i]
                p1.calculateEndVelocity(dt)

        # Now calculate total Energy after that would result if we take this move

        totalKE2 = total_end_kinetic_energy(world)
        totalPE2 = total_end_potential_energy(world)

        # energyDiff = (totalKE2 + totalPE2) - (totalKE + totalPE) # last move error only

        # Ah, a computational problem showed up.  When one of the two is many
        # orders of magnitude different from the other, the accuracy is all lost
        # in the difference if do the above way vs the way below!

        energyDiff = (totalKE2 - totalKE) + (
                    totalPE2 - totalPE)  # last move error only
        energyDiffStart = totalKE2 - startingTotalKE  # Error from start
        energyDiffStart += totalPE2 - startingTotalPE

        # print "old KE was", totalKE
        # print "old PE was", totalPE
        # print "new totalKE after move is", totalKE2
        # print "new totalPE after move is", totalPE2
        print("Energy error for move is %12.4e" % energyDiff, end=' ')
        if energyDiffMax is not None:
            print("max error %12.4e" % energyDiffMax, end=' ')
        print()
        # print "energy diff from start is", energyDiffStart

        cycleCount += 1

        # if False and cycleCount > 300:
        # sys.exit(1)

        ######################################################################
        # Dynamically change dt to maintain error and maximise simulation speed
        ######################################################################

        if dtAdjust and totalKE2 != 0.0:
            # print "==DO DTADJUST cycleCount", cycleCount, "abs(energyDiff)", abs(energyDiff),
            # print "totalKE2", totalKE2, "percent", abs(energyDiff) / totalKE2
            if dt < dtMax and cycleCount > 3 and abs(
                    energyDiff) / totalKE2 < 0.0001:
                # print "SPEEDUP -- increase DT abs(diff)/total is", abs(energyDiff) / totalKE2
                dt *= 2.0
                dt = min(dt, dtMax)
                # cycleCount = 0
                # continue
                # No need to restart -- take this move as fine but use a larger dt
                # for the next move

            elif dt > dtMin and abs(energyDiff) / totalKE2 > 0.001:
                # print "SLOWDOWN -- reduce DT abs(diff)/total is", abs(energyDiff) / totalKE2
                dt /= 2.0
                dt = max(dt, dtMin)
                for p1 in world:
                    p1.resetState()
                cycleCount = 0
                continue

        if energyDiffMax is None:
            energyDiffMax = energyDiff
        energyDiffMax = max(energyDiffMax, abs(energyDiff))

        break  # End of DT adjust loop

    #####################################################################
    # Make the move -- move current ending state to be the current
    # Also leaves ending force set to same as starting force
    #####################################################################

    for p1 in world:
        p1.move()

    # print
    # print "MOVE MADE -------------------------------------------"

    now += dt

    ######################################################################
    # Fix total energy error
    # Adjust all velocities to remove energy error of energyDiffStart
    ######################################################################

    if energyFix:
        if energyDiffStart > totalKE2:
            print("energy error greater than total KE error:", energyDiffStart,
                  "total", totalKE2)
            energyDiffStart = totalKE2 * 0.5
            # Only try to take half of total KE out of the system
            # time.sleep(1)
            # sys.exit(1)

        for i in range(len(world)):
            # print
            # print "Adjust KE for particle", i, "energyDiffStart error is", energyDiffStart
            p1 = world[i]
            ke = p1.keneticEnergy()
            # print "Current ke is", ke
            # print "percent of of total is", ke/totalKE2*100.0, "%"
            # print "amount added of total energyDiffStart is", ke/totalKE2*energyDiffStart
            newKe = ke - ke / totalKE2 * energyDiffStart
            # print "new should be", newKe
            p1.setKeneticEnergy(newKe)
            # print "new KE is now", p1.keneticEnergy()

        # totalKE3 = total_kinetic_energy(world)
        # totalPE3 = total_potential_energy(world)
        # energyDiff3 = (totalKE3 + totalPE3) - startingTotalE

        # print "--- after energy adjustments ----"
        # print "new totalKE after adjust is", totalKE3
        # print "new totalPE after adjust is", totalPE3
        # print "energy diff from start after adjust is", energyDiff3

pygame.quit()
sys.exit()
