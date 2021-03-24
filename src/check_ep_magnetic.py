#!/usr/bin/python

"""
    2021-03-23 check_ep_magnetic.ph

    New line of thinking.  See if I can make formals for a magnetic field
    work of we assume the source of the magnetic field is an EP pair.

    Standard math with currents and changing E fields don't seem to explain
    magnetism forces when we reduce down to this level of individual particles.

    But if we make the assumption that it takes an EP pair in orbit where
    the sum of the effect of the two particles on an external particle explains
    the magnetic force, maybe we can them make the math work!

    With some quick diagrams and logic, it looks like this might be sort
    of possible. But I now need to write the code, and test what the combined
    force vectors look like!  Can we make the combined force vectors from
    a C and from an P work like Magnetism with the idea that he EP pair is
    generating a magnetic field!

    What the back of the envelop thinking yield was the idea that the strength
    of the field needed to be relative the MOMENTUM of the EP particle.  Not
    the V.   The P needs to be 1800 times stronger effect for the same V.

    And strength of mag F needs to be relative to norm(dV).  And direction
    is 90 deg to motion of dV.

    Oh, and direction flips when sign flips.

    So similar to normal mag force formulas but using momentum instead of
    current is the big difference.

"""

import physics_sim as ps
import numpy as np
from numpy import ndarray
from numpy.linalg import norm


def main():
    print("Mag test")
    sim = ps.Simulation()
    e, p = sim.add_ep_a((0.0, 0.0, 0.0))
    e2 = sim.add_e_a((10000.0, 0.0, 0.0), v=(0.0, 0.0, 0.0))
    print(e.cur_state)
    print(p.cur_state)
    vp = np.abs(p.cur_state.v[1])
    print("======= v is zero:")
    do(e, p, e2)

    print("======= Match Proton speed:")
    e2.cur_state.v[1] = -vp
    do(e, p, e2)

    print("======= Match -Proton speed:")
    e2.cur_state.v[1] = vp
    do(e, p, e2)

    print("======= Match Electron speed:")
    e2.cur_state.v[1] = e.cur_state.v[1]
    do(e, p, e2)

    print("======= Match -Electron speed:")
    e2.cur_state.v[1] = -e.cur_state.v[1]
    do(e, p, e2)


def do(e: ps.Particle, p: ps.Particle, e2: ps.Particle):
    print("e2 v is:", e2.cur_state.v)
    fe = force(e, e2)
    fp = force(p, e2)
    print("fe: ", fe)
    print("fp: ", fp)
    print("sum:", fe+fp)


def force(p1: ps.Particle, p2: ps.Particle):
    # print("--- Compute force ----")
    f = force_state(p1.cur_state, p2.cur_state)
    # print("Force is", f)
    # print("---")
    return f


def force_state(p1_state: ps.ParticleState,
                p2_state: ps.ParticleState) -> ndarray:
    dr: ndarray = p1_state.r - p2_state.r
    dv: ndarray = p1_state.v - p2_state.v
    r = norm(dr)
    v = norm(dv)
    if r == 0.0 or v == 0.0:
        return np.zeros(3)
    dr_hat = dr / r
    dv_hat = dv / v
    # print("V is", v, "R is", r)
    # print("dr_hat:", dr_hat)
    # print("dv_hat", dv_hat)
    # vf_vec always points down
    vf_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)
    vf_hat = vf_vec / norm(vf_vec)
    # print("vf_hat", vf_hat)
    # v force Scaler magnitude:
    qq_rr = p1_state.p.charge * p2_state.p.charge / (r * r)
    vf = v * p1_state.p.mass * 1e-7 * qq_rr

    return vf_hat * vf


if __name__ == '__main__':
    main()
