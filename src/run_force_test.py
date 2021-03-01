#!/usr/bin/python

"""
    run_force_test.py
        Test out a new idea for the v_force formula.

    2021-02-23 Curt Welch
"""

import physics_sim as ps
import numpy as np
from numpy import ndarray


def main():
    test_2pe()


def test_6h():
    """ Simple test of add_ep_a(). """
    sim = ps.Simulation(title="New V Force Test", total_force=combined_es_v_force, dt_max=1e-19)

    sim.add_p_a((0.3, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.5, 0.5))
    # sim.add_ep_a((0.0, 0.5, 0.2), radius=0.05)
    # sim.add_ep_a((1.0, 0.0, 0.2), radius=0.05)
    # sim.add_ep_a((1.0, 0.5, 0.2), radius=0.05)

    sim.run()


def test_2pe():
    """ Simple test of 2P and one E. """
    sim = ps.Simulation(title="2PE Test with new v_force",
                        pixels_per_angstrom=10000,
                        total_force=combined_es_v_force,
                        dt_max=5e-23)

    e: ps.Electron
    e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.001)
    # Swap y and z velocity to change orbit direction
    # e.cur_state.r[0], e.cur_state.r[2] = e.cur_state.r[2], e.cur_state.r[0]
    e.charge *= 1.1

    sim.add_p_a((0.01, 0.0, 0.0001))

    sim.run()


def total_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState):
    """ Compute 3D force vector on p1, because of p2. """

    # Electrostatic force
    # f: ndarray = p1_state.es_force(p2_state)

    # f += p1_state.v_force(p2_state)
    f = new_v_force(p1_state, p2_state)

    return f


def new_v_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState) -> ndarray:
    """
        This version of Velocity force acts to reduce ES force
        when particles move towards each other and increase
        as they move away.
    """
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = np.linalg.norm(dr)
    dr_hat = dr / r
    # f_vec = (dr_hat * dv.dot(dr_hat) * ps.CONST_KE *
    #          p1_state.p.charge *
    #          p2_state.p.charge / (r * r * ps.CONST_C))

    # es force is:
    # return dr_unit * CONST_KE * self.p.charge * p_state.p.charge / (r * r)
    es_vec = dr_hat * ps.CONST_KE * p1_state.p.charge * p2_state.p.charge / (r*r)
    new_f = es_vec * dv.dot(dr_hat) / ps.CONST_C
    # print("old combined is:", all_force)
    # print("new          is:", new_f)

    return new_f


def combined_es_v_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState) -> ndarray:
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = np.linalg.norm(dr)
    dr_hat = dr / r
    # f_vec = (dr_hat * dv.dot(dr_hat) * ps.CONST_KE *
    #          p1_state.p.charge *
    #          p2_state.p.charge / (r * r * ps.CONST_C))

    # es force is:
    # return dr_unit * CONST_KE * self.p.charge * p_state.p.charge / (r * r)
    es_vec = dr_hat * ps.CONST_KE * p1_state.p.charge * p2_state.p.charge / (r*r)
    # all_force = es_vec + f_vec

    # Lets simplify
    new_f = es_vec * (1.0 + dv.dot(dr_hat) / ps.CONST_C)
    # print("old combined is:", all_force)
    # print("new          is:", new_f)

    return new_f


if __name__ == '__main__':
    main()
