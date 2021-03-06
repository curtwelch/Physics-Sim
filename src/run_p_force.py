#!/usr/bin/python

"""
    run_p_force.py

    More new ideas for force formula.  Calling it p force because
    it will be perpendicular to the K force vector.

    Always apply force perpendicular to center line so we don't
    mess up the PE of the system created by the K force.
    But unlike magnetism concepts where force is applied at 90 deg
    to both K force and to V, we are going to apply the force in
    line with the VR plane so it changes the velocity and transfers
    KE from one particle to the other.

    But there are many ways to do the above. Not sure what might be
    useful or "right".

    Goals are:

    (1) see if this type of perpendicular force vector fixes the loss
    of energy issue that happens in v_force().

    (2) configure it so it turns elliptical orbits into circular oro its
    for a single EP pair.

    (3) Do the above, AND cause multiple EP pairs to synchronize orbits!
    I do not know if this is possible but this is what I'm trying to
    produce with these ideas.

    2021-03-05 Curt Welch
"""

import physics_sim as ps
import numpy as np
from numpy import ndarray


def main():
    do_1h()


def do_1h():
    """ Simple test of p force with a single EP pair. """
    sim = ps.Simulation(title="P Force 1H Test",
                        total_force=total_force,
                        dt_max=1e-20)

    # sim.add_p_a((0.3, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.5, 0.5))
    e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.1)
    # e.cur_state.v[1] *= 1.3
    sim.add_p_a((0.2, 0.0, 0.0))

    # sim.add_ep_a((1.0, 0.0, 0.2), radius=0.05)
    # sim.add_ep_a((1.0, 0.5, 0.2), radius=0.05)

    sim.run()


# def test_2pe():
#     """ Simple test of 2P and one E. """
#     sim = ps.Simulation(title="2PE Test with new v_force",
#                         pixels_per_angstrom=10000,
#                         total_force=combined_es_v_force,
#                         dt_max=5e-23)
#
#     e: ps.Electron
#     e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.001)
#     # Swap y and z velocity to change orbit direction
#     # e.cur_state.r[0], e.cur_state.r[2] = e.cur_state.r[2], e.cur_state.r[0]
#     e.charge *= 1.1
#
#     sim.add_p_a((0.01, 0.0, 0.0001))
#
#     sim.run()


def total_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState):
    """ Compute 3D force vector on p1, because of p2. """

    # Electrostatic force
    f: ndarray = p1_state.es_force(p2_state)

    # New p force
    f += p_force(p1_state, p2_state)
    # f = combined_es_p_force(p1_state, p2_state)

    return f


def p_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState) -> ndarray:
    """
        Perpendicular force.  Perpendicular to r vector but in the same plane
        with the r and v so that it speeds up or slows down the particles.

        like v_force we want to reduce dv/dr so the V vector is parallel to the
        center line as it is in a circular orbit, but we are doing it in a
        complex way by working with the orbital dynamics to make this work.  We
        speed it up as the electron is approaching the proton which forces it
        into a higher orbit on the near side, and we slow it down as it leaves,
        which lowers the orbit on the high side.

    """
    dv: ndarray = p1_state.v - p2_state.v
    # dr: ndarray = p1_state.r - p2_state.r
    # r = np.linalg.norm(dr)
    # dr_hat = dr / r
    # f_vec = (dr_hat * dv.dot(dr_hat) * ps.CONST_KE *
    #          p1_state.p.charge *
    #          p2_state.p.charge / (r * r * ps.CONST_C))
    #
    # # The below is rotating the force vector 90 deg in the
    # # dv dr plane.  Probably a better way to do this.
    # dv_cross_dr = np.cross(dv, dr)
    # len_dv_cross_dr = np.linalg.norm(dv_cross_dr)
    # if len_dv_cross_dr == 0.0:
    #     return np.zeros(3)  # blows up, just punt and return zero
    # dv_cross_dr_hat = dv_cross_dr / len_dv_cross_dr
    # answer = np.cross(f_vec, dv_cross_dr_hat)
    # print()
    # print("dr       ", dr, "r", r)
    # print("dr hat   ", dr_hat)
    # print("dv       ", dv)
    # print("f_vec    ", f_vec)
    # print("vxr had  ", dv_cross_dr_hat)
    # print("answer   ", answer)
    # print("answer dot f_vec should be zero", answer.dot(f_vec))
    # print("answer dot dr    should be zero", answer.dot(dr))
    # print("dv dot dr_hat", dv.dot(dr_hat))

    # # This proved conservation of energy wasn't working when
    # # I kept the force orthogonal to dr.
    # if dv.dot(dr_hat) < 0:
    #     return np.zeros(3)

    es: ndarray = p1_state.es_force(p2_state)
    m = np.cross(es, dv) * 0.0000004
    # M is a made up new force that is orthogonal to both dv
    # and es.  Testing if this messes up conservation of energy.
    # print("m is", m)

    return m


def combined_es_p_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState) -> ndarray:
    """
        not written yet.
    """
    # dv: ndarray = p1_state.v - p2_state.v
    # dr: ndarray = p1_state.r - p2_state.r
    # r = np.linalg.norm(dr)
    # dr_hat = dr / r
    # # f_vec = (dr_hat * dv.dot(dr_hat) * ps.CONST_KE *
    # #          p1_state.p.charge *
    # #          p2_state.p.charge / (r * r * ps.CONST_C))
    #
    # # es force is:
    # # return dr_unit * CONST_KE * self.p.charge * p_state.p.charge / (r * r)
    # es_vec = dr_hat * ps.CONST_KE * p1_state.p.charge * p2_state.p.charge / (r*r)
    # # all_force = es_vec + f_vec
    #
    # # Lets simplify
    # new_f = es_vec * (1.0 + dv.dot(dr_hat) / ps.CONST_C)
    # # print("old combined is:", all_force)
    # # print("new          is:", new_f)
    #
    # return new_f


if __name__ == '__main__':
    main()
