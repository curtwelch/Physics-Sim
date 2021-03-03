#!/usr/bin/python3

"""
    run_pn.py one Proton one Fake Neutron
    2021-03-03 Curt Welch
"""

import physics_sim as ps


def main():
    """ Simple test of a Proton and a fake Neutron. """
    sim = ps.Simulation(title="PN Test nc 1.1", dt_max=2e-20)

    v = .010 * ps.CONST_C   # C is the speed of light.

    sim.add_p_a((0.0, 0.0, 0.0), v=(0.0, -v/100, 0.0))
    sim.add_n_a((0.2, 0.0, 0.0), v=(0.0, 0.000, 0.0), nc=1.1)
    sim.add_p_a((0.4, 0.0, 0.0), v=(0.0, v/100, 0.0))

    # sim.add_e_a((0.8, 0.0, 0.0), v=(0.0, v, 0.0))
    # sim.add_e_a((-0.4, 0.0, 0.0), v=(0.0, -v, 0.0))

    sim.run()


if __name__ == '__main__':
    main()
