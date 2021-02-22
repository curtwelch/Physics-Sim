#!/usr/bin/python

"""
    run_4h.py 4 Hydrogen Atoms -- 4e and 4p experiment.
    2021-02-21 Curt Welch
"""

import physics_sim as ps


def main():
    sim = ps.Simulation(title="4H Test", dt_max=2e-20)

    v = .012 * ps.CONST_C

    sim.add_p_a((0.2, 0.0, 0.0))
    sim.add_e_a((0.2, -0.1, 0.0), v=(v, 0.0, 0.0))

    sim.add_p_a((0.5, 0.0, 0.0))
    sim.add_e_a((0.5, 0.1, 0.1), v=(-v*.8, 0.0, 0.0))

    sim.add_p_a((1.0, 0.1, 0.0))
    sim.add_e_a((1.0, 0.2, 0.0), v=(v*1.1, 0.0, 0.0))

    sim.add_p_a((0.5, 0.5, 1.0))
    sim.add_e_a((0.6, 0.6, 1.3), v=(0.0, v/2, v/2))

    sim.run()


if __name__ == '__main__':
    main()
