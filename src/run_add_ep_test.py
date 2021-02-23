#!/usr/bin/python

"""
    run_add_ep_test.py 1 Hydrogen Atoms -- one ep pair.
    2021-02-22 Curt Welch
"""

import physics_sim as ps


def main():
    """ Simple test of add_ep_a(). """
    sim = ps.Simulation(title="add_ep_a() Test", dt_max=2e-20)

    sim.add_ep_a((0.0, 0.0, 0.0))

    sim.run()


if __name__ == '__main__':
    main()
