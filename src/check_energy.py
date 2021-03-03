#!/usr/bin/python

"""
    check_energy.py
        Does an EP pair with a large radius hve the same
        total E=KE+PE as one with a small radius?
        I know the total system energy will stay the same
        but I don't know for sure if an EP pair has the same
        Energy no mater what configuration it is in.

        Instead of working out the math, I decided to just use
        the code I already have working to Calculate it for me.


    2021-03-02 Curt Welch
"""

import physics_sim as ps


def main():
    last_e = 0.0
    for r in range(1, 10):
        radius = r/10
        # radius /= ps.ANGSTROM
        sim = ps.Simulation()
        sim.add_ep_a((0.0, 0.0, 0.0), radius=radius)
        pe = sim.total_potential_energy()
        ke = sim.total_kinetic_energy()
        diff_e = ke+pe - last_e
        print(f"{radius=} PE:{pe:.3e} KE:{ke:.3e} Total:{pe+ke:.3e} Change from Last:{diff_e:.3e}")
        last_e = ke+pe

    sim = ps.Simulation(dt_max=1e-20)
    sim.add_ep_a((0.0, 0.0, 0.0), radius=1/ps.ANGSTROM)
    # sim.run()


if __name__ == '__main__':
    main()
