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

    Output is:

    radius=0.1 PE:-2.306e-17 KE:1.153e-17 Total:-1.153e-17 Change from Last:-1.153e-17
    -pe/ke=2.00000000000000
    radius=0.2 PE:-1.153e-17 KE:5.765e-18 Total:-5.765e-18 Change from Last:5.765e-18
    -pe/ke=2.00000000000000
    radius=0.3 PE:-7.686e-18 KE:3.843e-18 Total:-3.843e-18 Change from Last:1.922e-18
    -pe/ke=2.00000000000000
    radius=0.4 PE:-5.765e-18 KE:2.882e-18 Total:-2.882e-18 Change from Last:9.608e-19
    -pe/ke=2.00000000000000
    radius=0.5 PE:-4.612e-18 KE:2.306e-18 Total:-2.306e-18 Change from Last:5.765e-19
    -pe/ke=2.00000000000000
    radius=0.6 PE:-3.843e-18 KE:1.922e-18 Total:-1.922e-18 Change from Last:3.843e-19
    -pe/ke=2.00000000000000
    radius=0.7 PE:-3.294e-18 KE:1.647e-18 Total:-1.647e-18 Change from Last:2.745e-19
    -pe/ke=2.00000000000000
    radius=0.8 PE:-2.882e-18 KE:1.441e-18 Total:-1.441e-18 Change from Last:2.059e-19
    -pe/ke=2.00000000000000
    radius=0.9 PE:-2.562e-18 KE:1.281e-18 Total:-1.281e-18 Change from Last:1.601e-19
    -pe/ke=2.00000000000000

    Which says larger circular orbits have more total energy.
    And look at that, pe is always twice the ke!

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
        print(f"{-pe/ke=:.14f}")
        last_e = ke+pe

    sim = ps.Simulation(dt_max=1e-20)
    sim.add_ep_a((0.0, 0.0, 0.0), radius=1/ps.ANGSTROM)
    # sim.run()


if __name__ == '__main__':
    main()
