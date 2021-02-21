#!/usr/bin/python

"""
    run_4h.py 4 Hydrogen Atoms -- 4e and 4p experiment.
    2021-02-21 Curt Welch
"""

from typing import List
import physics_sim as ps


def main():
    world: List[ps.Particle] = [
        ps.Proton(0.2 * ps.Angstrom, 0.0, 0.0),
        ps.Proton(0.5 * ps.Angstrom, 0.0, 0.0),
        ps.Proton(1.0 * ps.Angstrom, 0.1 * ps.Angstrom, 0.0),
        ps.Proton(1.2 * ps.Angstrom, 0.1 * ps.Angstrom, 1.0 * ps.Angstrom),
    ]
    e1 = ps.Electron(0.0 * ps.Angstrom, 0.0, 0.0)
    world.append(e1)
    e1.cur_state.v[1] = .011 * ps.CONST_C
    e2 = ps.Electron(0.8 * ps.Angstrom, 0.0, 0.0)
    world.append(e2)
    e2.cur_state.v[1] = .011 * ps.CONST_C
    world.append(
        ps.Electron(0.0 * ps.Angstrom, 0.4 * ps.Angstrom, 0.2 * ps.Angstrom))
    world.append(
        ps.Electron(0.0 * ps.Angstrom, 0.1 * ps.Angstrom, 0.2 * ps.Angstrom))

    sim = ps.Simulation(world, title="4H Test")
    sim.run()


if __name__ == '__main__':
    main()
