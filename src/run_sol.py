#!/usr/bin/python

"""
    Speed of light experiment.

    2021-02-11 Curt Welch curt@kcwc.com

    Create a 1D string of electrons and make them stable by adding and extra
    artificial force to hold it in place. Then wiggle the first one, and watch
    how the rest respond. I'm hoping it will show a wave moving through the
    string even though the simulation uses no time delay for the Coulomb force
    moving them.

    First calculation.
    The 10 electrons, passed the energy of the first e to the 9th
    in this much time.  This is where the KE of the 9th maxed.
    spacing is 2/10 A so 9 times that is
    18/10 A in this much time:
    Time now is 0.0000000000000000243099999999994800664042

    Output is:
    Speed: speed=7404360.349, c=299792458.000 speed/c=0.025

    Conclusion:

    No where near C, but yet does show how KE can move through the particles at
    a speed below the infinite speed of force interactions.

    2021-02-21 Split this out as run_sol.py from the old code and
    rewrote to get working with the new code.

    Curt Welch curt@kcwc.com
"""

import physics_sim as ps


def main():
    """ Speed of Light experiment, """

    c = ps.CONST_C  # Speed of Light in m/s

    sim = ps.Simulation(title="SOL Test")

    num = 20
    spacing = (2 / num)
    for cnt in range(num):
        sim.add_e_a((cnt * spacing, 0.0, 0.0))

    # world[0].p.lock_in_place = True
    sim.world[4].cur_state.v[0] = .05*c

    sim.init_world()
    # The particle states, including the forces have been re-defined now.

    for p in sim.world:
        # Set static starting force to keep particles in place.
        # Trying to fake an infinitely long line of particles.
        p.static_f[:] = -p.cur_state.f

    sim.run()


if __name__ == '__main__':
    main()
