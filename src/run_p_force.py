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
from numpy.linalg import norm

Force_Title = ""
Total_force = None


def main():
    global Total_force
    # Total_force = total_force
    Total_force = combined_es_p_force2b
    init_title()

    # do_1h()
    # do_2h()
    # do_2pe()
    do_6h()


def init_title():
    # Create a fake sim and make it do a force calculation
    # which sets Force_Title
    # Stupid hack but it works.  Until it doesn't.
    sim = ps.Simulation(total_force=Total_force)
    sim.add_ep_a((0.0, 0.0, 0.0))
    sim.init_world()


def do_1h():
    """ Simple test of p force with a single EP pair. """
    sim = ps.Simulation(title="do_1h v_force" + Force_Title,
                        total_force=Total_force,
                        dt_max=1e-22,
                        # pixels_per_angstrom=50,
                        pixels_per_angstrom=20000,
                        )

    # Black hole in the middle! :)
    # n = sim.add_n((0, 0, 0), n=2)
    # n.lock_in_place = True

    # One EP in elliptical orbit.
    e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.001)
    e.cur_state.v[1] *= 1.3

    # f = 200
    #
    # # sim.add_p_a((0.3, 0.0, 0.0))
    # # sim.add_ep_a((0.5, 0.0, 0.0))
    # # sim.add_ep_a((0.5, 0.5, 0.5))
    # e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.1)
    # v = e.cur_state.v[1]  # v y of electron
    # e.cur_state.v[1] += v / f
    # p.cur_state.v[1] += v / f

    # e.cur_state.v[1] *= 1.3
    # sim.add_p_a((0.2, 0.0, 0.0))
    # e, p = sim.add_ep_a((0.6, 0.0, 0.1), radius=0.08)

    # # sim.add_ep_a((1.0, 0.0, 0.2), radius=0.05)
    # e, p = sim.add_ep_a((1.0, 0.5, 0.2), radius=0.2)
    # v = e.cur_state.v[1]  # v y of electron
    # e.cur_state.v[1] -= v / f
    # p.cur_state.v[1] -= v / f

    sim.run()


def do_2h():
    """ Simple test of p force with a single EP pair. """
    sim = ps.Simulation(title="do_2h " + Force_Title,
                        total_force=Total_force,
                        dt_max=1e-20,
                        # pixels_per_angstrom=50,
                        # pull_center_force=1e-10,
                        )

    # sim.add_p_a((0.3, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.0, 0.0))
    # sim.add_ep_a((0.5, 0.5, 0.5))
    e, p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.1)
    e.cur_state.v[1] *= 1.2

    # sim.add_p_a((0.2, 0.0, 0.0))
    # e, p = sim.add_ep_a((0.6, 0.0, 0.1), radius=0.08)

    # sim.add_ep_a((1.0, 0.0, 0.2), radius=0.05)
    _e, _p = sim.add_ep_a((1.0, 0.5, 0.2), radius=0.05)

    sim.run()


def do_2pe():
    """ Simple test of 2P and one E. """
    sim = ps.Simulation(title="do_2pe " + Force_Title,
                        pixels_per_angstrom=10000,
                        total_force=total_force,
                        dt_max=5e-23)

    # e: ps.Electron
    _e, _p = sim.add_ep_a((0.0, 0.0, 0.0), radius=0.001)
    # Swap y and z velocity to change orbit direction
    # e.cur_state.r[0], e.cur_state.r[2] = e.cur_state.r[2], e.cur_state.r[0]
    # e.charge *= 1.1

    sim.add_p_a((0.01, 0.0, 0.0001))

    sim.run()


def do_6h():
    """ Test of 6H (6P 6E) """
    sim = ps.Simulation(title="do_6h [combined pullc] " + Force_Title,
                        pixels_per_angstrom=10000,
                        total_force=Total_force,
                        dt_max=1e-20,
                        # dt_max=1e-30,
                        pull_center_force=4e-7,
                        )

    e, p = sim.add_ep_a((0.0, 0.0, 0.000), radius=0.001)
    print(p.cur_state)
    print(e.cur_state)
    e, p = sim.add_ep_a((0.009, 0.001, 0.01), radius=0.001454)
    print(p.cur_state)
    print(e.cur_state)
    _e, _p = sim.add_ep_a((0.02, 0.0, 0.03), radius=0.0008)
    _e, _p = sim.add_ep_a((0.0, 0.02, 0.01), radius=0.002)
    _e, _p = sim.add_ep_a((0.01, 0.02, 0.02), radius=0.0025)
    _e, _p = sim.add_ep_a((0.02, 0.015, 0.03), radius=0.002)

    sim.run()


def total_force(p1_state: ps.ParticleState, p2_state: ps.ParticleState):
    """ Compute 3D force vector on p1, because of p2. """

    if isinstance(p1_state.p, ps.Neutron):
        # these fake Neutrons only have forces with Protons
        # Not with each other, or not with Electrons.
        if not isinstance(p2_state.p, ps.Proton):
            return np.zeros(3)

    if isinstance(p1_state.p, ps.Neutron):
        # these fake Neutrons only have forces with Protons
        # Not with each other, or not with Electrons.
        if not isinstance(p2_state.p, ps.Proton):
            return np.zeros(3)

    # Electrostatic force
    f: ndarray = p1_state.es_force(p2_state)

    if isinstance(p1_state.p, ps.Neutron) or \
            isinstance(p1_state.p, ps.Neutron):
        # No p force for Neutrons
        return f

    # pf = p_force3(p1_state, p2_state)

    # New p force
    # f += pf
    # f = combined_es_p_force(p1_state, p2_state)

    f += p_force2b(p1_state, p2_state)

    return f


def p_force1(p1_state: ps.ParticleState,
             p2_state: ps.ParticleState) -> ndarray:
    """
        Perpendicular force.  Perpendicular to r vector but in the same plane
        with the r and v so that it speeds up or slows down the particles.

        like v_force we want to reduce dv/dr so the V vector is parallel to the
        center line as it is in a circular orbit, but we are doing it in a
        complex way by working with the orbital dynamics to make this work.  We
        speed it up as the electron is approaching the proton which forces it
        into a higher orbit on the near side, and we slow it down as it leaves,
        which lowers the orbit on the high side.

        RESULT:  This didn't work.  Perpendicular to R vector was not energy
        conserving.  Decided Perpendicular to V is safe!  See p_force2() for
        that test.

        This code below is hacked for testing ideas and doesn't implement the
        idea described above (stuff in comments does).

    """
    global Force_Title
    Force_Title = "p_force1"

    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    dr_hat = dr / r
    f_vec = (dr_hat * dv.dot(dr_hat) * ps.CONST_KE *
             p1_state.p.charge *
             p2_state.p.charge / (r * r * ps.CONST_C))
    #
    # # The below is rotating the force vector 90 deg in the
    # # dv dr plane.  Probably a better way to do this.
    # dv_cross_dr = np.cross(dv, dr)
    # len_dv_cross_dr = norm(dv_cross_dr)
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

    # es: ndarray = p1_state.es_force(p2_state)
    # m = np.cross(es, dv) * 0.0000004
    # M is a made up new force that is orthogonal to both dv
    # and es.  Testing if this messes up conservation of energy.
    # print("m is", m)

    # To test conservation of energy theory, make a force perpendicular
    # to v which is NOT perpendicular to r.  My current thinking is
    # leading me to believe that any force perpendicular to V will be
    # energy safe.
    # Ok, just use the m from above, and cross again to flip it 90deg.
    pv = dv.copy()
    pv[0], pv[1] = pv[1], -pv[0]  # swap x and y
    # print("dv and pv", dv, pv)
    pv = pv / norm(pv)  # Turn into unit vector
    # print("pv normalized", pv, "len is", norm(pv))
    pv *= norm(f_vec)  # match magnitude of f_vec

    # print("pv is", pv)
    # print("pv[z]", pv[2])
    # print("pv dot dv should be zero", pv.dot(dv))

    # This is a test to see if a force perpendicular to V is energy
    # conserving even if not perpendicular to r.  Answers seems to
    # be yes. But this random test produced very interesting and
    # odd results.  Turned ellipse into circular orbit very quickly
    # but then started to add tiny loops in the orbit -- all in 2D.
    # Seemed to be if the electron got too far away it would do
    # loop maybe?

    return pv * -1000.0


def p_force2(p1_state: ps.ParticleState,
             p2_state: ps.ParticleState) -> ndarray:
    """
        Perpendicular force.  Try 2. Perpendicular to V vector but in
        the same plane with the r and v so that it TURNS the particle
        without speeding up or slowing down directly. This seems to
        maintain conservation of energy. Try 1 did not.

        BUG: Using dv cross dr as I did in the code below means the
        vectors is strongest at parallel, and weakest at 90 deg but
        then I multiply by v which makes it strong at 90 and weak
        at parallel. But the two combined ends up with something like
        strong at 45 weak everywhere else.  The direction is always
        right but the magnitude is not how I wanted to code it. This
        same bug is in p_force3() and p_force4().  I'll have to fix
        and experiment now. (actually decided it wasn't wrong) [later
        wrote p_force2b() to try the fixed version -- works well]

        So the thinking here is to push the particles into a parallel
        path where dv/dr is zero. When the particles are approaching
        each other we turn them away, and when they are headed away, we
        turn them towards each other.  This should, maybe, force a
        single EP pair to form a circular orbit. When in a circular
        orbit with the particles having zero dr/dt, there will be no
        extra force applied.

        Strength of force is relative to inverse square like the K
        force.  Linearly stronger with speed.  For now, will scale it
        so it is the same strength as K at the speed of light.
        Also scale it relative to the projection of V on R.  So full
        strength if V points in line with R and zero when 90 deg to it.
        aka V dot R.

        Results: works well to turn elliptical orbits into circular but
        does not work well to share energy between atoms. When atoms
        get too close they tend transfer too much energy to one electron
        which tears it away from the proton and flies wild. This
        issue was addressed with p force 3.
    """
    global Force_Title
    Force_Title = "p_force2"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        return np.zeros(3)  # Blow up, just punt
    v = norm(dv)
    if v == 0.0:
        return np.zeros(3)  # no force in this case
    dr_hat = dr / r
    dv_hat = dv / v
    v = (dv.dot(dr) * np.abs(ps.CONST_KE *
                             p1_state.p.charge *
                             p2_state.p.charge / (r * r * r * ps.CONST_C)))
    # Below without abs() it flips sign for EE and PP which
    # makes them move non-parallel.  That was unstable:
    # v = -(dv.dot(dr) * ps.CONST_KE *
    #       p1_state.p.charge *
    #       p2_state.p.charge / (r * r * r * ps.CONST_C))
    # v is positive when going away, negative when approaching
    # same for all mixes of particles.
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat) * v

    return f_vec


def p_force3(p1_state: ps.ParticleState,
             p2_state: ps.ParticleState) -> ndarray:
    """
        Perpendicular force.  Try 3. Keeping force perpendicular to the
        v to conserve energy. PE turn parallel.  PP and EE turn towards
        each other.

        ERROR:  What I coded, and what turned out to work, was PP and EE
        turn AWAY from each other!  I had the sign backward and never
        realized it!  When I tried "turn towards" it failed to balance
        energy between atoms.

        Same as try 2 for PE, but changed PP and EE so they turn towards
        each other.  The idea being when headed together they drive their
        velocity towards zero.  This drives relative V to zero to help
        balance energy across atoms.  This worked.  On a 2H and 6H tests
        the atoms stayed well balanced and never did an electron get
        stripped off with too much energy.

        The math was hard to do this as I intended so I cheated and
        used v dot r which means strong pull at 90 deg. weak at 0 or
        180. The intent was strong at 180 weak at 0 with linear
        difference based on angle.  Could not code that with dot and
        cross.  Would have needed to calculate the angle and work with
        that which I did not bother with.

        Results: worked very well to both turn orbits circular and to
        keep energy balanced.  With 2H test KE of E was 34% 65% which
        shows this interesting tendency to balanced at multiples which
        is something real atoms seem to do.  And with 6H test the KE of
        the 6 electrons ranged from 13.2% to 19.8% (on one point when I
        stopped to sample) which is well balanced around 1/6 or 16.7%
        each.

        Result later: I think this stopped working when I modified the
        the display logic to auto adjust to track a frame rate which
        then changed how often it was checking for bounce which them,
        somehow, broke this. Ugh.  This is also messy with with the
        sign test to use different code for EP vs EE PP.  There must
        be something cleaner.
    """
    global Force_Title
    Force_Title = "p_force3"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        return np.zeros(3)  # Blow up, just punt
    v = norm(dv)
    if v == 0.0:
        return np.zeros(3)  # no force in this case
    dr_hat = dr / r
    dv_hat = dv / v
    sign = np.sign(p1_state.p.charge * p2_state.p.charge)
    if sign > 0:
        # EE or PP make it turn towards each other
        # WOW, I don't see a way to do this with dot and cross
        # products.  OK, I can cheat making it maximal strong a 90 deg
        # and weaker as it moves away from that.
        f_vec = np.cross(dv_hat, np.cross(dr_hat, dv_hat))
        es_mag = (ps.CONST_KE *
                  p1_state.p.charge *
                  p2_state.p.charge / (r * r))
        f_mag = es_mag * v / ps.CONST_C
        f_vec *= f_mag

        # # Debug stuff
        # f_hat = f_vec / norm(f_vec)
        # p = np.cross(dv_hat, dr_hat)
        # # print(f"dr_hat", dr_hat)
        # # print(f"dv_hat", dv_hat)
        # # print(f" f_vec", f_vec)
        # # print(f" f_hat", f_hat)
        # # print(f"v x r   ", p)
        # # print(f"fh dot r", f_hat.dot(dr_hat), "should be positive always")
        # # print(f"f dot vh", abs(round(f_hat.dot(dv_hat), 10)), "should be zero always")
        # # print("fh dot p ", abs(round(f_hat.dot(p), 10)), "should be zero shows f in plane with r v")
        # # print()
        # assert f_hat.dot(dr_hat) > 0, "f dot dr is negative"
        # assert abs(round(f_hat.dot(dv_hat), 10)) == 0.0, "f dot v not zero"
        # assert abs(round(f_hat.dot(p), 10)) == 0.0, "f dot p not zero"

        return f_vec

    # EP: turn towards parallel
    v = (dv.dot(dr) * np.abs(ps.CONST_KE *
                             p1_state.p.charge *
                             p2_state.p.charge / (r * r * r * ps.CONST_C)))
    # v = -(dv.dot(dr) * ps.CONST_KE *
    #       p1_state.p.charge *
    #       p2_state.p.charge / (r * r * r * ps.CONST_C))
    # v is positive when going away, negative when approaching
    # same for all mixes of particles.
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat) * v

    return f_vec


def p_force4(p1_state: ps.ParticleState,
             p2_state: ps.ParticleState) -> ndarray:
    """

        Perpendicular force.  Try 4. Same idea as 3, but for EE and PP we turn
        in line with r.  Instead of turning towards each other, turn away or
        towards, whichever is closer. 90 deg shift of what we do for PE.  Math
        is cleaner for this than Try 3.

        RESULTS: fails to distribute energy from atom to atom.

    """
    global Force_Title
    Force_Title = "p_force4"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        return np.zeros(3)  # Blow up, just punt
    v = norm(dv)
    if v == 0.0:
        return np.zeros(3)  # no force in this case
    dr_hat: ndarray = dr / r
    dv_hat: ndarray = dv / v
    sign = np.sign(p1_state.p.charge * p2_state.p.charge)
    if sign > 0:
        # EE or PP make it turn in-line with r.
        f_vec = np.cross(dv_hat, np.cross(dr_hat, dv_hat))
        f_vec *= np.sign(dv_hat.dot(dr_hat))
        # f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)
        es_mag = (ps.CONST_KE *
                  p1_state.p.charge *
                  p2_state.p.charge / (r * r))
        f_mag = es_mag * v / ps.CONST_C
        f_vec *= f_mag
        return f_vec

    # EP: turn towards parallel
    v = (dv.dot(dr) * np.abs(ps.CONST_KE *
                             p1_state.p.charge *
                             p2_state.p.charge / (r * r * r * ps.CONST_C)))
    # v = -(dv.dot(dr) * ps.CONST_KE *
    #       p1_state.p.charge *
    #       p2_state.p.charge / (r * r * r * ps.CONST_C))
    # v is positive when going away, negative when approaching
    # same for all mixes of particles.
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat) * v

    return f_vec


def p_force5(p1_state: ps.ParticleState,
             p2_state: ps.ParticleState) -> ndarray:
    """
        Perpendicular force.  Try 5.

        Back to same concept as 3, but try a new math.

        New math didn't work.  And when I coded "turn towards"
        it failed to share energy across atoms. This then led
        to me figuring out I coded try 3 backward!  I actually
        coded "turn away", by mistake but that is what worked!

        So this ended up being exactly the same as force 3.

        So, new idea.  what if PE is "Turn towards" instead of turn to
        parallel?  That would make the math simple with no special case
        for PP vs PE.  Will it have a side effect of making circular
        orbits as well?

        Result (turn away EE PP, towards PE): H2 no round orbits formed.
        H6 blew up when PE got too close.  Doesn't seem good but the
        blow up was really a simulation limit not a force error. But
        doesn't feel good.

        What about always turn away?

        Result (all turn away): doesn't turn elliptical into circular.
        Doesn't seem to share energy very well, but Does not blow up
        with electrons flying away so far as I've seen.  So it's doing
        something right.  H6 test which is super slow. Has such weird
        orbits that the DT is jumping up and down and generally running
        very slow.  So hard to tell what it's doing. Atoms seem to push
        away from each other.

    """
    global Force_Title
    Force_Title = "p_force5"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        # Two particles are at the same location in space.
        return np.zeros(3)  # Blows up, just punt
    v = norm(dv)
    if v == 0.0:
        # Zero relative velocity.  This is not a bad thing
        # and can happen at startup but is highly unlikely
        # to happen otherwise.
        return np.zeros(3)  # no force in this case
    dr_hat: ndarray = dr / r
    dv_hat: ndarray = dv / v
    p_sign = np.sign(p1_state.p.charge * p2_state.p.charge)
    if p_sign > 0:
        # EE or PP them turn towards each other.
        f_vec = np.cross(dv_hat, np.cross(dv_hat, dr_hat))
        # f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)
        es_mag = (ps.CONST_KE *
                  p1_state.p.charge *
                  p2_state.p.charge / (r * r))
        f_mag = es_mag * v / ps.CONST_C
        f_vec *= f_mag
        f_vec *= -1  # make them turn away from each other!!!!!
        # f_vec *= p_sign     # turn away for all particle pairs!

        # # Debug stuff
        # f_hat = f_vec / norm(f_vec)
        # p = np.cross(dv_hat, dr_hat)
        # # print(f"dr_hat", dr_hat)
        # # print(f"dv_hat", dv_hat)
        # # print(f" es_mag", es_mag)
        # # print(f" f_vec", f_vec, "norm", norm(f_vec))
        # # print(f" f_hat", f_hat)
        # # print(f"v x r   ", p)
        # # print(f"fh dot r", f_hat.dot(dr_hat), "should be positive always")
        # # print(f"f dot vh", abs(round(f_hat.dot(dv_hat), 10)), "should be zero always")
        # # print("fh dot p ", abs(round(f_hat.dot(p), 10)), "should be zero shows f in plane with r v")
        # # print()
        # # assert f_hat.dot(dr_hat) > 0, "f dot dr is negative"
        # assert abs(round(f_hat.dot(dv_hat), 10)) == 0.0, "f dot v not zero"
        # assert abs(round(f_hat.dot(p), 10)) == 0.0, "f dot p not zero"

        return f_vec

    # EP: turn towards parallel
    v = (dv.dot(dr) * np.abs(ps.CONST_KE *
                             p1_state.p.charge *
                             p2_state.p.charge / (r * r * r * ps.CONST_C)))
    # v = -(dv.dot(dr) * ps.CONST_KE *
    #       p1_state.p.charge *
    #       p2_state.p.charge / (r * r * r * ps.CONST_C))
    # v is positive when going away, negative when approaching
    # same for all mixes of particles.
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat) * v

    return f_vec


def m_force(p1_state: ps.ParticleState,
            p2_state: ps.ParticleState) -> ndarray:
    """
        Magnetic Force.

        2021-3-17 Test new understanding of mag force

        I think I grasp Maxwell's equations now.  And
        though I have confusion about how exactly to apply
        it to two particles, I have an idea worth trying.

        From the old code from a few years back:

        # Calculate the force created on self, by the magnetic
        # field generated by p.
        # B = (1e-7 q1 v1 x rHat) / r²
        # F = q2 V2 x B
        # F = q2 V2 X (1e-7 q1 v1 x rHat) / r²
        # rHat is the unit vector pointing from p to self

        The code above left me confused because if V1 and V2
        were in the same direction there should be a strong
        magnetic field for two parallel wires with currents
        going the same same way.  But Now I understand that
        a change in E is equivalent to current in Maxwell's
        equation and v x r_hat measures that.  And we need
        a double cross product to get the B force pointing
        the right way, and that gives us the need V² in the
        above equations to make B force the same as K force
        when V is c.

        Makes we wonder if the F vector should rotate based
        on the speed?  So instead of K plus the same size
        force at 90 deg, make K rotate based on speed?

        But first, just have B force as:

        F = 1e-7 * q1 * q2 * V X (V x rHat) / r²

        NO, NO, I was thinking V dot rHat.  But the formula needs to be
        cross to get the vector pointing the right way.  If we flip one
        of the signs or something, does the sum of the two cross
        products counter act each other and make it act like dot so
        that it's zero at 90 degrees and greater for others?  I have to
        play with the math tomorrow.  No, the second cross will just be
        a multiply and a angle twist since the first guarantees the
        second will be 90 deg.


    """
    global Force_Title
    Force_Title = "m_force in plane of paper"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        # Two particles are at the same location in space.
        return np.zeros(3)  # Blows up, just punt
    v = norm(dv)
    if v == 0.0:
        # Zero relative velocity.  This is not a bad thing
        # and can happen at startup but is highly unlikely
        # to happen otherwise.
        return np.zeros(3)  # no force in this case
    dr_hat: ndarray = dr / r
    dv_hat: ndarray = dv / v

    #     # # Debug stuff
    #     # f_hat = f_vec / norm(f_vec)
    #     # p = np.cross(dv_hat, dr_hat)
    #     # # print(f"dr_hat", dr_hat)
    #     # # print(f"dv_hat", dv_hat)
    #     # # print(f" es_mag", es_mag)
    #     # # print(f" f_vec", f_vec, "norm", norm(f_vec))
    #     # # print(f" f_hat", f_hat)
    #     # # print(f"v x r   ", p)
    #     # # print(f"fh dot r", f_hat.dot(dr_hat), "should be positive always")
    #     # # print(f"f dot vh", abs(round(f_hat.dot(dv_hat), 10)), "should be zero always")
    #     # # print("fh dot p ", abs(round(f_hat.dot(p), 10)), "should be zero shows f in plane with r v")
    #     # # print()
    #     # # assert f_hat.dot(dr_hat) > 0, "f dot dr is negative"
    #     # assert abs(round(f_hat.dot(dv_hat), 10)) == 0.0, "f dot v not zero"
    #     # assert abs(round(f_hat.dot(p), 10)) == 0.0, "f dot p not zero"
    #
    #     return f_vec

    # # Old p_force2():
    # # v = -(dv.dot(dr) * ps.CONST_KE *
    # #       p1_state.p.charge *
    # #       p2_state.p.charge / (r * r * r * ps.CONST_C))
    # # v is positive when going away, negative when approaching
    # # same for all mixes of particles.
    # # f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat) * v
    #
    # This is the actual old stuff with the abs():
    # # # EP: turn towards parallel
    # # v = (dv.dot(dr) * np.abs(ps.CONST_KE *
    # #                          p1_state.p.charge *
    # #                          p2_state.p.charge / (r * r * r * ps.CONST_C)))

    # 3-17 YEAH THIS OLD code was close to correct!  But I need it to
    # be full strength with V is up or down! This didn't do that!
    # The new code above isn't right.  IT needs to be rotated 90 deg
    # so that's it's zero Magnetic when V is 90 deg to right.
    # Will figure this out tomorrow.

    # 3-18 Right, so the old code I first created in p_force2() is close
    # what I want now, but without the going to zero at effect at 180.
    # So I'm just going to fix that by normalizing the cross products
    # and calculating the correct magnitude separately as I was doing
    # in v_force2()

    # This is the normal out of paper vector:
    # f_vec = np.cross(dr_hat, dv_hat)

    # This is on the paper vector:
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)

    f_vec_hat = f_vec / norm(f_vec)

    # f magnitude = v·r_hat 1e-7 qqv²/r²
    m_mag = (1e-7 * dv.dot(dr_hat) * v *
             p1_state.p.charge * p2_state.p.charge / (r * r))
    m_f = f_vec_hat * (m_mag * -1)

    # fdv = m_f.dot(dv)
    # global fdv_max
    # fdv_max = max(np.abs(fdv), fdv_max)
    # if fdv > 1e-17:
    #     print("f dot v should be 0 is:", fdv, "max:", fdv_max)

    return m_f


def v2_force2(p1_state: ps.ParticleState,
              p2_state: ps.ParticleState) -> ndarray:
    """
        Speed limit Force.

        2021-3-21 Test new idea.

        Same idea as older v_force().  We apply force
        to reduce speed towards zero.

        But v_force only worked in the Y dimension (in line with
        dr).  This will test in line with V.  And this will use
        v^2 which v_force() did not.

        ResultL total disaster.  Causing V to drop to zero just causes
        EP pairs to collapse.  Stupid.
    """
    global Force_Title
    Force_Title = "v2_force2 [v² force in line with V] "
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        # Two particles are at the same location in space.
        return np.zeros(3)  # Blows up, just punt

    f_vec = dv * np.abs(dv) * (-1e-7 * np.abs(p1_state.p.charge * p2_state.p.charge) / (r * r))
    return f_vec


def p_force2b(p1_state: ps.ParticleState,
              p2_state: ps.ParticleState) -> ndarray:
    """
        2021-3-21 force 2 again, but fix the math.

        Turn towards parallel to make relative V drop to 0.
        Do the same for all particle combinations.
        The larger dr/dt, the more the force.  This is what was
        coded in force2, but there I made the turning force drop
        to zero at 0 and 180 by accident of code.
    """
    global Force_Title
    Force_Title = "p_force2b [turn towards parallel for all] "
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        # Two particles are at the same location in space.
        return np.zeros(3)  # Blows up, just punt
    v = norm(dv)
    if v == 0.0:
        # Zero relative velocity.
        return np.zeros(3)  # no force in this case
    dr_hat = dr / r
    dv_hat = dv / v
    # f_vec always points down
    f_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)
    f_vec /= norm(f_vec)
    # Scaler force:
    dr_dt = dv.dot(dr_hat)   # positive goes away
    f = dr_dt * np.abs(dr_dt * 1e-7 * p1_state.p.charge * p2_state.p.charge / (r * r))
    return f_vec * f


def combined_es_p_force2b(p1_state: ps.ParticleState,
                          p2_state: ps.ParticleState) -> ndarray:
    """
        Combined P-fore2b and es force.

        Runs a little aster if we calculate them both at the same
        time.  Only about 10% faster. Not a huge deal.
    """
    global Force_Title
    Force_Title = "combined p_force2b and es"
    dv: ndarray = p1_state.v - p2_state.v
    dr: ndarray = p1_state.r - p2_state.r
    r = norm(dr)
    if r == 0.0:
        return np.zeros(3)  # blow up - abort
    dr_hat = dr / r
    qq_rr = p1_state.p.charge * p2_state.p.charge / (r * r)
    # es force is:
    es_vec = dr_hat * ps.CONST_KE * qq_rr

    v = norm(dv)
    if v == 0.0:
        # Zero relative velocity.
        return es_vec  # no v force in this case

    dv_hat = dv / v
    # vf_vec always points down
    vf_vec = np.cross(np.cross(dr_hat, dv_hat), dv_hat)
    vf_hat = vf_vec / norm(vf_vec)
    # v force Scaler magnitude:
    dr_dt = dv.dot(dr_hat)   # positive goes away
    vf = dr_dt * np.abs(dr_dt * 1e-7 * qq_rr)

    return es_vec + vf_hat * vf


if __name__ == '__main__':
    main()
