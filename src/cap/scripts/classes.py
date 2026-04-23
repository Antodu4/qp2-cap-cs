#!/usr/bin/env python3
import sys


class Energy:
    """Stores the CAP energy for a single state at a single eta value.

    Attributes:
        eta   : CAP strength parameter value.
        re    : Real part of the CAP energy.
        im    : Imaginary part of the CAP energy.
        c_re  : Real part of the first-order energy correction.
        c_im  : Imaginary part of the first-order energy correction.
    """

    def __init__(self, eta, E_re, E_im, C_re, C_im):
        if not isinstance(eta, float) or not isinstance(E_re, float) \
           or not isinstance(E_im, float) or not isinstance(C_re, float) \
           or not isinstance(C_im, float):
            print("Invalid type when creating Energy class object.")
            sys.exit()

        self.eta = eta
        self.re = E_re
        self.im = E_im
        self.c_re = C_re
        self.c_im = C_im


class Energies:
    """Collection of Energy objects covering a range of eta values for one state.

    Attributes:
        n_eta    : Number of eta values.
        energies : List of Energy objects, one per eta value.
    """

    def __init__(self, n_eta, l_energy):
        if not isinstance(l_energy, list):
            print("Invalid type when creating Energies class object.")
            sys.exit()

        for e in l_energy:
            if not isinstance(e, Energy):
                print("Invalid type, Energies must take a list of Energy objects.")
                sys.exit()

        if len(l_energy) != n_eta:
            print("Wrong number of eta values")
            sys.exit()

        self.n_eta = n_eta
        self.energies = l_energy


class State:
    """Stores data for one electronic state at one CIPSI iteration.

    Attributes:
        n_det    : Number of determinants at this iteration.
        E        : Real CI energy (Hartree).
        pt2      : Stochastic PT2 correction.
        energies : Energies object with CAP energies over the eta range.
        state    : State index (1-based).
    """

    def __init__(self, n_det, E, pt2, energies, state):
        if not isinstance(n_det, int):
            print("Invalid type for n_det when creating State class object.")
            sys.exit()
        if not isinstance(E, float):
            print("Invalid type for E when creating State class object.")
            sys.exit()
        if not isinstance(pt2, float):
            print("Invalid type for pt2 when creating State class object.")
            sys.exit()
        if not isinstance(state, int):
            print("Invalid type for state when creating State class object.")
            sys.exit()

        self.n_det = n_det
        self.E = E
        self.pt2 = pt2
        self.energies = energies
        self.state = state


class Step:
    """Represents one CIPSI iteration, containing data for all electronic states.

    Attributes:
        n_states : Number of electronic states.
        states   : List of State objects, one per state.
    """

    def __init__(self, n_states, states):
        if not isinstance(states, list):
            print("Invalid type when creating Step class object.")
            sys.exit()

        if len(states) != n_states:
            print("Wrong number of states")
            sys.exit()

        self.n_states = n_states
        self.states = states


class Cipsi:
    """Stores all CIPSI iterations for a single calculation.

    Attributes:
        steps : List of Step objects, one per CIPSI iteration.
    """

    def __init__(self, steps):
        if not isinstance(steps, list):
            print("Invalid type, Cipsi class objects are creating with a list of Step.")
            sys.exit()

        self.steps = steps

    def append(self, step):
        """Append a new CIPSI iteration step.

        Args:
            step : Step object to append.
        """
        if not isinstance(step, Step):
            print("Invalid type, append methods in Cipsi must take an Step class object.")
            sys.exit()

        self.steps.append(step)
