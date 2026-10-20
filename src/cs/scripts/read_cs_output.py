#!/usr/bin/env python3
from classes import *


def get_lines_match(f_path, match):
    """Return all lines from a file that contain the given substring.

    Args:
        f_path : Path to the file to search.
        match  : Substring to search for.

    Returns:
        List of stripped lines containing match.
    """
    acc = []
    with open(f_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if match in line:
                acc.append(line.strip())

    return acc


def go_to_match(buff, i, match):
    """Advance index i in buff until a line containing match is found.

    Exits the program if the end of buff is reached without finding match.

    Args:
        buff  : List of strings (file lines).
        i     : Starting index.
        match : Substring to search for.

    Returns:
        Index of the first line at or after i that contains match.
    """
    while i < len(buff):
        if match not in buff[i]:
            i += 1
            continue
        else:
            break

    if i >= len(buff):
        print("Reaches end of file.")
        sys.exit()

    return i


def try_go_to_match(buff, i, match):
    """Advance index i in buff until a line containing match is found.

    Unlike go_to_match, does not abort if the end of buff is reached.

    Args:
        buff  : List of strings (file lines).
        i     : Starting index.
        match : Substring to search for.

    Returns:
        Index of the matching line, or len(buff) if not found.
    """
    while i < len(buff):
        if match not in buff[i]:
            i += 1
            continue
        else:
            break

    return i


# Columns of the TRAJ records printed by cs_print_traj (cs_traj_out.irp.f)
TRAJ_CS_COLUMNS = ['iter', 'N_det', 'theta', 'state',
                   'Re(E_cs)', 'Im(E_cs)', 'Re(E_cs+PT2)', 'Im(E_cs+PT2)',
                   'Re(PT2)', 'Im(PT2)', 'variance']


def is_traj_cs_output(f_path):
    """Return True if the file contains TRAJ records (cs_fci output since the
    unified CS trajectory output), False for the legacy format.

    Args:
        f_path : Path to the file to inspect.
    """
    with open(f_path, 'r') as f:
        for line in f:
            s = line.split()
            if s and s[0] == "TRAJ":
                return True
    return False


def read_traj_cs(f_path):
    """Read the TRAJ records of a cs_fci output file.

    Args:
        f_path : Path to the output file.

    Returns:
        algo : Algorithm name from the CS SCAN SUMMARY banner
               ('fixed_ndet' or 'per_theta'), or 'unknown'.
        df   : pandas DataFrame with the TRAJ_CS_COLUMNS columns, one row
               per (state, theta, iteration).
    """
    import pandas as pd

    algo = 'unknown'
    rows = []
    with open(f_path, 'r') as f:
        for line in f:
            if "CS SCAN SUMMARY" in line and "algorithm =" in line:
                algo = line.split("algorithm =")[-1].strip()
                continue
            s = line.split()
            if len(s) == len(TRAJ_CS_COLUMNS) + 1 and s[0] == "TRAJ":
                rows.append([int(s[1]), int(s[2]), float(s[3]), int(s[4])] +
                            [float(x) for x in s[5:]])

    df = pd.DataFrame(rows, columns=TRAJ_CS_COLUMNS)
    return algo, df


def print_traj_cs_summary(f_path):
    """Print the full trajectory and the converged values per theta of a
    cs_fci output file in the TRAJ format.

    The converged value for each (state, theta) is the record with the largest
    iteration index: the last CIPSI iteration in fixed_ndet mode, the last
    iteration of the CIPSI of that theta in per_theta mode.

    In fixed_ndet mode the PT2 is only computed at the base theta (first theta
    of the grid); for the other theta values PT2 = 0 and E+PT2 = E.

    Final table columns: theta, then for each state
    Re(E), Im(E), Re(E+PT2), Im(E+PT2).

    Args:
        f_path : Path to the output file.
    """
    import sys
    import numpy as np

    algo, df = read_traj_cs(f_path)

    print("# file: {}   algorithm = {}".format(f_path, algo))
    if df.empty:
        print("# No TRAJ record found.")
        return

    print("# Full trajectory (one row per state, theta, iteration)")
    df.to_csv(sys.stdout, index=False)

    last = df.loc[df.groupby(['state', 'theta'])['iter'].idxmax()]
    states = sorted(last['state'].unique().tolist())
    val_theta = sorted(last['theta'].unique().tolist())

    print("# Converged values per theta: theta, then for each state "
          "Re(E), Im(E), Re(E+PT2), Im(E+PT2)")
    cols = ['Re(E_cs)', 'Im(E_cs)', 'Re(E_cs+PT2)', 'Im(E_cs+PT2)']
    data = np.full((len(val_theta), 1 + 4*len(states)), np.nan)
    data[:, 0] = val_theta
    for k, theta in enumerate(val_theta):
        for i, st in enumerate(states):
            row = last[(last['state'] == st) & (last['theta'] == theta)]
            if not row.empty:
                data[k, 1+i*4:1+(i+1)*4] = row[cols].to_numpy()[0]

    for row in data:
        print(', '.join(map(str, row)))
    print()


if __name__ == "__main__":
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    files = sys.argv[1:]

    # New format (TRAJ records printed by cs_print_traj)
    traj_files = [fname for fname in files if is_traj_cs_output(fname)]
    for fname in traj_files:
        print_traj_cs_summary(fname)

    # Legacy format (per-iteration "CS energies" tables printed by cs())
    files = [fname for fname in files if fname not in traj_files]
    if not files:
        sys.exit()

    l_cipsi = []
    for fname in files:
        cipsi = Cipsi([])
        with open(fname, 'r') as f:
            buff = f.readlines()
            i = 0

            n_states = int(get_lines_match(fname, "Number of states:")[0].split()[-1])
            n_theta = int(get_lines_match(fname, "Number of theta values:")[0].split()[-1])

            first = True
            while i < len(buff):

                # Jump to the cs energies
                i = go_to_match(buff, i, "Number of theta values:")
                i += 8
                buff_theta = buff[i:i+n_theta]
                a = []
                for line in buff_theta:
                    for val in line.split():
                        a.append(val)
                a = np.asarray(a, dtype=float, order='F')
                E_theta = a.reshape((n_theta, 2*n_states+1))

                # Jump to the energy corrections
                i = go_to_match(buff, i, "#### Energy corrections ###")
                i += 8
                buff_corr = buff[i:i+n_theta]
                a = []
                for line in buff_corr:
                    for val in line.split():
                        a.append(val)
                a = np.asarray(a, dtype=float, order='F')
                corr = a.reshape((n_theta, 2*n_states+1))

                # N_det
                i = try_go_to_match(buff, i, "Summary at N_det =")
                if i < len(buff):
                    n_det = int(buff[i].split()[-1])
                else:
                    n_det = 0

                # E
                i = try_go_to_match(buff, i, "# E ")
                if i < len(buff):
                    l_str_E = buff[i].split()[2:n_states+2]
                    l_E = [float(e) for e in l_str_E]
                else:
                    l_E = [0.0 for j in range(n_states)]

                # PT2
                i += 3
                if i < len(buff):
                    l_str_pt2 = buff[i].split()[2:2*n_states+2]
                    l_pt2 = [float(e) for e in l_str_pt2[::2]]
                else:
                    l_pt2 = [0.0 for j in range(n_states)]

                l_states = []
                l_theta = E_theta[:, 0]
                for j, E, pt2 in zip(range(1, 2*n_states+1, 2), l_E, l_pt2):
                    l_energy = []
                    l_re = E_theta[:, j]
                    l_im = E_theta[:, j+1]
                    l_corr_re = corr[:, j]
                    l_corr_im = corr[:, j+1]
                    for theta, re, im, c_re, c_im in zip(l_theta, l_re, l_im, l_corr_re, l_corr_im):
                        l_energy.append(Energy(float(theta), float(re), float(im), float(c_re), float(c_im)))

                    energies = Energies(n_theta, l_energy)
                    state = (j+1)//2

                    l_states.append(State(n_det, E, pt2, energies, state))

                step = Step(n_states, l_states)

                cipsi.append(step)

                i = try_go_to_match(buff, i, "Number of theta values:")

        l_cipsi.append(cipsi)

    data = []
    for cispi in l_cipsi:
        for step in cispi.steps:
            for state in step.states:
                print(state.n_det, state.E, state.pt2, state.state)
                for e in state.energies.energies:
                    print(e.theta, e.re, e.im, e.c_re, e.c_im)

                    data.append([state.state, e.theta, e.re, e.im, e.c_re, e.c_im, state.n_det, state.E, state.pt2])

    df = pd.DataFrame(data, columns=['state', 'theta', 'Re(E_cs)', 'Im(E_cs)', 'Re(Corr_E_cs)', 'Im(Corr_E_cs)', 'N_det', 'E', 'PT2'])

    print(df)
    N_det = df['N_det'].unique().tolist()
    states = df['state'].unique().tolist()
    val_theta = df['theta'].unique().tolist()

    print(N_det, states, val_theta)

    df = df[(df['N_det'] == N_det[-1])]

    data = np.zeros((n_theta, n_states*4+1), order='F')
    data[:, 0] = val_theta
    for i in range(n_states):
        df_s = df[(df['state'] == i+1)]
        cols = df_s[['Re(E_cs)', 'Im(E_cs)', 'Re(Corr_E_cs)', 'Im(Corr_E_cs)']].to_numpy()
        data[:, 1+i*4:1+(i+1)*4] = cols

    for row in data:
        row_str = ', '.join(map(str, row))
        print(row_str)
