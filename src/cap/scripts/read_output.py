#!/bin/python3
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

if __name__ == "__main__":
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    files = sys.argv[1:]

    l_cipsi = []
    for fname in files:
        cipsi = Cipsi([])
        with open(fname,'r') as f:
            buff = f.readlines()
            i = 0

            n_states = int(get_lines_match(fname,"Number of states:")[0].split()[-1])
            n_eta = int(get_lines_match(fname,"Number of eta values:")[0].split()[-1])

            first = True
            while i < len(buff):

                # Jump to the cap energies
                i = go_to_match(buff,i,"Number of eta values:")
                i += 8
                buff_eta = buff[i:i+n_eta]
                a = []
                for line in buff_eta:
                    for val in line.split():
                        a.append(val)
                a = np.asarray(a, dtype=float, order='F')
                E_eta = a.reshape((n_eta,2*n_states+1))

                # Jump to the energy corrections
                i = go_to_match(buff,i,"#### Energy corrections ###")
                i += 8
                buff_corr = buff[i:i+n_eta]
                a = []
                for line in buff_corr:
                    for val in line.split():
                        a.append(val)
                a = np.asarray(a, dtype=float, order='F')
                corr = a.reshape((n_eta,2*n_states+1))

                # N_det
                i = try_go_to_match(buff,i,"Summary at N_det =")
                if i < len(buff):
                    n_det = int(buff[i].split()[-1])
                else:
                    n_det = 0

                # E                    
                i = try_go_to_match(buff,i,"# E ")
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
                l_eta = E_eta[:,0]
                for j,E,pt2 in zip(range(1, 2*n_states+1, 2),l_E,l_pt2):
                    l_energy = []
                    l_re = E_eta[:,j]
                    l_im = E_eta[:,j+1]
                    l_corr_re = corr[:,j]
                    l_corr_im = corr[:,j+1]
                    for eta, re, im, c_re, c_im in zip(l_eta,l_re,l_im,l_corr_re,l_corr_im):
                        l_energy.append(Energy(float(eta),float(re),float(im),float(c_re),float(c_im)))
 
                    energies = Energies(n_eta,l_energy)
                    state = (j+1)//2

                    l_states.append(State(n_det,E,pt2,energies,state))

                step = Step(n_states, l_states)

                cipsi.append(step)

                i = try_go_to_match(buff,i,"Number of eta values:")

        l_cipsi.append(cipsi)

    data = []
    for cispi in l_cipsi:
        for step in cispi.steps:
            for state in step.states:
                print(state.n_det,state.E,state.pt2,state.state)
                for e in state.energies.energies:
                    print(e.eta, e.re, e.im, e.c_re, e.c_im)
                    
                    data.append([state.state, e.eta, e.re, e.im, e.c_re, e.c_im, state.n_det, state.E, state.pt2])

    df = pd.DataFrame(data, columns=['state','eta','Re(E_cap)','Im(E_cap)','Re(Corr_E_cap)','Im(Corr_E_cap)','N_det','E','PT2'])
    
    print(df)
    N_det = df['N_det'].unique().tolist()
    states = df['state'].unique().tolist()
    val_eta = df['eta'].unique().tolist()

    print(N_det, states, val_eta)

    df = df[(df['N_det'] == N_det[-1])]    

    #print(last_pt_df)

    #df.plot(kind = 'scatter', x = 'Re(E_cap)', y = 'Im(E_cap)')
    data = np.zeros((n_eta,n_states*4+1), order='F')
    data[:,0] = val_eta
    for i in range(n_states):
        df_s = df[(df['state'] == i+1)]
        cols = df_s[['Re(E_cap)','Im(E_cap)','Re(Corr_E_cap)','Im(Corr_E_cap)']].to_numpy()
        data[:,1+i*4:1+(i+1)*4] = cols

    for row in data:
        row_str = ', '.join(map(str, row))
        print(row_str)

