#!/usr/bin/env python3
"""
convert TREXIO file to EZFIO

Usage:
    qp_import_trexio [-o EZFIO_DIR] FILE

Options:
    -o --output=EZFIO_DIR    Produced directory
                             by default is FILE.ezfio

"""

import sys
import os
import numpy as np
from functools import reduce
from ezfio import ezfio
from docopt import docopt
import qp_bitmasks

try:
  import trexio
except ImportError:
    print("Error: trexio python module is not found. Try python3 -m pip install trexio")
    sys.exit(1)


try:
    QP_ROOT = os.environ["QP_ROOT"]
    QP_EZFIO = os.environ["QP_EZFIO"]
except KeyError:
    print("Error: QP_ROOT environment variable not found.")
    sys.exit(1)
else:
    sys.path = [QP_EZFIO + "/Python",
                QP_ROOT + "/install/resultsFile",
                QP_ROOT + "/install",
                QP_ROOT + "/scripts"] + sys.path

def uint64_to_int64(u):
    # Check if the most significant bit is set
    if u & (1 << 63):
        # Calculate the two's complement
        result = -int(np.bitwise_not(np.uint64(u))+1)
    else:
        # The number is already positive
        result = u
    return result

def _double_fact(n):
    """Double factorial: (-1)!!=1, 0!!=1, 1!!=1, 2!!=2, 3!!=3, ..."""
    if n <= 1:
        return 1
    return n * _double_fact(n - 2)


def _get_cart_to_sphe_T(l):
    """Return the cart_to_sphe_l matrix as a numpy array.

    T[j,k] is the contribution of the j-th Cartesian AO (generate_xyz /
    ao_power_index order) to the k-th unnormalized spherical AO
    (ordering: m=0, +1, -1, +2, -2, ..., +l, -l).
    Values copied from QP2's Fortran cart_to_sphe_l providers (Horton).
    Shape: (n_cart, n_sphe).
    """
    n_cart = (l + 1) * (l + 2) // 2
    n_sphe = 2 * l + 1
    T = np.zeros((n_cart, n_sphe))

    if l == 0:
        T[0, 0] = 1.0

    elif l == 1:
        # Cartesian: x(0), y(1), z(2)
        T[2, 0] = 1.0   # z  → m=0
        T[0, 1] = 1.0   # x  → m=+1
        T[1, 2] = 1.0   # y  → m=-1

    elif l == 2:
        # Cartesian: xx(0),xy(1),xz(2),yy(3),yz(4),zz(5)
        T[0, 0] = -0.5;                       T[3, 0] = -0.5;                       T[5, 0] = 1.0
        T[2, 1] = 1.0
        T[4, 2] = 1.0
        T[0, 3] = 0.86602540378443864676;     T[3, 3] = -0.86602540378443864676
        T[1, 4] = 1.0

    elif l == 3:
        # Cartesian: xxx(0),xxy(1),xxz(2),xyy(3),xyz(4),xzz(5),yyy(6),yyz(7),yzz(8),zzz(9)
        T[2, 0] = -0.67082039324993690892;    T[7, 0] = -0.67082039324993690892;    T[9, 0] = 1.0
        T[0, 1] = -0.61237243569579452455;    T[3, 1] = -0.27386127875258305673;    T[5, 1] = 1.0954451150103322269
        T[1, 2] = -0.27386127875258305673;    T[6, 2] = -0.61237243569579452455;    T[8, 2] = 1.0954451150103322269
        T[2, 3] = 0.86602540378443864676;     T[7, 3] = -0.86602540378443864676
        T[4, 4] = 1.0
        T[0, 5] = 0.790569415042094833;       T[3, 5] = -1.0606601717798212866
        T[1, 6] = 1.0606601717798212866;      T[6, 6] = -0.790569415042094833

    elif l == 4:
        # Cartesian: xxxx(0),xxxy(1),xxxz(2),xxyy(3),xxyz(4),xxzz(5),
        #            xyyy(6),xyyz(7),xyzz(8),xzzz(9),yyyy(10),yyyz(11),yyzz(12),yzzz(13),zzzz(14)
        T[0,  0] = 0.375;                     T[3,  0] = 0.21957751641341996535
        T[5,  0] = -0.87831006565367986142;   T[10, 0] = 0.375
        T[12, 0] = -0.87831006565367986142;   T[14, 0] = 1.0
        T[2,  1] = -0.89642145700079522998;   T[7,  1] = -0.40089186286863657703;   T[9,  1] = 1.19522860933439364
        T[4,  2] = -0.40089186286863657703;   T[11, 2] = -0.89642145700079522998;   T[13, 2] = 1.19522860933439364
        T[0,  3] = -0.5590169943749474241;    T[5,  3] = 0.9819805060619657157
        T[10, 3] = 0.5590169943749474241;     T[12, 3] = -0.9819805060619657157
        T[1,  4] = -0.42257712736425828875;   T[6,  4] = -0.42257712736425828875;   T[8,  4] = 1.1338934190276816816
        T[2,  5] = 0.790569415042094833;      T[7,  5] = -1.0606601717798212866
        T[4,  6] = 1.0606601717798212866;     T[11, 6] = -0.790569415042094833
        T[0,  7] = 0.73950997288745200532;    T[3,  7] = -1.2990381056766579701;    T[10, 7] = 0.73950997288745200532
        T[1,  8] = 1.1180339887498948482;     T[6,  8] = -1.1180339887498948482

    elif l == 5:
        # Cartesian: 21 functions (xxxxx..zzzzz in ao_power_index order)
        T[2,  0] = 0.625;                     T[7,  0] = 0.36596252735569994226
        T[9,  0] = -1.0910894511799619063;    T[16, 0] = 0.625
        T[18, 0] = -1.0910894511799619063;    T[20, 0] = 1.0
        T[0,  1] = 0.48412291827592711065;    T[3,  1] = 0.21128856368212914438
        T[5,  1] = -1.2677313820927748663;    T[10, 1] = 0.16137430609197570355
        T[12, 1] = -0.56694670951384084082;   T[14, 1] = 1.2909944487358056284
        T[1,  2] = 0.16137430609197570355;    T[6,  2] = 0.21128856368212914438
        T[8,  2] = -0.56694670951384084082;   T[15, 2] = 0.48412291827592711065
        T[17, 2] = -1.2677313820927748663;    T[19, 2] = 1.2909944487358056284
        T[2,  3] = -0.85391256382996653194;   T[9,  3] = 1.1180339887498948482
        T[16, 3] = 0.85391256382996653194;    T[18, 3] = -1.1180339887498948482
        T[4,  4] = -0.6454972243679028142;    T[11, 4] = -0.6454972243679028142;    T[13, 4] = 1.2909944487358056284
        T[0,  5] = -0.52291251658379721749;   T[3,  5] = 0.22821773229381921394
        T[5,  5] = 0.91287092917527685576;    T[10, 5] = 0.52291251658379721749;    T[12, 5] = -1.2247448713915890491
        T[1,  6] = -0.52291251658379721749;   T[6,  6] = -0.22821773229381921394
        T[8,  6] = 1.2247448713915890491;     T[15, 6] = 0.52291251658379721749;    T[17, 6] = -0.91287092917527685576
        T[2,  7] = 0.73950997288745200532;    T[7,  7] = -1.2990381056766579701;    T[16, 7] = 0.73950997288745200532
        T[4,  8] = 1.1180339887498948482;     T[11, 8] = -1.1180339887498948482
        T[0,  9] = 0.7015607600201140098;     T[3,  9] = -1.5309310892394863114;    T[10, 9] = 1.169267933366856683
        T[1, 10] = 1.169267933366856683;      T[6, 10] = -1.5309310892394863114;    T[15,10] = 0.7015607600201140098

    elif l == 6:
        # Cartesian: 28 functions
        T[0,  0] = -0.3125;                   T[3,  0] = -0.16319780245846672329
        T[5,  0] = 0.97918681475080033975;    T[10, 0] = -0.16319780245846672329
        T[12, 0] = 0.57335309036732873772;    T[14, 0] = -1.3055824196677337863
        T[21, 0] = -0.3125;                   T[23, 0] = 0.97918681475080033975
        T[25, 0] = -1.3055824196677337863;    T[27, 0] = 1.0
        T[2,  1] = 0.86356159963469679725;    T[7,  1] = 0.37688918072220452831
        T[9,  1] = -1.6854996561581052156;    T[16, 1] = 0.28785386654489893242
        T[18, 1] = -0.75377836144440905662;   T[20, 1] = 1.3816985594155148756
        T[4,  2] = 0.28785386654489893242;    T[11, 2] = 0.37688918072220452831
        T[13, 2] = -0.75377836144440905662;   T[22, 2] = 0.86356159963469679725
        T[24, 2] = -1.6854996561581052156;    T[26, 2] = 1.3816985594155148756
        T[0,  3] = 0.45285552331841995543;    T[3,  3] = 0.078832027985861408788
        T[5,  3] = -1.2613124477737825406;    T[10, 3] = -0.078832027985861408788
        T[14, 3] = 1.2613124477737825406;     T[21, 3] = -0.45285552331841995543
        T[23, 3] = 1.2613124477737825406;     T[25, 3] = -1.2613124477737825406
        T[1,  4] = 0.27308215547040717681;    T[6,  4] = 0.26650089544451304287
        T[8,  4] = -0.95346258924559231545;   T[15, 4] = 0.27308215547040717681
        T[17, 4] = -0.95346258924559231545;   T[19, 4] = 1.4564381625088382763
        T[2,  5] = -0.81924646641122153043;   T[7,  5] = 0.35754847096709711829
        T[9,  5] = 1.0660035817780521715;     T[16, 5] = 0.81924646641122153043
        T[18, 5] = -1.4301938838683884732
        T[4,  6] = -0.81924646641122153043;   T[11, 6] = -0.35754847096709711829
        T[13, 6] = 1.4301938838683884732;     T[22, 6] = 0.81924646641122153043
        T[24, 6] = -1.0660035817780521715
        T[0,  7] = -0.49607837082461073572;   T[3,  7] = 0.43178079981734839863
        T[5,  7] = 0.86356159963469679725;    T[10, 7] = 0.43178079981734839863
        T[12, 7] = -1.5169496905422946941;    T[21, 7] = -0.49607837082461073572
        T[23, 7] = 0.86356159963469679725
        T[1,  8] = -0.59829302641309923139;   T[8,  8] = 1.3055824196677337863
        T[15, 8] = 0.59829302641309923139;    T[17, 8] = -1.3055824196677337863
        T[2,  9] = 0.7015607600201140098;     T[7,  9] = -1.5309310892394863114;    T[16, 9] = 1.169267933366856683
        T[4, 10] = 1.169267933366856683;      T[11,10] = -1.5309310892394863114;    T[22,10] = 0.7015607600201140098
        T[0, 11] = 0.67169328938139615748;    T[3, 11] = -1.7539019000502850245
        T[10,11] = 1.7539019000502850245;     T[21,11] = -0.67169328938139615748
        T[1, 12] = 1.2151388809514737933;     T[6, 12] = -1.9764235376052370825;    T[15,12] = 1.2151388809514737933

    else:
        raise ValueError(f"cart_to_sphe not implemented for l={l} > 6")

    return T


def _compute_T_eff(T, cart_powers):
    """Compute the normalized transformation T_eff[j,k] = T[j,k] * N_cart[j] / N_sphe[k].

    With this matrix: phi_sphe_norm[k] = sum_j T_eff[j,k] * phi_cart_norm[j],
    so normalized spherical MO coefficients transform as C_cart = T_eff @ C_sphe.
    """
    n_cart, n_sphe = T.shape

    # Angular self-overlap ∝ (2px-1)!! * (2py-1)!! * (2pz-1)!!
    w_cart = np.array([float(_double_fact(2*p[0]-1) * _double_fact(2*p[1]-1) * _double_fact(2*p[2]-1))
                       for p in cart_powers])

    # Cross-overlap of unnormalized Cartesian GTOs (angular part only)
    S = np.zeros((n_cart, n_cart))
    for j, (px1, py1, pz1) in enumerate(cart_powers):
        for jp, (px2, py2, pz2) in enumerate(cart_powers):
            if (px1+px2) % 2 == 0 and (py1+py2) % 2 == 0 and (pz1+pz2) % 2 == 0:
                S[j, jp] = float(_double_fact(px1+px2-1) * _double_fact(py1+py2-1) * _double_fact(pz1+pz2-1))

    # Self-overlap of the k-th unnormalized spherical GTO: T[:,k]^T @ S @ T[:,k]
    w_sphe = np.array([T[:, k] @ S @ T[:, k] for k in range(n_sphe)])

    return T * np.sqrt(w_cart)[:, np.newaxis] / np.sqrt(w_sphe)[np.newaxis, :]


def _sphe_to_cart_mo_coefs(MoMatrix_sphe, ang_mom_per_shell):
    """Transform MO coefficients from spherical to Cartesian AO basis.

    TREXIO/EZFIO convention: shape (mo_num, ao_num) on both input and output.
    Spherical ordering: m = 0, +1, -1, +2, -2, ..., +l, -l
    Cartesian ordering: generate_xyz (= ao_power_index) order per shell

    Returns (MoMatrix_cart, ao_cart_num).
    """
    cart_blocks = []
    sphe_idx = 0

    for l in ang_mom_per_shell:
        n_cart = (l + 1) * (l + 2) // 2
        n_sphe = 2 * l + 1
        cart_powers = generate_xyz(l)
        T = _get_cart_to_sphe_T(l)
        T_eff = _compute_T_eff(T, cart_powers)
        # (mo_num, n_sphe) @ (n_sphe, n_cart) → (mo_num, n_cart)
        cart_blocks.append(MoMatrix_sphe[:, sphe_idx:sphe_idx + n_sphe] @ T_eff.T)
        sphe_idx += n_sphe

    ao_cart_num = sum((l + 1) * (l + 2) // 2 for l in ang_mom_per_shell)
    return np.hstack(cart_blocks), ao_cart_num


def generate_xyz(l):

    def create_z(x,y,z):
       return (x, y, l-(x+y))

    def create_y(accu,x,y,z):
       if y == 0:
          result = [create_z(x,y,z)] + accu
       else:
          result = create_y([create_z(x,y,z)] + accu , x, y-1, z)
       return result

    def create_x(accu,x,y,z):
       if x == 0:
          result = create_y([], x,y,z) + accu
       else:
          xnew = x-1
          ynew = l-xnew
          result = create_x(create_y([],x,y,z) + accu , xnew, ynew, z)
       return result

    result = create_x([], l, 0, 0)
    result.reverse()
    return result



def write_ezfio(trexio_filename, filename):

    try:
        trexio_file = trexio.File(trexio_filename,mode='r',back_end=trexio.TREXIO_TEXT)
    except:
        trexio_file = trexio.File(trexio_filename,mode='r',back_end=trexio.TREXIO_HDF5)

    ezfio.set_file(filename)
    ezfio.set_trexio_trexio_file(trexio_filename)

    print("Nuclei\t\t...\t", end=' ')

    charge = [0.]
    if trexio.has_nucleus(trexio_file):
        charge = trexio.read_nucleus_charge(trexio_file)
        ezfio.set_nuclei_nucl_num(len(charge))
        ezfio.set_nuclei_nucl_charge(charge)

        coord = trexio.read_nucleus_coord(trexio_file)
        coord = np.transpose(coord)
        ezfio.set_nuclei_nucl_coord(coord)

        label = trexio.read_nucleus_label(trexio_file)
        nucl_num = trexio.read_nucleus_num(trexio_file)

        # Transformt H1 into H
        import re
        p = re.compile(r'(\d*)$')
        label = [p.sub("", x).capitalize() for x in label]
        ezfio.set_nuclei_nucl_label(label)
        print("OK")

    else:
        ezfio.set_nuclei_nucl_num(1)
        ezfio.set_nuclei_nucl_charge([0.])
        ezfio.set_nuclei_nucl_coord([0.,0.,0.])
        ezfio.set_nuclei_nucl_label(["X"])
        print("None")



    print("Electrons\t...\t", end=' ')

    try:
        num_beta = trexio.read_electron_dn_num(trexio_file)
    except:
        num_beta = int(sum(charge))//2

    try:
        num_alpha = trexio.read_electron_up_num(trexio_file)
    except:
        num_alpha = int(sum(charge)) - num_beta

    if num_alpha == 0:
        print("\n\nError: There are zero electrons in the TREXIO file.\n\n")
        sys.exit(1)
    ezfio.set_electrons_elec_alpha_num(num_alpha)
    ezfio.set_electrons_elec_beta_num(num_beta)

    print(f"{num_alpha} {num_beta}")

    print("Basis\t\t...\t", end=' ')

    shell_num = 0
    try:
        basis_type = trexio.read_basis_type(trexio_file)

        print ("BASIS TYPE: ", basis_type.lower())
        if basis_type.lower() in ["gaussian", "slater"]:
            shell_num   = trexio.read_basis_shell_num(trexio_file)
            prim_num    = trexio.read_basis_prim_num(trexio_file)
            ang_mom     = trexio.read_basis_shell_ang_mom(trexio_file)
            nucl_index  = trexio.read_basis_nucleus_index(trexio_file)
            exponent    = trexio.read_basis_exponent(trexio_file)
            coefficient = trexio.read_basis_coefficient(trexio_file)
            shell_index = trexio.read_basis_shell_index(trexio_file)
            ao_shell    = trexio.read_ao_shell(trexio_file)

            ezfio.set_basis_basis("Read from TREXIO")
            ezfio.set_ao_basis_ao_basis("Read from TREXIO")
            ezfio.set_basis_shell_num(shell_num)
            ezfio.set_basis_prim_num(prim_num)
            ezfio.set_basis_shell_ang_mom(ang_mom)
            ezfio.set_basis_basis_nucleus_index([ x+1 for x in nucl_index ])
            ezfio.set_basis_prim_expo(exponent)
            ezfio.set_basis_prim_coef(coefficient)

            nucl_shell_num = []
            prev = None
            m = 0
            for i in ao_shell:
                if i != prev:
                   m += 1
                   if prev is None or nucl_index[i] != nucl_index[prev]:
                        nucl_shell_num.append(m)
                        m = 0
                prev = i
            assert (len(nucl_shell_num) == nucl_num)

            shell_prim_num = []
            prev = shell_index[0]
            count = 0
            for i in shell_index:
                if i != prev:
                   shell_prim_num.append(count)
                   count = 0
                count += 1
                prev = i
            shell_prim_num.append(count)

            assert (len(shell_prim_num) == shell_num)

            ezfio.set_basis_shell_prim_num(shell_prim_num)
            ezfio.set_basis_shell_index([x+1 for x in shell_index])
            ezfio.set_basis_nucleus_shell_num(nucl_shell_num)


            shell_factor = trexio.read_basis_shell_factor(trexio_file)
            prim_factor  = trexio.read_basis_prim_factor(trexio_file)

        elif basis_type.lower() == "numerical":

            shell_num   = trexio.read_basis_shell_num(trexio_file)
            prim_num    = shell_num
            ang_mom     = trexio.read_basis_shell_ang_mom(trexio_file)
            nucl_index  = trexio.read_basis_nucleus_index(trexio_file)
            exponent    = [1.]*prim_num
            coefficient = [1.]*prim_num
            shell_index = [i for i in range(shell_num)]
            ao_shell    = trexio.read_ao_shell(trexio_file)

            ezfio.set_basis_basis("None")
            ezfio.set_ao_basis_ao_basis("None")
            ezfio.set_basis_shell_num(shell_num)
            ezfio.set_basis_prim_num(prim_num)
            ezfio.set_basis_shell_ang_mom(ang_mom)
            ezfio.set_basis_basis_nucleus_index([ x+1 for x in nucl_index ])
            ezfio.set_basis_prim_expo(exponent)
            ezfio.set_basis_prim_coef(coefficient)

            nucl_shell_num = []
            prev = None
            m = 0
            for i in ao_shell:
                if i != prev:
                   m += 1
                   if prev is None or nucl_index[i] != nucl_index[prev]:
                        nucl_shell_num.append(m)
                        m = 0
                prev = i
            assert (len(nucl_shell_num) == nucl_num)

            shell_prim_num = []
            prev = shell_index[0]
            count = 0
            for i in shell_index:
                if i != prev:
                   shell_prim_num.append(count)
                   count = 0
                count += 1
                prev = i
            shell_prim_num.append(count)

            assert (len(shell_prim_num) == shell_num)

            ezfio.set_basis_shell_prim_num(shell_prim_num)
            ezfio.set_basis_shell_index([x+1 for x in shell_index])
            ezfio.set_basis_nucleus_shell_num(nucl_shell_num)

            shell_factor = trexio.read_basis_shell_factor(trexio_file)
            prim_factor  = [1.]*prim_num
        else:
           raise TypeError

        print(basis_type)
    except:
        print("None")
        ezfio.set_ao_basis_ao_cartesian(True)

    print("AOS\t\t...\t", end=' ')

    try:
        cartesian = trexio.read_ao_cartesian(trexio_file)
    except:
        cartesian = True

    ao_num = trexio.read_ao_num(trexio_file)
    ezfio.set_ao_basis_ao_num(ao_num)

    sphe_mo_ang_mom = None  # set in the spherical branch; used later for MO transform

    if cartesian and shell_num > 0:
        ao_shell    = trexio.read_ao_shell(trexio_file)
        at = [ nucl_index[i]+1 for i in ao_shell ]
        ezfio.set_ao_basis_ao_nucl(at)

        num_prim0 = [ 0 for i in range(shell_num) ]
        for i in shell_index:
            num_prim0[i] += 1

        coef = {}
        expo = {}
        for i,c in enumerate(coefficient):
            idx = shell_index[i]
            if idx in coef:
              coef[idx].append(c)
              expo[idx].append(exponent[i])
            else:
              coef[idx] = [c]
              expo[idx] = [exponent[i]]

        coefficient = []
        exponent    = []
        power_x     = []
        power_y     = []
        power_z     = []
        num_prim    = []

        for i in range(shell_num):
            for x,y,z in generate_xyz(ang_mom[i]):
                power_x.append(x)
                power_y.append(y)
                power_z.append(z)
                coefficient.append(coef[i])
                exponent.append(expo[i])
                num_prim.append(num_prim0[i])

        assert (len(coefficient) == ao_num)
        ezfio.set_ao_basis_ao_power(power_x + power_y + power_z)
        ezfio.set_ao_basis_ao_prim_num(num_prim)

        prim_num_max = max( [ len(x) for x in coefficient ] )

        for i in range(ao_num):
            coefficient[i] += [0. for j in range(len(coefficient[i]), prim_num_max)]
            exponent   [i] += [0. for j in range(len(exponent[i]), prim_num_max)]

        coefficient = reduce(lambda x, y: x + y, coefficient, [])
        exponent    = reduce(lambda x, y: x + y, exponent   , [])

        coef = []
        expo = []
        for i in range(prim_num_max):
            for j in range(i, len(coefficient), prim_num_max):
                coef.append(coefficient[j])
                expo.append(exponent[j])

#        ezfio.set_ao_basis_ao_prim_num_max(prim_num_max)
        ezfio.set_ao_basis_ao_coef(coef)
        ezfio.set_ao_basis_ao_expo(expo)

        print("OK")

    elif not cartesian and shell_num > 0:
        # Spherical basis (ORCA default): expand shells to Cartesian AOs so that
        # QP2 can compute integrals. ao_cartesian stays False (default), which
        # tells QP2 to use the cart→sphe transform internally.
        ao_shell = trexio.read_ao_shell(trexio_file)
        at = [ nucl_index[i]+1 for i in ao_shell ]

        num_prim0 = [ 0 for i in range(shell_num) ]
        for i in shell_index:
            num_prim0[i] += 1

        coef = {}
        expo = {}
        for i, c in enumerate(coefficient):
            idx = shell_index[i]
            if idx in coef:
                coef[idx].append(c)
                expo[idx].append(exponent[i])
            else:
                coef[idx] = [c]
                expo[idx] = [exponent[i]]

        coefficient_cart = []
        exponent_cart    = []
        power_x          = []
        power_y          = []
        power_z          = []
        num_prim         = []
        ao_nucl_cart     = []
        sphe_mo_ang_mom  = []   # angular momentum per shell, for MO transform

        for i in range(shell_num):
            l = ang_mom[i]
            sphe_mo_ang_mom.append(l)
            for x, y, z in generate_xyz(l):
                power_x.append(x)
                power_y.append(y)
                power_z.append(z)
                coefficient_cart.append(coef[i])
                exponent_cart.append(expo[i])
                num_prim.append(num_prim0[i])
                ao_nucl_cart.append(nucl_index[i] + 1)

        ao_num_cart = len(coefficient_cart)
        ezfio.set_ao_basis_ao_num(ao_num_cart)
        ezfio.set_ao_basis_ao_nucl(ao_nucl_cart)
        ezfio.set_ao_basis_ao_power(power_x + power_y + power_z)
        ezfio.set_ao_basis_ao_prim_num(num_prim)

        prim_num_max = max(len(x) for x in coefficient_cart)

        for i in range(ao_num_cart):
            coefficient_cart[i] += [0.] * (prim_num_max - len(coefficient_cart[i]))
            exponent_cart[i]    += [0.] * (prim_num_max - len(exponent_cart[i]))

        flat_coef = reduce(lambda x, y: x + y, coefficient_cart, [])
        flat_expo = reduce(lambda x, y: x + y, exponent_cart,    [])

        coef_out = []
        expo_out = []
        for i in range(prim_num_max):
            for j in range(i, len(flat_coef), prim_num_max):
                coef_out.append(flat_coef[j])
                expo_out.append(flat_expo[j])

        ezfio.set_ao_basis_ao_coef(coef_out)
        ezfio.set_ao_basis_ao_expo(expo_out)
        print("OK (spherical basis, Cartesian AO storage)")

    else:
        print("None: integrals should be also imported using qp run import_trexio_integrals")


    #                _
    # |\/|  _   _   |_)  _.  _ o  _
    # |  | (_) _>   |_) (_| _> | _>
    #

    print("MOS\t\t...\t", end=' ')

    labels = { "Canonical" : "Canonical",
               "RHF" : "Canonical",
               "BOYS" : "Localized",
               "ROHF" : "Canonical",
               "UHF" : "Canonical",
               "Natural": "Natural" }
    try:
      label = labels[trexio.read_mo_type(trexio_file)]
    except:
      label = "None"
    ezfio.set_mo_basis_mo_label(label)
    ezfio.set_determinants_mo_label(label)

    try:
      clss = trexio.read_mo_class(trexio_file)
      core     = [ i for i in clss if i.lower() == "core" ]
      inactive = [ i for i in clss if i.lower() == "inactive" ]
      active   = [ i for i in clss if i.lower() == "active" ]
      virtual  = [ i for i in clss if i.lower() == "virtual" ]
      deleted  = [ i for i in clss if i.lower() == "deleted" ]
    except trexio.Error:
      pass

    try:
      mo_num = trexio.read_mo_num(trexio_file)
      ezfio.set_mo_basis_mo_num(mo_num)

      MoMatrix = trexio.read_mo_coefficient(trexio_file)
      if sphe_mo_ang_mom is not None:
          MoMatrix, _ = _sphe_to_cart_mo_coefs(MoMatrix, sphe_mo_ang_mom)
      ezfio.set_mo_basis_mo_coef(MoMatrix)

      mo_occ = [ 0. for i in range(mo_num) ]
      for i in range(num_alpha):
         mo_occ[i] += 1.
      for i in range(num_beta):
         mo_occ[i] += 1.
      ezfio.set_mo_basis_mo_occ(mo_occ)
      print("OK")
    except:
      print("None")



    print("Pseudos\t\t...\t", end=' ')

    ezfio.set_pseudo_do_pseudo(False)

    if trexio.has_ecp_ang_mom(trexio_file):
        ezfio.set_pseudo_do_pseudo(True)
        max_ang_mom_plus_1 = trexio.read_ecp_max_ang_mom_plus_1(trexio_file)
        z_core = trexio.read_ecp_z_core(trexio_file)
        ang_mom = trexio.read_ecp_ang_mom(trexio_file)
        nucleus_index = trexio.read_ecp_nucleus_index(trexio_file)
        exponent = trexio.read_ecp_exponent(trexio_file)
        coefficient = trexio.read_ecp_coefficient(trexio_file)
        power = trexio.read_ecp_power(trexio_file)

        lmax = max( max_ang_mom_plus_1 ) - 1
        ezfio.set_pseudo_pseudo_lmax(lmax)
        ezfio.set_pseudo_nucl_charge_remove(z_core)

        prev_center = None
        ecp = {}
        for i in range(len(ang_mom)):
            center = nucleus_index[i]
            if center != prev_center:
               ecp[center] = { "lmax": max_ang_mom_plus_1[center],
                               "zcore": z_core[center],
                               "contr": {} }
               for j in range(max_ang_mom_plus_1[center]+1):
                    ecp[center]["contr"][j] = []

            ecp[center]["contr"][ang_mom[i]].append( (coefficient[i], power[i], exponent[i]) )
            prev_center = center

        ecp_loc = {}
        ecp_nl  = {}
        kmax    = 0
        klocmax    = 0
        for center in ecp:
            ecp_nl [center] = {}
            for k in ecp[center]["contr"]:
                if k == ecp[center]["lmax"]:
                    ecp_loc[center] = ecp[center]["contr"][k]
                    klocmax = max(len(ecp_loc[center]), klocmax)
                else:
                    ecp_nl [center][k] = ecp[center]["contr"][k]
                    kmax = max(len(ecp_nl [center][k]), kmax)

        ezfio.set_pseudo_pseudo_klocmax(klocmax)
        ezfio.set_pseudo_pseudo_kmax(kmax)

        pseudo_n_k = [[0  for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_v_k = [[0. for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_dz_k = [[0. for _ in range(nucl_num)] for _ in range(klocmax)]
        pseudo_n_kl = [[[0  for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        pseudo_v_kl = [[[0. for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        pseudo_dz_kl = [[[0. for _ in range(nucl_num)] for _ in range(kmax)] for _ in range(lmax+1)]
        for center in ecp_loc:
            for k in range( len(ecp_loc[center]) ):
                v, n, dz = ecp_loc[center][k]
                pseudo_n_k[k][center] = n
                pseudo_v_k[k][center] = v
                pseudo_dz_k[k][center] = dz

        ezfio.set_pseudo_pseudo_n_k(pseudo_n_k)
        ezfio.set_pseudo_pseudo_v_k(pseudo_v_k)
        ezfio.set_pseudo_pseudo_dz_k(pseudo_dz_k)

        for center in ecp_nl:
            for l in range( len(ecp_nl[center]) ):
                for k in range( len(ecp_nl[center][l]) ):
                    v, n, dz = ecp_nl[center][l][k]
                    pseudo_n_kl[l][k][center] = n
                    pseudo_v_kl[l][k][center] = v
                    pseudo_dz_kl[l][k][center] = dz

        ezfio.set_pseudo_pseudo_n_kl(pseudo_n_kl)
        ezfio.set_pseudo_pseudo_v_kl(pseudo_v_kl)
        ezfio.set_pseudo_pseudo_dz_kl(pseudo_dz_kl)
        print("OK")

    else:
        print("None")

    print("Determinant\t...\t", end=' ')
    alpha = [ i for i in range(num_alpha) ]
    beta  = [ i for i in range(num_beta) ]
    if trexio.has_mo_spin(trexio_file):
       spin = trexio.read_mo_spin(trexio_file)
       if max(spin) == 1:
         alpha = [ i for i in range(len(spin)) if spin[i] == 0 ]
         alpha = [ alpha[i] for i in range(num_alpha) ]
         beta  = [ i for i in range(len(spin)) if spin[i] == 1 ]
         beta  = [ beta[i] for i in range(num_beta) ]
         print("Warning -- UHF orbitals --", end=' ')
    alpha_s = ['0']*mo_num
    beta_s  = ['0']*mo_num
    for i in alpha:
      alpha_s[i] = '1'
    for i in beta:
      beta_s[i] = '1'
    alpha_s = ''.join(alpha_s)[::-1]
    beta_s = ''.join(beta_s)[::-1]
    def conv(i):
      try:
        result = np.int64(i)
      except:
        result = np.int64(i-2**63-1)
      return result

    alpha = [ uint64_to_int64(int(i,2)) for i in qp_bitmasks.string_to_bitmask(alpha_s) ][::-1]
    beta  = [ uint64_to_int64(int(i,2)) for i in qp_bitmasks.string_to_bitmask(beta_s ) ][::-1]
    ezfio.set_determinants_bit_kind(8)
    ezfio.set_determinants_n_int(1+mo_num//64)
    ezfio.set_determinants_n_det(1)
    ezfio.set_determinants_n_states(1)
    ezfio.set_determinants_psi_det(alpha+beta)
    ezfio.set_determinants_psi_coef([[1.0]])
    print("OK")




def get_full_path(file_path):
    file_path = os.path.expanduser(file_path)
    file_path = os.path.expandvars(file_path)
    return file_path


if __name__ == '__main__':
    ARGUMENTS = docopt(__doc__)

    FILE = get_full_path(ARGUMENTS['FILE'])
    trexio_filename = FILE

    if ARGUMENTS["--output"]:
        EZFIO_FILE = get_full_path(ARGUMENTS["--output"])
    else:
        EZFIO_FILE = "{0}.ezfio".format(FILE)

    write_ezfio(trexio_filename, EZFIO_FILE)
    sys.stdout.flush()

