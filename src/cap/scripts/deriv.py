import math
import numpy as np
import sys


def deriv(x, y):
    """Compute the numerical derivative of y with respect to x using finite differences.

    Uses a second-order forward difference at the first point, centered differences
    for interior points, and a second-order backward difference at the last point.

    Args:
        x : 1-D array of evenly spaced abscissa values.
        y : 1-D array of function values at x.

    Returns:
        df : 1-D array of derivative values, same length as y.
    """
    dx = x[1] - x[0]
    df = np.zeros((len(y)))

    df[0] = 1.0/(2.0*dx) * (-3.0 * y[0] + 4.0 * y[1] - y[2])
    for i in range(1, len(y)-1):
        df[i] = 1.0/(2.0 * dx) * (-y[i-1] + y[i+1])
    df[-1] = 1.0/(2.0*dx) * (y[-3] - 4.0 * y[-2] + 3.0 * y[-1])

    return df


def deriv_complex(x, y_re, y_im):
    """Print eta * dE/d(eta) for a complex energy E = y_re + i*y_im.

    For each eta value prints: eta, eta*d(Re E)/d(eta), eta*d(Im E)/d(eta),
    and eta * |dE/d(eta)|.

    Args:
        x    : 1-D array of eta values.
        y_re : 1-D array of real parts of the energy.
        y_im : 1-D array of imaginary parts of the energy.
    """
    d_re = deriv(x, y_re)
    d_im = deriv(x, y_im)

    for i in range(len(x)):
        print("{:8.3e} {:12.6f} {:12.6f} {:12.4f}".format(
            x[i], x[i] * d_re[i], x[i] * d_im[i],
            x[i] * abs(math.sqrt(d_re[i]**2 + d_im[i]**2))))


def read_text_file(file_path):
    """Read a whitespace-delimited text file into a list of float rows.

    Blank lines are skipped.

    Args:
        file_path : Path to the input file.

    Returns:
        data : List of lists of floats, one inner list per non-blank line.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    data = []
    for line in lines:
        if len(line.strip()) != 0:
            data.append([float(i) for i in line.split()])

    return data


def main():
    """Entry point: read a CAP energy file and print eta * dE/d(eta).

    Usage: python deriv.py <file_path>

    The file must have columns: eta  Re(E)  Im(E)  ...
    """
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]

    data = read_text_file(file_path)
    data = np.asarray(data, dtype=float)

    print(data)
    deriv_complex(data[:, 0], data[:, 1], data[:, 2])


if __name__ == "__main__":
    main()
