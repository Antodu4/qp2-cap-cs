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
    """Print theta * dE/d(theta) for a complex CS energy E = y_re + i*y_im.

    For each theta value prints: theta, theta*d(Re E)/d(theta), theta*d(Im E)/d(theta),
    and theta * |dE/d(theta)|.

    Args:
        x    : 1-D array of theta values.
        y_re : 1-D array of real parts of the CS energy.
        y_im : 1-D array of imaginary parts of the CS energy.
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
    """Entry point: read a CS energy file and print theta * dE/d(theta).

    Usage: python deriv.py <file_path>

    The file must have columns: theta  Re(E)  Im(E)  ...
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
