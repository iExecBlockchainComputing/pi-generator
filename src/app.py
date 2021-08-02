import os
import sys
import json
import math
from decimal import *

iexec_out = os.environ['IEXEC_OUT']


def sqrt(n, one):
    """
    Return the square root of n as a fixed point number with the one
    passed in.  It uses a second order Newton-Raphson convergence.  This
    doubles the number of significant figures on each iteration.
    """
    # Use floating point arithmetic to make an initial guess
    floating_point_precision = 10**16
    n_float = float((n * floating_point_precision) // one) / floating_point_precision
    x = (int(floating_point_precision * math.sqrt(n_float)) * one) // floating_point_precision
    n_one = n * one
    while 1:
        x_old = x
        x = (x + n_one // x) // 2
        if x == x_old:
            break
    return x


def pi_chudnovsky_bs(digits):
    """
    Compute int(pi * 10**digits)

    This is done using Chudnovsky's series with binary splitting
    """
    c = 640320
    c3_over_24 = c**3 // 24

    def bs(a, b):
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        if b - a == 1:
            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                pab = qab = 1
            else:
                pab = (6*a-5)*(2*a-1)*(6*a-1)
                qab = a*a*a*c3_over_24
            tab = pab * (13591409 + 545140134*a)  # a(a) * p(a)
            if a & 1:
                tab = -tab
        else:
            # Recursively compute P(a,b), Q(a,b) and T(a,b)
            # m is the midpoint of a and b
            m = (a + b) // 2
            # Recursively calculate P(a,m), Q(a,m) and T(a,m)
            pam, qam, tam = bs(a, m)
            # Recursively calculate P(m,b), Q(m,b) and T(m,b)
            pmb, qmb, tmb = bs(m, b)
            # Now combine
            pab = pam * pmb
            qab = qam * qmb
            tab = qmb * tam + pam * tmb
        return pab, qab, tab
    # how many terms to compute
    digits_per_term = math.log10(c3_over_24/6/2/6)
    n = int(digits/digits_per_term + 1)
    # Calclate P(0,N) and Q(0,N)
    p, q, t = bs(0, n)
    one = 10**digits
    sqrt_c = sqrt(10005*one, one)
    return (q*426880*sqrt_c) // t / Decimal(10 ** digits)


if len(sys.argv) != 2 or not sys.argv[1].isdigit():
    print('Please input the digit ordinal you need as an integer.')
    exit(1)

digits_number = int(sys.argv[1])
getcontext().prec = digits_number

pi = pi_chudnovsky_bs(digits_number)
print(pi)

# Append some results in /iexec_out/
with open(iexec_out + '/result.txt', 'w+') as fout:
    fout.write(str(pi))

# Declare everything is computed
with open(iexec_out + '/computed.json', 'w+') as f:
    json.dump({"deterministic-output-path": iexec_out + '/result.txt'}, f)
