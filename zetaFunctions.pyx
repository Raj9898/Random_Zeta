"""
Name: Rajesh Rao
Riemann Zeta Function
Professor: Louis-Pierre Arguin
Version 1.0
--
Cython programs for efficiently computing the various forms of the zeta function
for use in research of the large value conjectures
"""

cimport numpy as nc

from mpmath import zeta
from math import cos, sqrt, pi, log
from random import uniform


def zetaNorm(nc.complex128_t point) -> double:
    """
    Computes the zeta norm at a given complex point
    :param i: Provide a complex number to comute the zeta norm
    :return: Returns a float, defined as the zeta norm
    """
    cdef double zetaN = abs(zeta(point))

    return zetaN


def logZetaNorm(nc.complex128_t i) -> double:
    """
    Computes the log zeta norm at a given complex point
    :param i: Provide a complex number to comute the log zeta norm
    :return: Returns a float, defined as the log zeta norm
    """
    cdef nc.complex128_t point = i

    return log(abs(zeta(point)))


def zetaStochastic(double h, nc.ndarray prime_list) -> double:
    """
    An iterative algorithm for computing the stochastic Zeta function defined X_t(h)
    :param h: Provide a floating interval range to observe the Zeta function
    :param prime list: Provide an array of prime numbers
    :return: Returns a value for the stochastic Zeta function 
    """
    cdef double ret_val = 0.0
    cdef long i = 0
    cdef double upper = 2 * pi
    
    # iterate through the primes using the specified model
    for p in prime_list: 
        ret_val += ( 1/sqrt(p) * cos(uniform(0, upper) - h*log(p)) )
   
    return ret_val


def zetaUniformNorm(int T, int h) -> double:
    """
    An algorithm for computing the Zeta function with random uniform component
    :param T: Provide a number (presumably large) to start computing values for
    :param h: Provide an integer number to modify the imaginary component
    :return: Returns an value for the normalized values for Zeta function with random uniform component
    """
    cdef nc.complex128_t point = complex(0, 0)
    cdef double upper = 2 * T

    # computing the Zeta function for a single point h for a given T
    point = zeta(complex(real=0.5, imag=uniform(T, upper) + h))

    return abs(point)