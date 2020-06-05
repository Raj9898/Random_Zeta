cimport numpy as nc

from random import uniform
from mpmath import zeta
from numpy import zeros, array
from math import log, sqrt, pi, cos


def zetaNormRange(int N , int deltaN) -> nc.ndarray:
    """
    An itterative algorithm for computing the normalized range for the Zeta function
    :param N: Provide a number (presumably large) to start computing values for
    :param deltaN: Provide an integer number to cap off the consecutive sums of prime
    :return: Returns an array of the normalized values for Zeta 
    """
    cdef nc.complex128_t point = complex(0, 0)
    cdef nc.complex128_t zeta_value = point
    cdef int t = N

    # storage for the normal zeta 
    cdef nc.ndarray zeta_array = zeros(shape=deltaN+1)
    
    # computing the Zeta function over a small interval defined by deltaN
    for t in range(N, N+deltaN+1):
        point = complex(real=0.5, imag=t)
  
        # computes the normalized zeta function -> refer to absoulte value code 
        zeta_array[t-N] = abs(zeta(point))

    return zeta_array


def zetaUniformNorm(int T, int h) -> double:
    """
    An algorithm for computing the Zeta function with random uniform component
    :param T: Provide a number (presumably large) to start computing values for
    :param h: Provide an integer number to modify the imaginary component
    :return: Returns an array of the normalized values for Zeta function with random uniform component
    """
    cdef nc.complex128_t point = complex(0, 0)
    cdef double upper = 2 * T

    # computing the Zeta function for a single point h for a given T
    point = zeta(complex(real=0.5, imag=uniform(T, upper) + h))

    return abs(point)


def zetaStochastic(double h, nc.ndarray prime_list) -> double:
    """
    An itterative algorithm for computing the stochastic Zeta function defined X_t(h)
    :param h: Provide a floating interval range to observe the Zeta function
    :param prime list: Provide an array of prime numbers
    :return: Returns a value for the stochastic Zeta function 
    """
    cdef double ret_val = 0.0
    cdef long i = 0
    cdef double upper = 2 * pi
    
    # itterate through the primes using the specified model 
    for p in prime_list: 
        ret_val += ( 1/sqrt(p) * cos(uniform(0, upper) - h*log(p)) )
   
    return ret_val
