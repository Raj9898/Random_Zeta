
cimport numpy as nc
from libc.math cimport log, sqrt, pi, cos

from random import uniform

# Model is developed above 
def zetaStochastic(double h, long T, nc.ndarray prime_list) -> double:
    """
    An itterative algorithm for computing the stochastic Zeta function defined X_t(h)
    :param h: Provide a floating interval range to observe the Zeta function
    :param T: Provide an integer number to cap off the consecutive sum of primes
    :param prime list: Provide an array of prime numbers
    :return: Returns a value for the stochastic Zeta function 
    """
    assert T < 10000000000, 'Our prime list does not exceed the value 10^10'
    
    cdef double ret_val = 0.0
    cdef long p = 0
    cdef double upper = 2 * pi
    
    # filter the primes up to and including the value of T
    cdef nc.ndarray primes = prime_list[prime_list <= T]
    
    # itterate through the primes using the specified model 
    for p in primes: 
        ret_val += ( 1/sqrt(p) * cos(uniform(0, upper) - h*log(p)) )
   
    return ret_val