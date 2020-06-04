
cimport numpy as nc
from libc.math cimport log, sqrt, pi, cos

from random import uniform
from mpmath import zeta
from numpy import zeros, array
from math import log, sqrt, pi, cos


def zetaNormRange(int N , int deltaN) -> nc.ndarray:
	"""
	An iterative algorithm for computing the normalized range for the Zeta function
	:param N: Provide a number (presumably large) to start computing values for
	:param deltaN: Provide an integer number to cap off the consecutive sums of prime
	:return: Returns an array of the normalized values for Zeta 
	"""
	assert N < int(1e12), 'Warning: computational efficiency slows beyond this point'
	assert deltaN < int(1e3), 'Warning: computational efficiency slows beyond this point'

	cdef nc.complex128_t point = complex(0, 0)
	cdef nc.complex128_t zeta_value = point
	cdef int t = N

	# storage for the normal zeta 
	cdef nc.ndarray zeta_array = zeros(shape = deltaN + 1)

	# computing the Zeta function over a small interval defined by deltaN
	for t in range(N, N + deltaN + 1):
	    point = complex(real=0.5, imag=t)

	    # computes the normalized zeta function -> refer to absoulte value code 
	    zeta_array[t - N] = abs(zeta(point))

	return zeta_array


def zetaStochastic(double h, long T, nc.ndarray prime_list) -> double:
	"""
	An iterative algorithm for computing the stochastic Zeta function defined X_t(h)
	:param h: Provide a floating interval range to observe the Zeta function
	:param T: Provide an integer number to cap off the consecutive sum of primes
	:param prime list: Provide an array of prime numbers
	:return: Returns a value for the stochastic Zeta function 
	"""
	assert T < int(1e10), 'Our prime list does not exceed the value 10^10'

	cdef double ret_val = 0.0
	cdef long p = 0
	cdef double upper = 2 * pi

	# filter the primes up to and including the value of T
	cdef nc.ndarray primes = prime_list[prime_list <= T]

	# iterate through the primes using the specified model
	for p in primes: 
	    ret_val += ( 1/sqrt(p) * cos(uniform(0, upper) - h*log(p)) )

	return ret_val


def zetaUniformNorm(int T, int h) -> double:
	"""
	Determines the values of the Zeta function with uniform parameters 
	:param T: Provide a number (presumably large) to start computing values for
	:param h: Provide a value h to be used in the complex axis
	:return: Returns a value from the zeta function at a given complex point 
	"""
	assert T < int(1e12), 'Warning: computational efficiency slows beyond this point'

	cdef nc.complex128_t point = complex(0, 0)
	cdef double upper = 2 * T

	# computing the Zeta function for a single point h for a given T
	point = zeta(complex(real=0.5, imag=uniform(T, upper) + h))

	return abs(point)