#########################################################################################
#########################################################################################
# Zeta Function 
# NSF CAREER grant #1653602 : Statistics of Extrema in Complex and Disordered Systems
# Prof. Louis-Pierre Arguin
# ---------------------------------------------------------------------------------------
# Authors: Rajesh Rao, Kwokching (Kelvin) Hui, Eli Amzallag 
#########################################################################################
#########################################################################################


#################
# Imports
#################
import numpy as np
import mpmath as mp
from typing import Union


#################
# Function
#################
def zetaNaive(x: complex) -> complex: 
  """
  Naive implementation of the Riemann Zeta Function with the restricted sigma condition at a single point
  :param x: Provide a complex number, with real component strictly greater than 1
  :return: Returns a complex number from the Zeta function 
  """
  assert x.real > 1, "Real component of the complex component must be strictly greater than 1"
  zetaSum = 0.0
  
  # computes the zeta function defined Z = 1^-x + 2^-x + 3^-x + ...
  for i in range(1, int(1e6)):
      zetaSum += 1/(i**x)
      
  return zetaSum


def zetaStochastic(h : float, T: int, prime_list: np.array = None, upper: double = np.pi) -> float:
  """
  An itterative algorithm for computing the stochastic Zeta function defined X_t(h)
  :param h: Provide a floating interval range to observe the Zeta function
  :param T: Provide an integer number to cap off the consecutive sum of primes
  :return: Returns a value for the stochastic Zeta function 
  """
  assert T < int(1e10), 'Our prime list does not exceed the value 10^7'
  ret_val = 0.0
  
  # filter the primes up to and including the value of T
  primes = prime_list[prime_list <= T]
  ret_val = sum([p**(-0.5) * np.cos(np.random.uniform(0, 2*upper) - h*np.log(p)) for p in primes])

  return ret_val


def zetaNormRange(N : int, deltaN : int) -> np.array:
  """
  An itterative algorithm for computing the stochastic Zeta function defined X_t(h)
  :param N: Provide a number (presumably large) to start computing values for
  :param deltaN: Provide an integer number to cap off the consecutive sums of prime
  :return: Returns an array of the normalized values for Zeta 
  """
  point = complex(0, 0)
  zeta_value = point

  # creating the interval range t inclusive [N, N+delta]
  t_range=np.arange(N, N+deltaN)

  # storage for the normal zeta 
  norm_zeta_value = np.array([0.0]*deltaN)

  # computing the Zeta function over a small interval defined by deltaN
  for t in t_range:
      point = complex(real=0.5, imag=t)
      zeta_value = mp.zeta(point)
      
      # computes the normalized zeta function -> refer to absoulte value code 
      norm_zeta_value[t-N] = abs(zeta_value)

  return norm_zeta_value