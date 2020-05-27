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


def zetaStochastic(h : float, T: int, prime_list: Union[list, np.array] = None) -> float:
  """
  An itterative algorithm for computing the stochastic Zeta function defined X_t(h)
  :param h: Provide a floating interval range to observe the Zeta function
  :param T: Provide an integer number to cap off the consecutive sums of prime
  :return: Returns a value for the stochastic Zeta function 
  """
  assert T < int(1e7), 'Our prime list does not exceed the value 10^10'
  ret_val = 0.0
  p = 0.0

  for i in range(len(prime_list)): 
    p = prime_list[i]
    
    # add the sequential terms of the primes up to some integer T 
    if (p <= T):
      ret_val += 1/np.sqrt(p) * np.cos(np.random.uniform(0, 2*np.pi) - h*np.log(p))
    else:
      break

  return ret_val