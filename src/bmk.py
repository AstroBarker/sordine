#!/usr/bin/env python3
"""
 Purpose: Blandford McKee solution
 Author: Brandon L. Barker
"""

import numpy as np

class BMK:
  """
  Class for constructing Blandford-McKee (BMK) solution
  See https://ui.adsabs.harvard.edu/abs/1976PhFl...19.1130B/abstract

  Args:
    k (float) : Density power law index
    W (float) : Lorentz factor
    t (float) : time
    P0 (float) : Pressure scale
    t_shock (float) : time of shock formation (an offset for t)
    rad (float OR np.array) : positions to evaluate
  """

  def __init__(self, k, W, t, P0, t_shock, rad):
    self.k = k
    self.W = W
    self.t = t + t_shock
    self.P0 = P0

    self.m = 3 - k

    self.R_shock = self.R_()

    self.r = rad

    self.chi = self.chi_()
    self.g = self.g_()
    self.f = self.f_()
    self.h = self.h_()

    self.p = self.p_()

  # End __init__

  def R_(self):
    """
    Shock radius (Eq 26)
    """

    mp1 = self.m + 1.0
    W2 = self.W * self.W
    val = 1.0 - (1.0) / (2.0 * mp1 * W2)
    return self.t * val

  # End R_

  def chi_(self):
    """
    chi self similiar variable (Eq 27)
    """

    mp1 = self.m + 1.0
    term1 = 1.0 + 2.0 * mp1 * self.W * self.W
    term2 = 1.0 - self.r / self.t
    return term1 * term2

  # End chi_

  def g_(self):
    """
    Self similiar variable g (Eq 65)
    """

    return np.power(self.chi, -1)

  # End g_

  def f_(self):
    """
    Self similar variable f (Eq 66)
    """
    exponent = (4.0 * self.k - 17.0) / (12.0 - 3 * self.k)
    return np.power(self.chi, exponent)

  # End f_

  def h_(self):
    """
    Self similiar variable h (Eq 67)
    """

    exponent = (2.0 * self.k - 7) / (4.0 - self.k)
    return np.power(self.chi, exponent)

  # End h_

  def gamma_square_(self):
    """
    Self similiar function gamma^2 (Eq 29)
    """
    return 0.5 * self.W * self.W * self.g

  # End gamma_square_

  def p_(self):
    """
    Self similiar pressure function (Eq 30)
    """
    return self.P0 * self.f

  # End p_


# End BMK

if __name__ == "__main__":
  # example:

  t = 0.5
  t_shock = 0.25
  W = 0.5
  k = 0.0
  P0 = 1.0e-4
  r = np.linspace(0.0, 1.0, 1000)
  bmk = BMK(k, W, t, P0, t_shock, r)

# End main
