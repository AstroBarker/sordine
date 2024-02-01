#!/usr/bin/env python3
"""
 Purpose: Sedov solution
 Based on work by James R. Kamm (https://cococubed.com/papers/kamm_2000.pdf)
 Author: Brandon L. Barker
"""


import numpy as np
import scipy.integrate as integrate  # quad
from scipy import optimize  # brentq


class Sedov:
  """
  Class for constructing Sedov solution
  See https://cococubed.com/papers/kamm_2000.pdf

  Assumes density profile of form rho(r) = rho0 * r^(-w)

  Args:
    j (float) : Dimension index (1: planar, 2: cylindrical. 3: spherical)
    w (float) : Density power law index
    E (float) : Explosion energy
    rho0 (float) : Ambient density
    gamma (float) : adiabatic index
    t_end (float) : end time
    r_model (float) : radial extent

  Attributes:
    self.r (np.array) : radial grid
    self.rho_sol (np.array) holds solution density profile

  Usage:
    >>> sedov = Sedov(j, w, Eexp, rho0, gamma, t_end, r_model)
    >>> rho_solution = sedov.rho_sol

  TODO:
    - Add more quantities of interest.
  """

  def __init__(self, j, w, E, rho0, gamma, t_end, r_model):
    assert (
      j == 1 or j == 2 or j == 3
    ), "Please select an appropriate geometry (j=1,2,3)"
    assert E > 0.0, "Explosion energy must be positive"

    self.j = j
    self.w = w
    self.E = E
    self.rho0 = rho0
    self.gamma = gamma
    self.t_end = t_end
    self.r_model = r_model

    # set up grid
    npoints_r = 32768
    self.r = np.linspace(0.0, r_model, npoints_r)
    self.rho_sol = np.zeros(npoints_r)
    self.t = t_end

    # other stuff
    j2w = j + 2.0 - w
    self.j2w = j2w

    # Equations (33)-(37)
    self.a = j2w * (gamma + 1.0) / 4.0
    self.b = (gamma + 1.0) / (gamma - 1.0)
    self.c = j2w * gamma / 2.0
    self.d = (j2w * (gamma + 1.0)) / (
      j2w * (gamma + 1.0) - 2.0 * (2.0 + j * (gamma - 1.0))
    )
    self.e = (2.0 + j * (gamma - 1.0)) / 2.0

    # Equations (42)-(47)
    self.alpha0 = 2.0 / j2w
    self.alpha2 = -(gamma - 1.0) / (2.0 * (gamma - 1.0) + j - gamma * w)
    self.alpha1 = (
      (gamma * j2w)
      / (2.0 + j * (gamma - 1.0))
      * ((2.0 * (j * (2.0 - gamma) - w)) / (gamma * j2w**2) - self.alpha2)
    )
    self.alpha3 = (j - w) / (2.0 * (gamma - 1.0) + j - gamma * w)
    self.alpha4 = self.alpha1 * (j - w) * j2w / (j * (2.0 - gamma) - w)
    self.alpha5 = (w * (gamma + 1.0) - 2.0 * j) / (j * (2.0 - gamma) - w)

    self.V2 = 4.0 / (j2w * (gamma + 1.0))  # after equation (17)
    self.V0 = 2.0 / (gamma * j2w)  # equation (23)

    eps = np.finfo(float).eps
    result1, err1 = integrate.quad(self.Integrand1_, self.V0 + eps, self.V2)
    result2, err2 = integrate.quad(self.Integrand2_, self.V0 + eps, self.V2)

    self.J1 = result1  # equation (67)
    self.J2 = result2  # equation (68)

    factor = 1.0 if (j == 1) else np.pi
    self.alpha = (
      2.0 ** (j - 2.0) * factor * result1
      + (2.0 ** (j - 1.0)) / (gamma - 1.0) * factor * result2
    )  # equation (66)

    self.r_sh = self.r_shock_()

    # Do the solve
    self.Solve_()

  # End __init__

  def __str__(self):
    print("Sedov problem parameters:")
    print(f"Geometric index   : {self.j}")
    print(f"Explosion energy  : {self.E}")
    print(f"rho0              : {self.rho0}")
    print(f"Density power law : {self.w}")
    print(f"Adiabatic index   : {self.gamma}")
    print(f"t_end             : {self.t_end}")
    print(f"r_model           : {self.r_model}")
    return ""

  # End __str__

  def r_shock_(self):
    """
    Return shock radius (Eq 14)
    """
    return (self.E / (self.rho0 * self.alpha)) ** (1.0 / self.j2w) * self.t ** (
      2.0 / self.j2w
    )  # equation (14)

  def Integrand1_(self, V):  # equation (73)
    """
    Energy integrand 1 (Equation 73)
    """
    x1 = self.a * V
    x2 = self.b * (self.c * V - 1)
    x3 = self.d * (1 - self.e * V)
    x4 = self.b * (1 - self.c * V / self.gamma)
    return (
      -(self.gamma + 1.0)
      / (self.gamma - 1.0)
      * V**2
      * (
        self.alpha0 / V
        + self.alpha2 * self.c / (self.c * V - 1.0)
        - self.alpha1 * self.e / (1.0 - self.e * V)
      )
      * ((x1**self.alpha0) * (x2**self.alpha2) * (x3**self.alpha1))
      ** (-self.j2w)
      * x2**self.alpha3
      * x3**self.alpha4
      * x4**self.alpha5
    )

  # End Integrand1_

  def Integrand2_(self, V):  # equation (74)
    """
    Energy integrand 2 (Equation 74)
    """
    x1 = self.a * V
    x2 = self.b * (self.c * V - 1.0)
    x3 = self.d * (1.0 - self.e * V)
    x4 = self.b * (1.0 - self.c * V / self.gamma)
    return (
      -(self.gamma + 1.0)
      / (2.0 * self.gamma)
      * V**2
      * (self.c * V - self.gamma)
      / (1.0 - self.c * V)
      * (
        self.alpha0 / V
        + self.alpha2 * self.c / (self.c * V - 1.0)
        - self.alpha1 * self.e / (1.0 - self.e * V)
      )
      * (x1**self.alpha0 * x2**self.alpha2 * x3**self.alpha1) ** (-self.j2w)
      * x2**self.alpha3
      * x3**self.alpha4
      * x4**self.alpha5
    )

  # End Integrand2

  def target_r_(self, V, rx):
    """
    Root find r = r_sh lambda(V) for V
    """
    x1 = self.a * V
    x2 = self.b * (self.c * V - 1.0)
    x3 = self.d * (1.0 - self.e * V)
    return (
      self.r_sh
      * x1 ** (-self.alpha0)
      * x2 ** (-self.alpha2)
      * x3 ** (-self.alpha1)
      - rx
    )

  # End target_r_

  def rho_(self, V, rho2):  # equation (40)
    """
    Returns density
    """
    x1 = self.a * V
    x2 = self.b * (self.c * V - 1.0)
    x3 = self.d * (1.0 - self.e * V)
    x4 = self.b * (1.0 - self.c * V / self.gamma)
    return (
      rho2
      * x1 ** (self.alpha0 * self.w)
      * x2 ** (self.alpha3 + self.alpha2 * self.w)
      * x3 ** (self.alpha4 + self.alpha1 * self.w)
      * x4**self.alpha5
    )

  # End rho_

  def rho2(self, rho1):
    """
    Eq 13
    """
    return rho1 * (self.gamma + 1.0) / (self.gamma - 1.0)  # equation (13)

  def Solve_(self):
    """
    Solve main state
    """
    rho1 = self.rho0 * self.r_sh ** (-self.w)
    rho2 = self.rho2(rho1)  # equation (13)
    CUTOFF_ = 0.001  # below 0.3 * r_shock, set density to 0.0
    for i in range(len(self.r)):
      r = self.r[i]
      if r < CUTOFF_ * self.r_sh:
        self.rho_sol[i] = 0.0
      elif r >= CUTOFF_ * self.r_sh and r < self.r_sh:  # do work here
        V_x = optimize.brentq(self.target_r_, self.V0, self.V2, args=(r))
        self.rho_sol[i] = self.rho_(V_x, rho2)
      else:  # unshocked regio
        self.rho_sol[i] = rho1

  # End Solve_


# End Sedov

if __name__ == "__main__":
  # example:
  j = 1
  w = 0.0
  E = 1.0
  rho0 = 1.0
  gamma = 1.4
  t_end = 0.5
  r_model = 1.0
  sedov = Sedov(j, w, E, rho0, gamma, t_end, r_model)
  print(sedov)

# End main
