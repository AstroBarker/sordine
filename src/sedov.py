#!/usr/bin/env python3
"""
 Purpose: Sedov solution
 Based on work by James R. Kamm (https://cococubed.com/papers/kamm_2000.pdf)
 Author: Brandon L. Barker
"""

from enum import Enum

import numpy as np
import scipy.integrate as integrate  # quad
from scipy import optimize  # brentq, root


class Family(Enum):
  """
  Enum class to act as type of solution
  Names:
    standard
    singular
    vacuum
  """

  standard = 0
  singular = 1
  vacuum = 2


# End Family


class Singularity(Enum):
  """
  Enum class to hold singularity type
  Names:
    none
    omega2
    omega3
  """

  none = 0
  omega2 = 1
  omega3 = 2


# End Singularity


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
    p0 (float) : Ambient pressure
    gamma (float) : adiabatic index
    t_end (float) : end time
    r_model (float) : radial extent

  Attributes:
    self.r (np.array) : radial grid
    self.rho_sol (np.array) holds solution density profile

  Usage:
    >>> sedov = Sedov(j, w, Eexp, rho0, p0, gamma, t_end, r_model)
    >>> rho_solution = sedov.rho_sol

  TODO:
    - Add velocity solution
  """

  def __init__(self, j, w, E, rho0, p0, gamma, t_end, r_model):
    # Various checks for a physical solution or realistic
    assert (
      j == 1 or j == 2 or j == 3
    ), "Please select an appropriate geometry (j=1,2,3)"
    assert E > 0.0, "Explosion energy must be positive"
    assert w >= 0.0, "Density power law index w must be non-negative"
    assert w < j, "Must have w < j for finite mass."

    self.j = j
    self.w = w
    self.E = E
    self.rho0 = rho0
    self.p0 = p0
    self.gamma = gamma
    self.t = t_end
    self.r_model = r_model

    # set up grid
    npoints_r = 2048
    r_min = 1.0e-5  # avoid r = 0
    self.r = np.linspace(r_min, r_model, npoints_r)
    self.rho_sol = np.zeros(npoints_r)  # density solution
    self.p_sol = np.zeros(npoints_r)  # pressure solution

    # other stuff
    j2w = j + 2.0 - w
    self.j2w = j2w

    # First, check solution family
    self.V2 = 4.0 / (j2w * (gamma + 1.0))  # after equation (17)
    self.Vstar = 2.0 / (self.j * (gamma - 1.0) + 2.0)  # Eq 19
    self.V0 = 2.0 / (gamma * j2w)  # equation (23)

    # Following Kamm & Timmes, introduce TOL in logic for determining solution family
    # limits exponents later to not blow up
    TOL = 1.0e-5
    self.family = Family.standard
    if self.V2 < self.Vstar - TOL:
      self.family = Family.standard
    elif self.V2 > self.Vstar + TOL:
      self.family = Family.vacuum
    elif abs(self.V2 - self.Vstar) < TOL:
      self.family = Family.singular
    else:
      raise ValueError(
        "Something weird has happened: initial condition does not correspond to any solution family"
      )

    # Check for removable singularities
    self.singularity = Singularity.none
    #w1 = (3.0 * j - 2.0 + gamma * (2.0 - j)) / (gamma + 1.0)
    self.w2 = (2.0 * (gamma - 1.0) + j) / gamma
    self.w3 = j * (2.0 - gamma)
    if abs(w - self.w2) < TOL:
      self.singularity = Singularity.omega2
    if abs(w - self.w3) < TOL:
      self.singularity = Singularity.omega3

    # Equations (33)-(37)
    self.a = j2w * (gamma + 1.0) / 4.0
    self.b = (gamma + 1.0) / (gamma - 1.0)
    self.c = j2w * gamma / 2.0
    if self.family == Family.singular:
      self.d = 1.0  # unused in singular case but avoids div by 0
    else:
      self.d = (j2w * (gamma + 1.0)) / (
        j2w * (gamma + 1.0) - 2.0 * (2.0 + j * (gamma - 1.0))
      )
    self.e = (2.0 + j * (gamma - 1.0)) / 2.0

    # Equations (42)-(47)
    # The if logic below avoids divide by zero in the omega2, omega3 singularity cases.
    # We avoid these singularities so no harm is done by modifying alpha2, alpha4, alpha5
    denom2 = (
      gamma * (self.w2 - w) if self.singularity != Singularity.omega2 else TOL
    )
    denom3 = (self.w3 - w) if self.singularity != Singularity.omega3 else TOL
    self.alpha0 = 2.0 / j2w
    self.alpha2 = -(gamma - 1.0) / denom2
    self.alpha1 = (
      (gamma * j2w)
      / (2.0 + j * (gamma - 1.0))
      * ((2.0 * (j * (2.0 - gamma) - w)) / (gamma * j2w**2) - self.alpha2)
    )
    self.alpha3 = (j - w) / denom2
    self.alpha4 = self.alpha1 * (j - w) * j2w / denom3
    self.alpha5 = (w * (gamma + 1.0) - 2.0 * j) / denom3

    # perform energy entegral, stores alpha, J1, J2
    self.energy_integral_()

    # shock position
    self.r_sh = self.r_shock_()

    # deal with vacuum radius
    self.r_vac = 0.0
    if self.family == Family.vacuum:
      V_vac = 2.0 / self.j2w
      self.r_vac = optimize.brentq(
        self.target_v_, self.r[0], self.r_sh, args=(V_vac), xtol=1.0e-20
      )

    # Do the solve
    self.solve_()

  # End __init__

  def __str__(self):
    print("Sedov problem parameters:")
    print(f"Geometric index   : {self.j}")
    print(f"Explosion energy  : {self.E}")
    print(f"rho0              : {self.rho0}")
    print(f"p0                : {self.p0}")
    print(f"Density power law : {self.w}")
    print(f"Adiabatic index   : {self.gamma}")
    print(f"t                 : {self.t}")
    print(f"r_model           : {self.r_model}")
    print(f"Solution family   : {self.family}")
    print(f"Singularity       : {self.singularity}")
    return ""

  # End __str__

  def r_shock_(self):
    """
    Return shock radius (Eq 14)
    """
    return (self.t * self.t * self.E / (self.rho0 * self.alpha)) ** (
      1.0 / self.j2w
    )  # equation (14)

  def Integrand1_(self, V):  # equation (55) of K&T
    """
    Energy integrand 1 (Equation 73)
    """
    return (
      ((self.gamma + 1.0) / (self.gamma - 1.0))
      * self.lambda_(V) ** (self.j + 1.0)
      * self.g_(V)
      * self.dlam_dV_(V)
      * V
      * V
    )

  # End Integrand1_

  def Integrand2_(self, V):  # equation (56) of K&T
    """
    Energy integrand 2 (Equation 74)
    """

    return (
      (8.0 / ((self.gamma + 1.0) * self.j2w * self.j2w))
      * self.lambda_(V) ** (self.j - 1.0)
      * self.h_(V)
      * self.dlam_dV_(V)
    )

  # End Integrand2

  def energy_integral_(self):
    """
    return J1 and J2 energy integrals
    # the ennergy integrals can contain singularities at the lower bound.
    # for now, offset by machine epsilon to avoid.
    # TODO: Implement singularity removal ala https://cococubed.com/papers/la-ur-07-2849.pdf
    # TODO: set Vmin appropriately
    """
    self.J1 = 0.0
    self.J2 = 0.0
    self.alpha = 1.0

    if self.family == Family.singular:  # singular case, integration is analytic
      self.J1 = (self.gamma + 1.0) / (
        self.j * ((self.gamma - 1.0) * self.j + 2.0) ** 2.0
      )
      self.J2 = 2.0 * self.J1 / (self.gamma - 1.0)
      self.alpha = (
        ((self.gamma + 1.0) / (self.gamma - 1.0))
        * 2.0 ** (self.j)
        / (self.j * ((self.gamma - 1.0) * self.j + 2.0) ** 2.0)
      )
      if self.j != 1.0:
        self.alpha *= np.pi

    else:  # standard and vacuum, integrate numerically
      eps = self.V0 * np.finfo(float).eps
      vmin = 2.0 / self.j2w if self.family == Family.vacuum else self.V0
      self.J1, err1 = integrate.quad(
        self.Integrand1_, vmin + eps, self.V2, epsabs=eps
      )
      self.J2, err2 = integrate.quad(
        self.Integrand2_, vmin + eps, self.V2, epsabs=eps
      )

      if self.family == Family.singular:
        self.alpha = np.pi * self.J2 * 2.0 ** (self.j - 1.0)
      else:
        if self.j == 1:
          self.alpha = self.J1 + 2.0 * self.J2 / (self.gamma - 1.0)
        else:
          self.alpha = (
            (self.j - 1.0)
            * np.pi
            * (self.J1 + 2.0 * self.J2 / (self.gamma - 1.0))
          )

  # End energy_integral_

  def target_r_(self, V, rx):
    """
    Root find r = r_sh lambda(V) for V
    """
    # fac accounts for a factor of r in lambda for the singular case
    fac = 1.0 if self.family != Family.singular else rx
    return self.r_sh * fac * self.lambda_(V) - rx

  # End target_r_

  def target_v_(self, rx, V_vac):
    """
    Solving for vacuum radius
    Root find r_vac = r_sh lambda(V_vac) for r_vac
    """
    # fac accounts for a factor of r in lambda for the singular case
    fac = 1.0 if self.family != Family.singular else rx
    return self.r_sh * fac * self.lambda_(V_vac) - rx

  # End target_v_

  def lambda_(self, V):
    """
    self similar function lambda
    depends on solution family
    See Kamm & Timmes section 2
    NOTE: The singular case is missing a factor of radius, because it's a pain to
    pull through all of these functions. It's accounted for in other odd places.
    """
    eps = 1.0e-60
    x1 = self.a * V
    x2 = self.b * max(eps, self.c * V - 1.0)
    x3 = self.d * (1.0 - self.e * V)

    val = 0.0
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity != Singularity.omega2:
      val = x1 ** (-self.alpha0) * x2 ** (-self.alpha2) * x3 ** (-self.alpha1)
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega2:
      val = (
        x1 ** (-self.alpha0)
        * x2 ** ((self.gamma - 1.0) / (2.0 * self.e))
        * np.exp(
          (self.gamma - 1.0)
          * (1.0 - x1)
          / ((2.0 * self.e) * (x1 - (self.gamma + 1) / (2.0 * self.gamma)))
        )
      )
    if self.family == Family.singular:
      val = 1.0 / self.r_sh  # scale by r elsewhere
    return val

  # End lambda_

  def dlam_dV_(self, V):
    """
    Compute d lambda / dV for energy integral
    See Kamm & Timmes section 2
    """
    eps = 1.0e-60
    x1 = self.a * V
    x2 = self.b * max(eps, self.c * V - 1.0) + eps
    x3 = self.d * (1.0 - self.e * V)
    x4 = self.b * (1.0 - self.c * V / self.gamma)
    dx1dv = self.a
    dx2dv = self.b * self.c
    dx3dv = -self.d * self.e
    dx4dv = -self.b * self.c / self.gamma
    lam = self.lambda_(V)
    dlamdv = 0.0
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.none:
      dlamdv = -(
        (self.alpha0 * dx1dv / x1)
        + (self.alpha2 * dx2dv / x2)
        + (self.alpha1 * dx3dv / x3)
      )
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega2:
      term1 = self.alpha0 * dx1dv / x1
      term2 = (self.gamma - 1.0) * dx2dv / (2.0 * self.e * x2)
      term3 = -(self.gamma + 1.0) * dx1dv / (2.0 * self.e)
      term4 = 1.0 / (x1 - (self.gamma + 1.0) / (2.0 * self.gamma))
      term5 = 1.0 + (1.0 - x1) / (x1 - (self.gamma + 1.0) / (2.0 * self.gamma))
      dlamdv = -(term1 + term2 + term3 * term4 * term5)
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega3:
      dlamdv = -(
        (self.alpha0 * dx1dv / x1)
        + (self.alpha2 * dx2dv / x2)
        + (self.alpha1 * dx4dv / x4)
      )
    if self.family == Family.singular:
      dlamdv = 0.0
    return lam * dlamdv

  # End dlam_dV_

  def f_(self, V):
    """
    self similar function f
    depends on solution family
    See Kamm & Timmes section 2
    """
    x1 = self.a * V
    val = 0.0
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity != Singularity.omega2:
      val = x1 * self.lambda_(V)
    if self.family == Family.singular:
      val = self.lambda_(V)
    return val

  # End f_

  def g_(self, V):
    """
    self similar function g
    depends on solution family
    See Kamm & Timmes section 2
    """
    eps = 1.0e-60
    x1 = self.a * V
    x2 = self.b * max(eps, self.c * V - 1.0)
    x3 = self.d * (1.0 - self.e * V)
    x4 = self.b * (1.0 - self.c * V / self.gamma)

    val = 0.0
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.none:
      val = (
        x1 ** (self.alpha0 * self.w)
        * x2 ** (self.alpha3 + self.alpha2 * self.w)  # round off protection..
        * x3 ** (self.alpha4 + self.alpha1 * self.w)
        * x4 ** (self.alpha5)
      )
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega2:
      val = (
        x1 ** (self.alpha0 * self.w)
        * x2 ** (4.0 - self.j - 2.0 * self.gamma / (2.0 * self.e))
        * x4 ** (self.alpha5)
        * np.exp(
          (self.gamma + 1.0)
          * (1.0 - x1)
          / (self.e * (x1 - (self.gamma + 1) / (2.0 * self.gamma)))
        )
      )
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega3:
      val = (
        x1 ** (self.alpha0 * self.w)
        * x2 ** (self.alpha3 + self.alpha2 * self.w)
        * x4 ** (1.0 - 2.0 / self.e)
        * np.exp(
          -self.j
          * self.gamma
          * (self.gamma + 1)
          * (1.0 - x1)
          / ((2.0 * self.e) * ((self.gamma + 1.0) / 2.0 - x1))
        )
      )
    if self.family == Family.singular:
      val = self.lambda_(V) ** (self.j - 2.0)
    return val

  # End g_

  def h_(self, V):
    """
    self similar function h
    depends on solution family
    See Kamm & Timmes section 2
    """
    x1 = self.a * V
    x3 = self.d * (1.0 - self.e * V)
    x4 = self.b * (1.0 - self.c * V / self.gamma)

    val = 0.0
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.none:
      val = (
        x1 ** (self.alpha0 * self.j)
        * x3 ** (self.alpha4 + self.alpha1 * (self.w - 2.0))
        * x4 ** (1.0 + self.alpha5)
      )
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega2:
      val = (
        x1 ** (self.alpha0 * self.j)
        * x3 ** (-self.j * self.gamma / (2.0 * self.e))
        * x4 ** (1.0 + self.alpha5)
      )
    if (
      self.family == Family.standard or self.family == Family.vacuum
    ) and self.singularity == Singularity.omega3:
      val = (
        x1 ** (self.alpha0 * self.j)
        * x4 ** ((self.j * (self.gamma - 1.0) - self.gamma) / self.e)
        * np.exp(
          -self.j
          * self.gamma
          * (self.gamma + 1)
          * (1.0 - x1)
          / ((2.0 * self.e) * ((self.gamma + 1.0) / 2.0 - x1))
        )
      )
    if self.family == Family.singular:
      val = self.lambda_(V) ** (self.j)
    return val

  # End h_

  def rho_(self, V, rho2):  # equation (40)
    """
    Returns post shock density (Eq 40)
    """
    return rho2 * self.g_(V)

  # End rho_

  def p_(self, V, p2):
    """
    Calculate post shock pressure (Eq 41)
    """
    return p2 * self.h_(V)

  # End p_

  def rho2(self, rho1):
    """
    Eq 13
    """
    return rho1 * (self.gamma + 1.0) / (self.gamma - 1.0)  # equation (13)

  def solve_(self):
    """
    Solve main state for standard case
    """
    rho1 = self.rho0 * self.r_sh ** (-self.w)
    rho2 = self.rho2(rho1)  # equation (13)
    v_sh = 2.0 * self.r_sh / (self.j2w * self.t)
    p2 = 2.0 * rho1 * v_sh * v_sh / (self.gamma + 1.0)

    for i in range(len(self.r)):
      r = self.r[i]
      if self.family == Family.vacuum and r < self.r_vac:
        continue
      if r >= 0.0 and r < self.r_sh:  # shocked region
        V_x = self.Vstar  # singular case
        if self.family != Family.singular:
          vmin = 0.9 * self.V0 if self.family == Family.standard else self.V2
          vmax = (
            1.0 * self.V2 if self.family == Family.standard else 2.0 / self.j2w
          )
          V_x = optimize.brenth(
            self.target_r_, vmin, vmax, args=(r), xtol=1.0e-20
          )

        # update post-shock state
        # fac adjusts for missing r in singular lambda
        fac = r**(self.j - 2.0) if self.family == Family.singular else 1.0
        self.rho_sol[i] = fac * self.rho_(V_x, rho2)
        fac = r**(self.j) if self.family == Family.singular else 1.0
        self.p_sol[i] = fac * self.p_(V_x, p2)
      else:  # unshocked region
        self.rho_sol[i] = self.rho0 * r ** (-self.w)
        self.p_sol[i] = self.p0

  # End solve_


# End Sedov

if __name__ == "__main__":
  # example:
  j = 2
  w = 1.99
  E = 1.0
  rho0 = 1.0
  p0 = 1.0e-5
  gamma = 1.4
  t_end = 1.0
  r_model = 1.2
  sedov = Sedov(j, w, E, rho0, p0, gamma, t_end, r_model)
  print(sedov)

# End main
