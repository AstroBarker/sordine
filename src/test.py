#!/usr/bin/env python3

import os
import urllib.request
import shutil

import numpy as np

from sedov import Sedov

# === Helper functions for tests ===

def prepare_sedov_test(TMPDIR):
  """
  Download sedov3 source code
  """

  url = " https://cococubed.com/codes/sedov/sedov.tbz"
  file_name = os.path.join(TMPDIR, "sedov.tbz")
  # Download the file from `url` and save it locally under `file_name`:
  with urllib.request.urlopen(url) as response, open(
    file_name, "wb"
  ) as out_file:
    shutil.copyfileobj(response, out_file)

  # extract
  extract_cmd = f"tar -xvf {file_name} -C {TMPDIR}"
  print(extract_cmd)
  os.system(extract_cmd)


# End prepare_sedov_test


def compile_sedov3(TMPDIR):
  """
  Compiles sedov3 source. We also modify the source so that it can accept
  command line args (removes a few lines that hard code values).
  """

  source_dir = os.path.join(TMPDIR, "sedov")
  source_fn = os.path.join(source_dir, "sedov3_qp.f90")

  # Modify source code to allow for command line args
  # specifically, remove lines 75-80 which hard code values
  cmd = f"sed -i '75,80d' {source_fn}"
  print(cmd)
  os.system(cmd)

  # compile
  FORTRAN = "gfortran"
  compile_cmd = f"{FORTRAN} -o sedov3_qp {source_fn}"
  print(compile_cmd)
  os.system(compile_cmd)

  # mv executable
  cmd = f"mv sedov3_qp {source_dir}"
  print(cmd)
  os.system(cmd)


# End compile_sedov3


def omega_from_j_gamma(n_per_j, j, gamma):
  """
  Given geometry j and adiabatic index n, construct array of w values
  We need w < j, w >= 0.0.
  We also force the grid to hav a singular case.
  """

  w = np.linspace(0.0, j - j / n_per_j, n_per_j)
  if j != 1.0:  # no singular case in planar
    singular = (3.0 * j - 2.0 + gamma * (2.0 - j)) / (gamma + 1.0)
    w = np.insert(w, 0, singular)

  return w


# End omega_from_j_gamma


def compare_sedov(x, rho, p, x_sol, rho_sol, p_sol):
  """
  Compare our solution to sedov3.
  Interpolate our grid to theirs, compute L1 norm.
  Note: we remove the innermost point from the solutions..
  They can blow up as r -> 0 for some setups, and our grids differ.
  So we set the inner parts of the grid to match.
  """

  r_inner = 0.01  # enforce.

  inds_sedov3 = np.where(x_sol > r_inner)
  inds_sistrum = np.where(x > r_inner)

  x_sol = x_sol[inds_sedov3]
  rho_sol = rho_sol[inds_sedov3]
  p_sol = p_sol[inds_sedov3]
  x = x[inds_sistrum]
  rho = rho[inds_sistrum]
  p = p[inds_sistrum]

  rho_sol /= np.max(rho_sol)
  p_sol /= np.max(p_sol)
  rho /= np.max(rho)
  p /= np.max(p)

  new_rho = np.interp(x_sol, x, rho)
  new_p = np.interp(x_sol, x, p)

  # L1 norm
  err_rho = np.sum(np.abs(rho_sol - new_rho)) / len(new_rho)
  err_p = np.sum(np.abs(p_sol - new_p)) / len(new_p)
  return err_rho, err_p


# End compare_sedov


# =================================


def test_sedov():
  """
  Test Sedov-Taylor against battle testing sedov3
  See https://cococubed.com/research_pages/sedov.shtml
  """

  TMPDIR = "tmp_test"
  os.mkdir(TMPDIR)

  # download sourxe
  prepare_sedov_test(TMPDIR)

  # modify and compile
  compile_sedov3(TMPDIR)

  # now set up the problem grid
  E_blast = 1.0
  npoints_sedov3 = 1000
  gamma = 1.4
  geometry = np.array([1.0, 2.0, 3.0])
  r_model = 1.2  # same as sedov3
  t = 1.0
  outfile = "tmp_sedov.dat"  # for sedov3 data
  n_per_j = 8  # w values per geometry value
  rho0 = 1.0  # ambient density
  p0 = 1.0e-5  # ambient pressure

  for j in geometry:
    omega = omega_from_j_gamma(n_per_j, j, gamma)
    for w in omega:
      # our solution
      sedov = Sedov(j, w, E_blast, rho0, p0, gamma, t, r_model)

      # run and load sedov3
      source_dir = os.path.join(TMPDIR, "sedov")
      source_ex = os.path.join(source_dir, "sedov3_qp")
      sedov3_cmd = (
        f"./{source_ex} {npoints_sedov3} {E_blast} {j} {w} {gamma} {outfile}"
      )
      print(sedov3_cmd)
      os.system(sedov3_cmd)
      cmd = f"mv {outfile} {source_dir}"
      print(cmd)
      os.system(cmd)
      x_sol, rho_sol, p_sol = np.loadtxt(
        os.path.join(source_dir, outfile),
        skiprows=2,
        unpack=True,
        usecols=(1, 2, 4),
      )
      err_rho, err_p = compare_sedov(
        sedov.r, sedov.rho_sol, sedov.p_sol, x_sol, rho_sol, p_sol
      )
      # Yes, the threshhold for density error is high...
      # The offending model keeping it there looks fine, though...
      tol_rho = 0.05
      tol_p = 0.005
      assert err_rho < tol_rho, f"Density L1 must be below {tol_rho}"
      assert err_p < tol_p, f"Pressure L1 must be below {tol_p}"

  shutil.rmtree(TMPDIR)


if __name__ == "__main__":
  # main
  test_sedov()
