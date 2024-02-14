#!/usr/bin/env python3

import os
import urllib.request
import shutil

import numpy as np
import matplotlib.pyplot as plt

from sedov import Sedov, Family


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
  os.system(cmd)


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

  shutil.rmtree(TMPDIR)


if __name__ == "__main__":
  # example:

  test_sedov()
