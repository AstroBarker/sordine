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

  url =" https://cococubed.com/codes/sedov/sedov.tbz"
  file_name = os.path.join(TMPDIR, "sedov.tbz")
  # Download the file from `url` and save it locally under `file_name`:
  with urllib.request.urlopen(url) as response, open(file_name, 'wb') as out_file:
    shutil.copyfileobj(response, out_file)

  # extract
  source_dir = os.path.join(TMPDIR, "sedov3")
  os.mkdir(source_dir)
  extract_cmd = f"tar -xvf {file_name} -C {source_dir}" 

def test_sedov():
  """
  Test Sedov-Taylor against battle testing sedov3
  See https://cococubed.com/research_pages/sedov.shtml
  """

  TMPDIR = "tmp_test"
  os.mkdir(TMPDIR)

  prepare_sedov_test(TMPDIR)

  #shutil.rmtree(TMPDIR)

if __name__ == "__main__":
  # example:

  test_sedov()

