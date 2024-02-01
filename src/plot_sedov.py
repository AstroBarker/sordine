#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from sedov import Sedov

def plot_density():
  """
  plot density
  """
  w = np.linspace(1.5, 1.99, 10)
  w = np.insert(w, 0, 0.0)

  j = 3 # spherical
  gamma = 5.0 / 3.0
  t_end = 0.25
  r_model = 1.0
  Eexp = 1.0
  rho0 = 1.0

  print(w)

  fig, ax = plt.subplots()
  for i in range(len(w)):
    rho = rho0**(-w[i])
    sedov = Sedov(j, w[i], Eexp, rho, gamma, t_end, r_model)
    ax.plot( sedov.r / sedov.r_sh, sedov.rho_sol / sedov.rho2(rho0*sedov.r_sh**(-w[i])), color = "teal")
  ax.set(xlim = [0.0, 1.0], ylim = [0.0, 1.0], xlabel = "position", ylabel = "Scaled Density")
  plt.savefig("sedov.png")
    
plot_density()
