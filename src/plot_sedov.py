#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from sedov import Sedov, Family

def plot_density():
  """
  hacky plot density stuff.
  """
  w = np.linspace(1.0, 2.75, 24)
  w = np.insert(w, 0, 0.0)
  w = np.insert(w, 0, 2.0)


  j = 3.0 # spherical
  gamma = 5.0 / 3.0
  t_end = 1.0
  r_model = 1.8
  Eexp = 5.45670 #0.851072
  rho0 = 1.0
  P0 = 1.0e-5


  fig, ax = plt.subplots()
  for i in range(len(w)):
    print(w[i])
    sedov = Sedov(j, w[i], Eexp, rho0, P0, gamma, t_end, r_model)
    rho1 = rho0 * sedov.r_sh**(-w[i])

    # options
    color = "teal"
    style = "-"
    label = ""
    if sedov.family == Family.standard:
      color = "teal"
      label = r"Standard, $\omega <$ j/$\gamma$"
      style = "--"
      if w[i] >= j / gamma:
        color = "slateblue"
        label = r"Standard, $\omega \geq$ j/$\gamma$"
    elif sedov.family == Family.singular:
      color = "#808000"
      style = "-"
      label = r"Singular, $\omega = (7-\gamma)/(\gamma + 1)$"
    else:
      color = "#800080"
      style = "-."
      label = r"Vacuum, $\omega \leq 2j/(\gamma + 1)$"
      if w[i] > 2.0 * j /(gamma + 1.0):
        color = "orchid"
        style = "-."
        label = r"Vacuum, $\omega > 2j/(\gamma + 1)$"

    ax.plot( sedov.r / sedov.r_sh, sedov.rho_sol / sedov.rho2(rho1), color = color, ls = style, lw = 0.75, label=label)

    #ax.plot(x_s3[1:], rho_s3[1:], color="k", ls=" ", marker=".")
    #ax.plot( sedov.r, sedov.rho_sol, color = "teal")


  ax.set(xlim = [0.0, 1.0], ylim = [0.0, 1.0], xlabel = "position", ylabel = "Scaled Density")
  handles, labels = plt.gca().get_legend_handles_labels()
  by_label = dict(zip(labels, handles))
  ax.legend(by_label.values(), by_label.keys(), prop={'size': 8})
  ax.set(xlim = [0.0, 1.0], ylim = [0.0, 1.0], xlabel = "position", ylabel = "Scaled Density")
  plt.savefig("sedov.png")
    
plot_density()
