# SOftwaRe valiDatIonN suitE (SORDINE)
A musical mute and also a repository of Python 3.x [self-similar solutions](https://en.wikipedia.org/wiki/Self-similar_solution) that I've needed, and so written.
Contributions welcome.
A minimally working [environment file](environment.yml) is provided.

# [Sedov-Taylor blast wave](src/sedov.py)
Self-similar solution for the Sedov-Taylor non-relativistic blast wave.
Based on work by [James R. Kamm](https://cococubed.com/papers/kamm_2000.pdf) and [Kamm & Timmes](https://cococubed.com/papers/la-ur-07-2849.pdf).
Supports density profiles of the from $\rho(r) = \rho_{0} r^{-w}$.
Supports planar (`j = 1`), cylindrical (`j = 2`), and spherical (`j = 3`) geometries.
There are constraints on the density power law index. 
See the paper.

## Usage
```python
from sedov import Sedov
sedov_solution = Sedov(j, w, E, rho0, gamma, t_end, r_model)
r = sedov.r
rho_sol = sedov.rho_sol
# plot away
```

## Testing
We include `test.py` that implements tests for the Sedov-Taylor solutions.
For this test we set up a grid in geometry, density power law space.
To test, we download and compile the `sedov3` [code](https://cococubed.com/research_pages/sedov.shtml) 
provided by Frank Timmes.
We compare our solution to his.
This is triggered on every push / PR, and you can invoke by, e.g., `pytest src/test.py`.
In the future it should test other solutions.

## TODO:
- Some solutions, such as velocity, still need implementation.

![Self similar density profiles](sedov.png "Sedov density profiles")

# [Blandford McKee Blast Wave](src/bmk.py)
See [the paper](https://ui.adsabs.harvard.edu/abs/1976PhFl...19.1130B/abstract).
Computes post-shock similarity variables and pressure, etc.
Assumes density profiles of the from $\rho(r) \propto r^{-k}$.

## Usage
```python
from bmk import BMK
bmk_solution = BMK(k, W, t, P0, t_shock, rad)
rad = bmk_solution.r
pressure = bmk_solution.p # see also bmk_solution.f for similarity variable
```

# Code Style
Code linting and formatting is done with [ruff](https://docs.astral.sh/ruff/).
Rules are listed in [ruff.toml](ruff.toml).
To check all python in the current directory, you may `ruff .`.
To format a given file according to `ruff.toml`, run `ruff format file.py`.
Checks for formatting are performed on each push / PR.
We (try to) adhere to [semantic versioning](https://semver.org/) and [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/).
