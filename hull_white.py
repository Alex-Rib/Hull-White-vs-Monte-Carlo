"""Hull-White lattice pricer for arithmetic-average Asian call options

Implements the Hull-White (1993) approach on a recombining binomial tree: at
each node a grid of representative running averages is carried, bounded by the
minimum and maximum attainable average (forward induction), and the option
value is rolled back by backward induction with linear interpolation between
grid averages. Reference: Van der Hoek & Elliott, Binomial Models in Finance (2006)
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HullWhiteResult:
    """Outcome of a Hull-White pricing run"""

    price: float
    n_averages: int
    elapsed: float


def _interpolate(avg_grid: np.ndarray, target: float, child_values: np.ndarray) -> float:
    """Linear interpolation of the child option value at a target average

    Given the grid of averages and the matching option values at a child node,
    finds the bracketing interval for target and interpolates linearly. Falls
    back to the bracket value when the interval collapses
    """
    valid = ~np.isnan(avg_grid) & ~np.isnan(child_values)
    grid = avg_grid[valid]
    children = child_values[valid]

    order = np.argsort(grid)
    grid = grid[order]
    children = children[order]

    tolerance = np.sqrt(np.finfo(float).eps)
    idx = np.searchsorted(grid, target, side="right") - 1
    idx = max(0, min(idx, len(grid) - 2))

    width = grid[idx + 1] - grid[idx]
    if abs(width) < tolerance:
        return float(children[idx])

    weight = (grid[idx + 1] - target) / width
    return float(weight * children[idx] + (1.0 - weight) * children[idx + 1])


class HullWhiteAsianPricer:
    """Hull-White lattice pricer for an arithmetic-average Asian call

    n_averages controls the number of representative averages carried per node.
    With n_averages = 4 the grid uses the book's A_min, two interior thirds and
    A_max; otherwise the averages are spread uniformly between the bounds
    """

    def __init__(self, params, *, n_averages: int = 4) -> None:
        """Store the market parameters and the averaging-grid resolution"""
        self.params = params
        self.n_averages = n_averages

        p = params
        self.dt = p.maturity / p.n_steps
        self.up = np.exp(p.sigma * np.sqrt(self.dt))
        self.down = 1.0 / self.up
        self.growth = np.exp(p.rate * self.dt)
        self.q = (self.growth - self.down) / (self.up - self.down)

    def _average_bounds(self, i: int, j: int) -> tuple[float, float]:
        """Smallest and largest running average attainable at node (i, j)

        Computed by forward induction over the i steps, where j is the number
        of down moves so far
        """
        p = self.params
        s0, h, b = p.spot, self.up, self.down

        # minimum average: down moves first, then up moves
        term1 = s0 * j if abs(b - 1) < 1e-10 else s0 * b * (1 - b**j) / (1 - b)
        term2 = s0 * (i - j) if abs(h - 1) < 1e-10 else s0 * b**j * h * (1 - h ** (i - j)) / (1 - h)
        a_min = (s0 + term1 + term2) / (i + 1)

        # maximum average: up moves first, then down moves
        term1 = s0 * (i - j) if abs(h - 1) < 1e-10 else s0 * h * (1 - h ** (i - j)) / (1 - h)
        term2 = s0 * j if abs(b - 1) < 1e-10 else s0 * h ** (i - j) * b * (1 - b**j) / (1 - b)
        a_max = (s0 + term1 + term2) / (i + 1)

        return a_min, a_max

    def _build_average_grids(self) -> list[np.ndarray]:
        """Forward induction: averages grid A[i] of shape (n_averages, i + 1)"""
        p = self.params
        m = self.n_averages
        grids: list[np.ndarray] = [np.zeros((m, i + 1)) for i in range(p.n_steps + 1)]

        for i in range(p.n_steps + 1):
            for j in range(i + 1):
                a_min, a_max = self._average_bounds(i, j)
                if m == 4:
                    a2 = (2 * a_min + a_max) / 3
                    a3 = (a_min + 2 * a_max) / 3
                    grids[i][:, j] = [a_min, a2, a3, a_max]
                else:
                    grids[i][:, j] = np.linspace(a_min, a_max, m)
        return grids

    def price(self, strike: float) -> HullWhiteResult:
        """Price the Asian call for a given strike via backward induction"""
        import time

        p = self.params
        s0, h, b = p.spot, self.up, self.down
        start = time.perf_counter()

        avg = self._build_average_grids()
        value: list[np.ndarray] = [np.zeros((self.n_averages, i + 1)) for i in range(p.n_steps + 1)]

        # terminal payoff on the averages grid
        value[p.n_steps] = np.maximum(avg[p.n_steps] - strike, 0.0)

        for i in range(p.n_steps - 1, -1, -1):
            for j in range(i + 1):
                s_up = s0 * h ** (i - j + 1) * b**j
                s_down = s0 * h ** (i - j) * b ** (j + 1)

                for k in range(self.n_averages):
                    prev_avg = avg[i][k, j]
                    avg_up = ((i + 1) * prev_avg + s_up) / (i + 2)
                    avg_down = ((i + 1) * prev_avg + s_down) / (i + 2)

                    value_up = _interpolate(avg[i + 1][:, j], avg_up, value[i + 1][:, j])
                    value_down = _interpolate(
                        avg[i + 1][:, j + 1], avg_down, value[i + 1][:, j + 1]
                    )

                    value[i][k, j] = (self.q * value_up + (1 - self.q) * value_down) / self.growth

        elapsed = time.perf_counter() - start
        return HullWhiteResult(
            price=float(value[0][0, 0]), n_averages=self.n_averages, elapsed=elapsed
        )
