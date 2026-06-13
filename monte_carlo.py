"""Monte Carlo pricer for arithmetic-average Asian call options

Simulates geometric Brownian motion paths under the risk-neutral measure and
prices a fixed-strike arithmetic Asian call as the discounted expected payoff
on the path average. Vectorised over paths; reproducible via a seeded RNG
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class MonteCarloResult:
    """Outcome of a Monte Carlo pricing run"""

    price: float
    std_error: float
    elapsed: float


class MonteCarloAsianPricer:
    """Monte Carlo pricer for an arithmetic-average Asian call

    The path average includes the spot and the n simulated steps, matching the
    discrete averaging convention used by the Hull-White lattice
    """

    def __init__(self, params, *, n_paths: int = 1_000_000, seed: int = 42) -> None:
        """Store the market parameters and the simulation settings

        params : MarketParams
            Spot, vol, maturity, rate and number of steps
        n_paths : int
            Number of simulated paths
        seed : int
            RNG seed for reproducibility
        """
        self.params = params
        self.n_paths = n_paths
        self._rng = np.random.default_rng(seed)

    def simulate_paths(self) -> np.ndarray:
        """Simulate price paths under the risk-neutral GBM

        Returns an array of shape (n_paths, n_steps + 1) including the spot at
        column 0, using the exact log-Euler step for GBM
        """
        p = self.params
        dt = p.maturity / p.n_steps
        drift = (p.rate - 0.5 * p.sigma**2) * dt
        diffusion = p.sigma * np.sqrt(dt)

        z = self._rng.standard_normal((self.n_paths, p.n_steps))
        log_increments = drift + diffusion * z
        log_paths = np.cumsum(log_increments, axis=1)
        log_paths = np.hstack([np.zeros((self.n_paths, 1)), log_paths])
        return p.spot * np.exp(log_paths)

    def price(self, strike: float) -> MonteCarloResult:
        """Price the Asian call for a given strike

        Returns the discounted mean payoff, the standard error of the estimate,
        and the elapsed time
        """
        import time

        p = self.params
        start = time.perf_counter()

        paths = self.simulate_paths()
        average = paths.mean(axis=1)
        payoffs = np.maximum(average - strike, 0.0)
        discount = np.exp(-p.rate * p.maturity)

        price = discount * payoffs.mean()
        std_error = discount * payoffs.std(ddof=1) / np.sqrt(self.n_paths)
        elapsed = time.perf_counter() - start

        return MonteCarloResult(price=float(price), std_error=float(std_error), elapsed=elapsed)
