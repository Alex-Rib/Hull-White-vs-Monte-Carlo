"""Tests for the Monte Carlo and Hull-White Asian-option pricers"""

import numpy as np

from hull_white import HullWhiteAsianPricer
from market_data import MarketParams
from monte_carlo import MonteCarloAsianPricer

PARAMS = MarketParams(spot=100.0, sigma=0.20, maturity=1.0, rate=0.05, n_steps=50)


def test_monte_carlo_is_reproducible():
    """Same seed gives the same price"""
    a = MonteCarloAsianPricer(PARAMS, n_paths=50_000, seed=7).price(100.0)
    b = MonteCarloAsianPricer(PARAMS, n_paths=50_000, seed=7).price(100.0)
    assert a.price == b.price


def test_asian_cheaper_than_vanilla():
    """An arithmetic Asian call is cheaper than the vanilla call

    Averaging lowers the variance of the terminal quantity, so the Asian call
    has less optionality and a lower price than the European call on the spot
    """
    pricer = MonteCarloAsianPricer(PARAMS, n_paths=200_000, seed=1)
    paths = pricer.simulate_paths()
    discount = np.exp(-PARAMS.rate * PARAMS.maturity)

    asian = discount * np.maximum(paths.mean(axis=1) - 100.0, 0.0).mean()
    vanilla = discount * np.maximum(paths[:, -1] - 100.0, 0.0).mean()
    assert asian < vanilla


def test_payoff_is_non_negative():
    """A call price is non-negative for any strike"""
    pricer = MonteCarloAsianPricer(PARAMS, n_paths=20_000, seed=3)
    for strike in (80.0, 100.0, 130.0):
        assert pricer.price(strike).price >= 0.0


def test_price_decreases_with_strike():
    """Call price is non-increasing in strike"""
    pricer = MonteCarloAsianPricer(PARAMS, n_paths=100_000, seed=2)
    p_low = pricer.price(90.0).price
    p_high = pricer.price(110.0).price
    assert p_low > p_high


def test_hull_white_average_bounds_ordered():
    """A_min does not exceed A_max at every node"""
    pricer = HullWhiteAsianPricer(PARAMS, n_averages=4)
    for i in range(PARAMS.n_steps + 1):
        for j in range(i + 1):
            a_min, a_max = pricer._average_bounds(i, j)
            assert a_min <= a_max + 1e-9


def test_hull_white_converges_to_monte_carlo():
    """Hull-White approaches the MC price as the averaging grid refines

    The coarse grid (M=4) is far off; refining to M=48 must move the price
    much closer to the Monte Carlo benchmark
    """
    mc = MonteCarloAsianPricer(PARAMS, n_paths=300_000, seed=5).price(100.0).price
    coarse = HullWhiteAsianPricer(PARAMS, n_averages=4).price(100.0).price
    fine = HullWhiteAsianPricer(PARAMS, n_averages=48).price(100.0).price
    assert abs(fine - mc) < abs(coarse - mc)
