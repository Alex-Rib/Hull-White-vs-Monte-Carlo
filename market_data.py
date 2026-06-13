"""Load the underlying price snapshot and derive the pricing parameters

The CSV snapshot (committed for reproducibility) holds the daily close of the
underlying. From it we compute the spot, the annualised volatility from log
returns, and the maturity implied by the number of trading days
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

TRADING_DAYS = 252


@dataclass(frozen=True)
class MarketParams:
    """Model parameters derived from the price snapshot"""

    spot: float  # S0, most recent close
    sigma: float  # annualised volatility from log returns
    maturity: float  # T in years (n_days / 252)
    rate: float  # risk-free rate r
    n_steps: int  # number of time steps (one per trading day)


def load_params(path: str, rate: float = 0.04387, n_days: int = 90) -> MarketParams:
    """Load the snapshot and build the MarketParams

    Keeps the last n_days closes, takes the latest as spot, annualises the
    log-return volatility, and sets the maturity to n_days / 252. Non-numeric
    rows (e.g. a stray ticker header from a malformed export) are dropped
    """
    df = pd.read_csv(path)
    closes = pd.to_numeric(df["Close"], errors="coerce").dropna().to_numpy(dtype=float)
    closes = closes[-n_days:]

    if len(closes) < 2:
        raise ValueError(f"snapshot at {path} has too few valid closes ({len(closes)})")

    spot = float(closes[-1])
    log_returns = np.diff(np.log(closes))
    sigma = float(np.std(log_returns) * np.sqrt(TRADING_DAYS))
    n_steps = len(closes)
    maturity = n_steps / TRADING_DAYS

    return MarketParams(spot=spot, sigma=sigma, maturity=maturity, rate=rate, n_steps=n_steps)
