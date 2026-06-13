"""Fetch the AAPL daily close from yfinance and save a snapshot to CSV

Run once to refresh data/aapl_snapshot.csv. The pipeline reads the snapshot,
not the live feed, so results stay reproducible
"""

from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import yfinance as yf


def main(ticker: str = "AAPL", n_days: int = 90) -> None:
    """Download recent closes and write a clean two-column snapshot

    yfinance returns a MultiIndex column frame for a single ticker, so the
    close series is extracted explicitly and written as a flat Date,Close CSV
    """
    end = datetime.now()
    start = end - timedelta(days=n_days * 2)  # buffer for non-trading days
    raw = yf.download(ticker, start=start, end=end, auto_adjust=True, progress=False)

    if raw.empty:
        raise RuntimeError(f"no data returned for {ticker} (check ticker or connection)")

    # handle both the MultiIndex (single-ticker) and flat layouts
    close = raw["Close"]
    if isinstance(close, pd.DataFrame):
        close = close.iloc[:, 0]

    snapshot = close.tail(n_days).rename("Close").reset_index().rename(columns={"index": "Date"})
    snapshot.columns = ["Date", "Close"]

    Path("data").mkdir(exist_ok=True)
    out = Path("data") / f"{ticker.lower()}_snapshot.csv"
    snapshot.to_csv(out, index=False)
    print(f"saved {len(snapshot)} rows -> {out}")


if __name__ == "__main__":
    main()
