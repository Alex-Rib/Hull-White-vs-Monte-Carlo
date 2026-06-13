"""End-to-end Asian-option pricing: Monte Carlo vs Hull-White

Loads the AAPL snapshot, prices an arithmetic Asian call at two strikes (ATM
and 5% OTM) with both methods, then studies the convergence of the Hull-White
lattice towards the Monte Carlo benchmark as the number of averages grows
"""

import matplotlib.pyplot as plt
import pandas as pd

from hull_white import HullWhiteAsianPricer
from market_data import load_params
from monte_carlo import MonteCarloAsianPricer

SNAPSHOT_PATH = "data/aapl_snapshot.csv"
N_PATHS = 1_000_000
CONVERGENCE_AVERAGES = [4, 8, 16, 32, 64, 128, 256]


def comparison_table(params, strikes: dict[str, float]) -> pd.DataFrame:
    """Price each strike with both methods and tabulate price and time"""
    rows = []
    mc_pricer = MonteCarloAsianPricer(params, n_paths=N_PATHS)
    for label, strike in strikes.items():
        mc = mc_pricer.price(strike)
        rows.append(
            {"Method": "Monte-Carlo", "Strike": label, "Price": mc.price, "Time": mc.elapsed}
        )
    for label, strike in strikes.items():
        hw = HullWhiteAsianPricer(params, n_averages=4).price(strike)
        rows.append(
            {"Method": "Hull-White", "Strike": label, "Price": hw.price, "Time": hw.elapsed}
        )
    return pd.DataFrame(rows)


def convergence_table(params, strike: float) -> pd.DataFrame:
    """Price with Hull-White for a range of averaging resolutions"""
    rows = []
    for m in CONVERGENCE_AVERAGES:
        hw = HullWhiteAsianPricer(params, n_averages=m).price(strike)
        rows.append({"Averages": m, "Price": hw.price, "Time": hw.elapsed})
    return pd.DataFrame(rows)


def plot_convergence(conv: pd.DataFrame, mc_price: float) -> None:
    """Plot Hull-White price vs number of averages against the MC benchmark"""
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(conv["Averages"], conv["Price"], "o-", color="navy", label="Hull-White")
    ax.axhline(mc_price, color="black", ls="--", lw=0.8, label="Monte Carlo")
    ax.set_xlabel("number of averages")
    ax.set_ylabel("price")
    ax.set_xticks(CONVERGENCE_AVERAGES)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()


def main() -> None:
    """Run the full comparison and convergence study"""
    params = load_params(SNAPSHOT_PATH)
    print("Model parameters")
    print(
        f"  spot={params.spot:.2f} sigma={params.sigma:.4f} "
        f"T={params.maturity:.4f} n_steps={params.n_steps} rate={params.rate:.4f}"
    )

    strikes = {"ATM": params.spot, "ATM + 5%": params.spot * 1.05}

    print("\nMonte Carlo vs Hull-White (n_averages=4)")
    table = comparison_table(params, strikes)
    print(table.to_string(index=False))

    print("\nHull-White convergence to Monte Carlo (ATM)")
    mc_atm = MonteCarloAsianPricer(params, n_paths=N_PATHS).price(params.spot)
    conv = convergence_table(params, params.spot)
    print(conv.to_string(index=False))
    print(f"\nMonte Carlo benchmark (ATM): {mc_atm.price:.4f} +/- {mc_atm.std_error:.4f}")

    plot_convergence(conv, mc_atm.price)
    plt.show()


if __name__ == "__main__":
    main()
