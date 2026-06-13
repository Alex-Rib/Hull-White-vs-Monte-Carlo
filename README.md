# Asian Option Pricing: Monte Carlo vs Hull-White

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Finance](https://img.shields.io/badge/Finance-Path--Dependent%20Options-green)
![Tests](https://img.shields.io/badge/Tests-pytest-purple)
![Lint](https://img.shields.io/badge/Lint-ruff-orange)
![Status](https://img.shields.io/badge/Status-Educational-orange)

## 📊 Description

Two pricing methods for fixed-strike **arithmetic-average Asian calls**, compared on real Apple (AAPL) data:

- **Monte Carlo** simulation of geometric Brownian motion paths
- **Hull-White** recombining binomial lattice with a grid of running averages

The study prices the option at two strikes (ATM and 5% OTM) and analyses how the Hull-White lattice converges to the Monte Carlo benchmark as its averaging grid is refined, illustrating the accuracy-versus-runtime trade-off.

## 🎯 Pipeline

| Stage | Module | Role |
|-------|--------|------|
| Market data | `market_data.py` | Load the AAPL snapshot; derive spot, annualised vol and maturity |
| Monte Carlo | `monte_carlo.py` | GBM path simulation and discounted-payoff pricing |
| Hull-White | `hull_white.py` | Lattice with forward/backward induction on the averages grid |
| Orchestration | `run.py` | Comparison table, convergence study and plot |

## 🔬 Data

The underlying is **AAPL** over its 90 most recent trading days, taken from a committed CSV snapshot for reproducibility. The spot is the latest close, the volatility is annualised from log returns, the maturity is 90/252 years, and the risk-free rate is 4.387%.

## 📐 Methods

### 1. Monte Carlo

Under the risk-neutral measure the spot follows a geometric Brownian motion:

$$dS_t = r\,S_t\,dt + \sigma\,S_t\,dW_t$$

Applying Ito's lemma to the log-price and discretising with an exact log-Euler step:

$$\log S_{t_{k+1}} = \log S_{t_k} + \left(r - \tfrac{1}{2}\sigma^2\right)\Delta t + \sigma\sqrt{\Delta t}\,Z_k, \qquad Z_k \sim \mathcal{N}(0, 1)$$

The Asian call value is the discounted mean payoff on the path average, with averaging over the spot and the simulated steps:

$$C = e^{-rT}\,\mathbb{E}\!\left[\max\!\left(\bar{S} - K,\ 0\right)\right], \qquad \bar{S} = \frac{1}{N+1}\sum_{k=0}^{N} S_{t_k}$$

The estimator is vectorised over paths and reproducible through a seeded RNG; the standard error is reported alongside the price.

### 2. Hull-White lattice

On a recombining binomial tree with up factor $h = e^{\sigma\sqrt{\Delta t}}$, down factor $b = 1/h$ and risk-neutral probability $q = (R - b)/(h - b)$ where $R = e^{r\Delta t}$, each node $(n, j)$ carries a grid of representative running averages

$$A_1(n, j) < A_2(n, j) < \dots < A_M(n, j)$$

bounded by the minimum and maximum attainable average. With $M = 4$ the grid uses the two bounds and two interior thirds; otherwise the averages are spread uniformly.

**Forward induction** propagates the averages. With $j_h$ an up move and $j_b$ a down move:

$$A_i(n+1, j_\bullet) = \frac{1}{n+2}\Bigl((n+1)\,A_i(n, j) + S(n+1, j_\bullet)\Bigr)$$

and for $M = 4$ the interior points are rebuilt as thirds of the refreshed bounds:

$$A_2 = \tfrac{1}{3}\bigl(2 A_1 + A_4\bigr), \qquad A_3 = \tfrac{1}{3}\bigl(A_1 + 2 A_4\bigr)$$

**Terminal payoff** on the averages grid:

$$C_i(N, j) = \max\bigl(A_i(N, j) - K,\ 0\bigr)$$

**Backward induction with linear interpolation.** Rolling back from $(n+1)$ to $(n, j)$, the post-move averages

$$M(\uparrow) = \frac{(n+1)A_i(n,j) + S(n+1, j_h)}{n+2}, \qquad M(\downarrow) = \frac{(n+1)A_i(n,j) + S(n+1, j_b)}{n+2}$$

generally fall between grid points of the child node, so the child option value is linearly interpolated. For the up move, if $M(\uparrow) \in [A_k, A_{k+1}]$ then

$$\lambda_\uparrow = \frac{A_{k+1} - M(\uparrow)}{A_{k+1} - A_k}, \qquad C_i(n,j)^\uparrow = \lambda_\uparrow C_k + (1 - \lambda_\uparrow) C_{k+1}$$

and symmetrically for the down move. The node value discounts the risk-neutral expectation:

$$C_i(n, j) = \frac{1}{R}\Bigl(q\,C_i(n,j)^\uparrow + (1 - q)\,C_i(n,j)^\downarrow\Bigr)$$

## 📈 Results

On the committed AAPL snapshot (spot ≈ 291, vol ≈ 24.4%, T ≈ 0.36y), the Hull-White price converges towards the Monte Carlo benchmark as the number of averages grows:

| Averages | Hull-White price |
|---------:|-----------------:|
| 4 | 25.71 |
| 8 | 18.40 |
| 16 | 14.42 |
| 32 | 12.19 |
| 64 | 11.10 |
| 128 | 10.91 |
| 256 | 10.86 |

Monte Carlo benchmark (ATM, 1M paths): **10.85 ± 0.02**. At M=256 the Hull-White price (10.86) falls inside the Monte Carlo 95% confidence interval: the two independent methods agree, which is the project's main correctness check. The coarse grid (M=4) is far off because the averaging grid is too sparse to track the path-dependent state; refining it closes the gap, at a runtime that roughly doubles each time M doubles. This is the accuracy-versus-cost trade-off the project sets out to show: Monte Carlo is simpler and fast here, while the lattice needs a fine averages grid to match it.

## 🚀 Usage

```bash
python run.py
```

## ✅ Tests & Lint

```bash
pytest -q
ruff check .
```

Tests cover Monte Carlo reproducibility, the Asian-below-vanilla inequality, monotonicity in strike, ordered Hull-White average bounds, and Hull-White convergence to the Monte Carlo price.

## 📚 Reference

Van der Hoek, J. & Elliott, R. J. (2006). *Binomial Models in Finance*. Springer Finance.

## 👨‍💻 Author

Alexandre R. - Université Paris Cité
