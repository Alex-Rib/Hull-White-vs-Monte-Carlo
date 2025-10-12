# Pricing d'Options Asiatiques : Monte Carlo vs Hull & White

## 📊 Description

Ce projet compare deux méthodes de pricing pour les calls asiatiques :
- **Monte Carlo** 
- **Hull & White** 

L'étude utilise les données réelles d'**Apple (AAPL)** pour calibrer les paramètres du modèle.

## 🎯 Objectifs

- Implémenter et comparer deux approches de pricing d'options asiatiques
- Analyser la convergence d'Hull & White vers Monte Carlo
- Évaluer le trade-off entre précision et temps de calcul


## 🔬 Méthodologie

### Données
- **Actif sous-jacent** : Apple (AAPL)
- **Période** : 90 derniers jours de bourse
- **Source** : Yahoo Finance (via `yfinance`)

### Paramètres du modèle
- **S₀** : Prix de clôture le plus récent
- **K** : Strike ATM et OTM (+5%)
- **T** : Maturité (90 jours / 252 jours ouvrés)
- **σ** : Volatilité annualisée (calculée sur les rendements logarithmiques)
- **r** : Taux sans risque (4.387%)

### Méthodes implémentées

#### 1. Monte Carlo

- **Hypothèse** : mouvement brownien géométrique (Black-Scholes) : 
```math
$$dS_t = r\,S_t\,dt + \sigma\,S_t\,dW_t$$
```
- **Lemme d'Itô** appliqué au log-prix : 
```math
$$d\!\bigl(\log S_t\bigr) = \Bigl(r - \tfrac{1}{2}\sigma^{2}\Bigr)\,dt + \sigma\,dW_t$$
```
- **Méthode d'Euler** pour la simulation : 
```math
$$\log S_{t_{k+1}} = \log S_{t_k} + \Bigl(r - \tfrac{1}{2}\sigma^{2}\Bigr)\,\Delta t + \sigma\sqrt{\Delta t}\,Z_k, \qquad Z_k \sim \mathcal{N}(0,1)$$
```
#### 2. Hull & White

- **Arbre binomial recombinant** avec grille de moyennes :

$$A_1(n,j) < A_2(n,j) < A_3(n,j) < A_4(n,j)$$

- **Calcul des $A_i$ par forward induction** : 
```math
$$\begin{aligned}
A_i(0,0) &= S(0,0) \quad \forall i \in \{1,2,3,4\}\\[1em]
A_1(n\!+\!1,j_h) &= \tfrac{1}{n+2}\bigl((n+1)A_1(n,j) + S(n+1,j_h)\bigr)\\[1em] 
A_1(n\!+\!1,j_b) &= \tfrac{1}{n+2}\bigl((n+1)A_1(n,j) + S(n+1,j_b)\bigr)\\[1em]
A_4(n\!+\!1,j_h) &= \tfrac{1}{n+2}\bigl((n+1)A_4(n,j) + S(n+1,j_h)\bigr)\\[1em]
A_4(n\!+\!1,j_b) &= \tfrac{1}{n+2}\bigl((n+1)A_4(n,j) + S(n+1,j_b)\bigr)\\[1em]
A_2(n\!+\!1,j) &= \tfrac{1}{3}\bigl(2\,A_1(n\!+\!1,j)+A_4(n\!+\!1,j)\bigr)\\[1em]
A_3(n\!+\!1,j) &= \tfrac{1}{3}\bigl(A_1(n\!+\!1,j)+2\,A_4(n\!+\!1,j)\bigr)
\end{aligned}$$
```
avec $j_h$ = hausse et $j_b$ = baisse.

- **Calcul des payoffs** : 
```math
$$C_i(N,j) = \max\bigl(A_i(N,j)-K,\ 0\bigr), \qquad i=1,\ldots,4$$
```
- **Backward induction et interpolation linéaire** :

Calcul des moyennes conditionnelles :
```math
$$\begin{aligned}
M(\uparrow) &= \tfrac{1}{n+2}\bigl((n+1)A_i(n,j) + S(n+1,j_h)\bigr)\\[1em]
M(\downarrow) &= \tfrac{1}{n+2}\bigl((n+1)A_i(n,j) + S(n+1,j_b)\bigr)
\end{aligned}$$
```
avec $M(\uparrow)$ = moyenne si hausse, $M(\downarrow)$ = moyenne si baisse.

- **En cas de hausse** :

S'il existe $k$ tel que $M(\uparrow) \in [A_k(n+1,j_h), A_{k+1}(n+1,j_h)]$, alors :

$$\lambda_{\uparrow} = \frac{A_{k+1}(n+1,j_h) - M(\uparrow)}{A_{k+1}(n+1,j_h) - A_k(n+1,j_h)}$$

- **En cas de baisse** :

S'il existe $t$ tel que $M(\downarrow) \in [A_t(n+1,j_b), A_{t+1}(n+1,j_b)]$, alors :

$$\lambda_{\downarrow} = \frac{A_{t+1}(n+1,j_b) - M(\downarrow)}{A_{t+1}(n+1,j_b) - A_t(n+1,j_b)}$$

- **Valeurs interpolées** :
```math
$$\begin{aligned}
C_i(n,j)^{\uparrow} &= \lambda_{\uparrow}C_k(n+1,j_h) + (1-\lambda_{\uparrow})C_{k+1}(n+1,j_h)\\[1em]
C_i(n,j)^{\downarrow} &= \lambda_{\downarrow}C_t(n+1,j_b) + (1-\lambda_{\downarrow})C_{t+1}(n+1,j_b)
\end{aligned}$$
```
- **Valeur actualisée** : 

$$C_i(n,j) = \tfrac{1}{R(n,j)}\bigl(\pi(n,j)C_i(n,j)^{\uparrow} + (1-\pi(n,j))C_i(n,j)^{\downarrow}\bigr)$$

avec $R(n,j)$ le facteur d'actualisation et $\pi(n,j) = \dfrac{R(n,j) - b(n,j)}{h(n,j) - b(n,j)}$ la probabilité risque-neutre.



## 📚 Références

- Hoek, J. van der, & Elliott, R. J. (2006). *Binomial Models in Finance*. Springer Finance


## 👨‍💻 Auteur

Alexandre R. - Master ISIFAR, Université Paris Cité


