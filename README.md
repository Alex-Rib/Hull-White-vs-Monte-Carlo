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

$$dS_t = r\,S_t\,dt + \sigma\,S_t\,dW_t$$

- **Lemme d'Itô** appliqué au log-prix : 

$$d\!\bigl(\log S_t\bigr) = \Bigl(r - \tfrac{1}{2}\sigma^{2}\Bigr)\,dt + \sigma\,dW_t$$

- **Méthode d'Euler** pour la simulation : 

$$\log S_{t_{k+1}} = \log S_{t_k} + \Bigl(r - \tfrac{1}{2}\sigma^{2}\Bigr)\,\Delta t + \sigma\sqrt{\Delta t}\,Z_k, \qquad Z_k \sim \mathcal{N}(0,1)$$

#### 2. Hull & White

- **Arbre binomial recombinant** avec grille de moyennes :

$$A_1(n,j) < A_2(n,j) < A_3(n,j) < A_4(n,j)$$

- **Calcul des $A_i$ par forward induction** : 

$$\begin{aligned}
A_i(0,0) &= S(0,0) \quad \forall i \in \{1,2,3,4\}\\[1em]
A_1(n\!+\!1,j_h) &= \tfrac{1}{n+2}\bigl((n+1)A_1(n,j) + S(n+1,j_h)\bigr)\\[1em] 
A_1(n\!+\!1,j_b) &= \tfrac{1}{n+2}\bigl((n+1)A_1(n,j) + S(n+1,j_b)\bigr)\\[1em]
A_4(n\!+\!1,j_h) &= \tfrac{1}{n+2}\bigl((n+1)A_4(n,j) + S(n+1,j_h)\bigr)\\[1em]
A_4(n\!+\!1,j_b) &= \tfrac{1}{n+2}\bigl((n+1)A_4(n,j) + S(n+1,j_b)\bigr)\\[1em]
A_2(n\!+\!1,j) &= \tfrac{1}{3}\bigl(2\,A_1(n\!+\!1,j)+A_4(n\!+\!1,j)\bigr)\\[1em]
A_3(n\!+\!1,j) &= \tfrac{1}{3}\bigl(A_1(n\!+\!1,j)+2\,A_4(n\!+\!1,j)\bigr)
\end{aligned}$$

avec $j_h$ = hausse et $j_b$ = baisse.

- **Calcul des payoffs** : 

$$C_i(N,j) = \max\bigl(A_i(N,j)-K,\ 0\bigr), \qquad i=1,\ldots,4$$

- **Backward induction et interpolation linéaire** :

Calcul des moyennes conditionnelles :

$$\begin{aligned}
M(\uparrow) &= \tfrac{1}{n+2}\bigl((n+1)A_i(n,j) + S(n+1,j_h)\bigr)\\[1em]
M(\downarrow) &= \tfrac{1}{n+2}\bigl((n+1)A_i(n,j) + S(n+1,j_b)\bigr)
\end{aligned}$$

avec $M(\uparrow)$ = moyenne si hausse, $M(\downarrow)$ = moyenne si baisse.

- **En cas de hausse** :

S'il existe $k$ tel que $M(\uparrow) \in [A_k(n+1,j_h), A_{k+1}(n+1,j_h)]$, alors :

$$\lambda_{\uparrow} = \frac{A_{k+1}(n+1,j_h) - M(\uparrow)}{A_{k+1}(n+1,j_h) - A_k(n+1,j_h)}$$

- **En cas de baisse** :

S'il existe $t$ tel que $M(\downarrow) \in [A_t(n+1,j_b), A_{t+1}(n+1,j_b)]$, alors :

$$\lambda_{\downarrow} = \frac{A_{t+1}(n+1,j_b) - M(\downarrow)}{A_{t+1}(n+1,j_b) - A_t(n+1,j_b)}$$

- **Valeurs interpolées** :

$$\begin{aligned}
C_i(n,j)^{\uparrow} &= \lambda_{\uparrow}C_k(n+1,j_h) + (1-\lambda_{\uparrow})C_{k+1}(n+1,j_h)\\[1em]
C_i(n,j)^{\downarrow} &= \lambda_{\downarrow}C_t(n+1,j_b) + (1-\lambda_{\downarrow})C_{t+1}(n+1,j_b)
\end{aligned}$$

- **Valeur actualisée** : 

$$C_i(n,j) = \tfrac{1}{R(n,j)}\bigl(\pi(n,j)C_i(n,j)^{\uparrow} + (1-\pi(n,j))C_i(n,j)^{\downarrow}\bigr)$$

avec $R(n,j)$ le facteur d'actualisation et $\pi(n,j) = \dfrac{R(n,j) - b(n,j)}{h(n,j) - b(n,j)}$ la probabilité risque-neutre.

## 📈 Résultats

Le programme génère :

1. **Tableau comparatif** : Prix et temps d'exécution pour les deux méthodes
2. **Graphiques de convergence** : Évolution du prix Hull & White selon le nombre de moyennes (M)
3. **Visualisations** :
   - Évolution du prix AAPL sur 90 jours
   - 100 trajectoires Monte Carlo
   - Convergence vers le prix Monte Carlo
   - Analyse du temps d'exécution

## 🚀 Installation

### Prérequis
- Python 3.8+
- pip

### Dépendances

```bash
pip install -r requirements.txt
```

Les packages nécessaires :
- `numpy` : calculs numériques
- `pandas` : manipulation de données
- `matplotlib` : visualisations
- `yfinance` : récupération des données financières

## 💻 Utilisation

```bash
python asian_options_pricing.py
```

Le script va :
1. Télécharger les données AAPL
2. Calculer les paramètres du modèle
3. Exécuter les simulations Monte Carlo
4. Calculer les prix via Hull & White
5. Afficher les résultats et graphiques

## 📚 Références

- Hull, J. C., & White, A. (1993). *Efficient Procedures for Valuing European and American Path-Dependent Options*. Journal of Derivatives, 1(1), 21-31.

- Hoek, J. van der, & Elliott, R. J. (2006). *Binomial Models in Finance*. Springer Finance. New York: Springer-Verlag. ISBN: 978-0-387-25898-0.

## 🔍 Analyse de convergence

Le programme teste la convergence de Hull & White avec M = [4, 8, 16, 32, 64, 128] moyennes par nœud.

**Observations** :
- Plus M augmente, plus le prix H&W converge vers Monte Carlo
- Le temps de calcul augmente de manière exponentielle avec M
- Un bon compromis est généralement M = 16 ou 32
## 👨‍💻 Auteur

Alexandre R. - Master ISIFAR, Université Paris Cité

## 📄 Licence

Ce projet est à usage académique et pédagogique.

## ⚠️ Avertissement

Ce code est destiné à des fins éducatives uniquement. Ne pas utiliser pour des décisions d'investissement réelles sans validation appropriée.

