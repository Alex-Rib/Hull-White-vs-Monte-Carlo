

# ============================================================================
# 1 IMPORTS
# ============================================================================
import pandas as pd
import yfinance as yf
import numpy as np
from datetime import datetime, timedelta
import time
import matplotlib.pyplot as plt


# ============================================================================
# 4 FONCTIONS SIMULATION ET PRICING
# ============================================================================

# Simulation de chemins de prix via la méthode de Monte Carlo
def simulate_monte_carlo_paths(S0, r, sigma, T, N, M):
    dt = T / N  # Taille du pas de temps
    Z = np.random.randn(M, N)  # M chemins, N pas de temps
    log_prices_acrroissement = (r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * Z  # Incréments log-normaux
    log_paths = np.cumsum(log_prices_acrroissement, axis=1)  # Chemins log-normaux
    Snew = np.hstack([np.zeros((M, 1)), log_paths]) # Ajouter S0 au début (log(1)=0)
    prices = S0 * np.exp(Snew)  # Convertir en prix
    return prices


# MONTE CARLO - PRICING D'UNE OPTION ASIATIQUE CALL
# avec M simulations
def price_asian_call_monte_carlo_all_in_one(S0, K, r, sigma, T, N, M):
    t0 = time.time()
    paths_matrix = simulate_monte_carlo_paths(S0, r, sigma, T, N, M)
    S_avg = np.mean(paths_matrix, axis=1) 
    payoffs = np.maximum(S_avg - K, 0)
    price = np.exp(-r * T) * np.mean(payoffs)
    t1 = time.time()
    return {
        'price': price,
        'time': t1 - t0
    }


# HULL & WHITE - PRICING D'UNE OPTION ASIATIQUE CALL
# avec M moyennes par nœud
def price_asian_call_HW_all_in_one(S0, K, r, sigma, T, N, M):
    def interpolate_safe(avg_grid, target_avg, child_vals):
        # On enlève les NA
        valid = ~np.isnan(avg_grid) & ~np.isnan(child_vals) # indices valides
        grid_vals = avg_grid[valid]  # valeurs A_k de la grille qui sont ok
        grid_children = child_vals[valid]  # valeurs C_k de la grille qui sont ok
        
        # On trie par ordre croissant
        order_idx = np.argsort(grid_vals) # indices de tri
        sorted_grid = grid_vals[order_idx]  # A_k triés
        sorted_children = grid_children[order_idx]  # C_k triés
        
        # Tolérance basée sur la précision machine
        tolerance = np.sqrt(np.finfo(float).eps)
        
        # Position de target_avg dans l'intervalle [sorted_grid[idx], sorted_grid[idx+1]]
        grid_index = np.searchsorted(sorted_grid, target_avg, side='right') - 1
        # On borne pour ne pas sortir de l'intervalle : [0, length-2]
        grid_index = max(0, min(grid_index, len(sorted_grid) - 2))
        
        # Distance entre deux points (longueur de l'intervalle [A_k, A_k+1])
        delta_A = sorted_grid[grid_index + 1] - sorted_grid[grid_index]
        
        if abs(delta_A) < tolerance:
            # Si les deux points sont presque identiques on retourne la valeur
            return sorted_children[grid_index]
        else:
            # Coef pour mesurer la distance entre la valeur et l'un des points
            lambda_coef = (sorted_grid[grid_index + 1] - target_avg) / delta_A # dans [0,1]
            return (lambda_coef * sorted_children[grid_index] + 
                    (1 - lambda_coef) * sorted_children[grid_index + 1]) 
    
    # Paramètres du modèle
    delta_t = T / N
    h = np.exp(sigma * np.sqrt(delta_t)) # facteur de hausse (u)
    b = 1 / h # facteur de baisse (d)
    R = np.exp(r * delta_t)
    Q = (R - b) / (h - b)
    
    t0 = time.time()  # Début du chrono
    
    # Initialisation des grilles de A_k et de C_k
    Agrid = [None] * (N + 1)  # Stock les moyennes
    Cgrid = [None] * (N + 1)  # Stock les valeurs de l'option
    
    # Forward induction pour les bornes Amin et Amax
    for i in range(N + 1):
        Agrid[i] = np.zeros((M, i + 1))  # Agrid[i][k, j] stocke le kème A_k(i,j)
        Cgrid[i] = np.zeros((M, i + 1))  # Cgrid[i][k, j] stocke C_k(i,j)
        
        for j in range(i + 1):
            # Calcul de la plus petite moyenne Amin
            if abs(b - 1) < 1e-10:
                term1 = S0 * j
            else:
                term1 = S0 * b * (1 - b**j) / (1 - b)
            
            if abs(h - 1) < 1e-10:
                term2 = S0 * (i - j)
            else:
                term2 = S0 * b**j * h * (1 - h**(i - j)) / (1 - h)
            
            A_min = (S0 + term1 + term2) / (i + 1)
            
            # Calcul de la plus grande moyenne Amax
            if abs(h - 1) < 1e-10:
                term1 = S0 * (i - j)
            else:
                term1 = S0 * h * (1 - h**(i - j)) / (1 - h)
            
            if abs(b - 1) < 1e-10:
                term2 = S0 * j
            else:
                term2 = S0 * h**(i - j) * b * (1 - b**j) / (1 - b)
            
            A_max = (S0 + term1 + term2) / (i + 1)
            
            if M == 4:
                # Points définis comme dans le livre
                A1 = A_min
                A4 = A_max
                A2 = (2 * A1 + A4) / 3
                A3 = (A1 + 2 * A4) / 3
                Agrid[i][:, j] = [A1, A2, A3, A4]
            else:
                # Répartition entre A_min et A_max (équidistance)
                Agrid[i][:, j] = np.linspace(A_min, A_max, M)
    
    # Calcul du payoff pour chacun des A_k(N,j)
    Cgrid[N] = np.maximum(Agrid[N] - K, 0)
    
    # Backward induction
    for i in range(N - 1, -1, -1):
        Cgrid[i] = np.zeros((M, i + 1))
        
        for j in range(i + 1):
            # Valeurs de S après une hausse ou une baisse
            S_h = S0 * h**(i - j + 1) * b**j  # S(n+1,j_h)
            S_b = S0 * h**(i - j) * b**(j + 1)  # S(n+1,j_b)
            
            for k in range(M):
                prev_avg = Agrid[i][k, j]  # Moyenne partielle au nœud (i,j)
                
                # Nouvelles moyennes après une hausse ou une baisse
                M_h = ((i + 1) * prev_avg + S_h) / (i + 2)
                M_b = ((i + 1) * prev_avg + S_b) / (i + 2)
                
                # Interpolation dans la colonne j (hausse) au temps i+1
                # value_h = C_i(n,j)^↑
                value_h = interpolate_safe(
                    avg_grid=Agrid[i + 1][:, j],
                    target_avg=M_h,
                    child_vals=Cgrid[i + 1][:, j]
                )
                
                # Interpolation dans la colonne j+1 (baisse) au temps i+1
                # value_b = C_i(n,j)^↓
                value_b = interpolate_safe(
                    avg_grid=Agrid[i + 1][:, j + 1],
                    target_avg=M_b,
                    child_vals=Cgrid[i + 1][:, j + 1]
                )
                
                # Valeur actualisée via l'espérance sous Q
                Cgrid[i][k, j] = (Q * value_h + (1 - Q) * value_b) / R
    
    t1 = time.time()  # On stop le chrono
    
    return {
        'price': Cgrid[0][0, 0],  # Valeur de l'option au temps 0
        'time': t1 - t0,
        'M': M
    }


def main():
    # ============================================================================
    # 2 IMPORTATION ET PRÉPARATION DES DONNÉES AAPL
    # ============================================================================

    # Import des données AAPL des 140 derniers jours depuis Yahoo Finance
    end_date = datetime.now()  # Date actuelle
    start_date = end_date - timedelta(days=140)  
    data = yf.download('AAPL', start=start_date, end=end_date)  
    data = data.tail(90).copy()  # On garde les 90 derniers jours uniquement
    data.reset_index(inplace=True)  # Réinitialise l'index pour avoir la date comme colonne
    data['Date'] = pd.to_datetime(data['Date'])  # S'assure que la colonne Date est au format datetime
    data['Day'] = range(1, len(data) + 1)  # Ajoute une colonne 'Day' pour le numéro du jour


    # ============================================================================
    # 3 CALCUL DES PARAMÈTRES DU MODÈLE
    # ============================================================================


    close_prices = data['Close'].values.flatten()  # Convertion en array 1D
    S0 = close_prices[-1]  # Dernier prix de clôture
    K_ATM = S0  # Strike ATM
    K_5 = S0 * 1.05  # Strike OTM en espérant une hausse
    N = len(data)  # Nombre de pas
    T = N / 252  # Maturité (une année boursière = 252 jours [marché US])

    # Volatilité annualisée
    log_returns = np.diff(np.log(close_prices))
    sigma = np.std(log_returns) * np.sqrt(252)

    # Taux d'intérêt sans risque
    r = 0.04387


    # ============================================================================
    # 5 PRICING DES OPTIONS ET DATAFRAMES 
    # ============================================================================

    # Prix via Monte Carlo - Strike ATM (1000000 simulations)
    resultMC_ATM = price_asian_call_monte_carlo_all_in_one(S0, K_ATM, r, sigma, T, N, M=1000000)
    # Création du DataFrame
    MC_ATM_Data_frame = pd.DataFrame([{
        'Type': 'Call Asiatique',
        'Price': resultMC_ATM['price'],
        'Time': resultMC_ATM['time']
    }])

    # Prix via Monte Carlo strike otm +5% (1000000 simulations)
    resultMC_5 = price_asian_call_monte_carlo_all_in_one(S0, K_5, r, sigma, T, N, M=1000000)
    # Création du DataFrame
    MC_5_Data_frame = pd.DataFrame([{
        'Type': 'Call Asiatique',
        'Price': resultMC_5['price'],
        'Time': resultMC_5['time']
    }])

    # Prix via Hull & White strike atm (M=4)
    resultHW_ATM = price_asian_call_HW_all_in_one(S0, K_ATM, r, sigma, T, N, M=4)
    # Création du DataFrame
    HW_ATM_Data_frame = pd.DataFrame([{
        'Type': 'Call Asiatique',
        'Price': resultHW_ATM['price'],
        'Time': resultHW_ATM['time']
    }])
    # Prix via Hull & White strike otm +5% (M=4)
    resultHW_5 = price_asian_call_HW_all_in_one(S0, K_5, r, sigma, T, N, M=4)
    # Création du DataFrame
    HW_5_Data_frame = pd.DataFrame([{
        'Type': 'Call Asiatique',
        'Price': resultHW_5['price'],
        'Time': resultHW_5['time']
    }]) 

    # Tableau comparatif des résultats
    MC_HW_combined_data_frame = pd.DataFrame({
        'Type': ['Call Asiatique'] * 4,
        'Method': ['Monte-Carlo', 'Monte-Carlo', 'Hull & White', 'Hull & White'],
        'Strike': ['ATM', 'ATM + 5%', 'ATM', 'ATM + 5%'],
        'Price': [
            resultMC_ATM['price'],
            resultMC_5['price'],
            resultHW_ATM['price'],
            resultHW_5['price']
        ],
        'Time': [
            resultMC_ATM['time'],
            resultMC_5['time'],
            resultHW_ATM['time'],
            resultHW_5['time']
        ]
    })

    # ============================================================================
    # 6 CONVERGENCE HULL & WHITE VERS MONTE CARLO
    # ============================================================================

    M_conv = [4, 8, 16, 32, 64, 128]  # les différentes moyennes à tester

    # DataFrame pour stocker les résultats en fonction de M (l'éxécution est lente)
    Convergence_HW = pd.DataFrame([
        {
            'Type': 'Call asiatique',
            'Method': 'Hull & White',
            'Strike': 'ATM',
            'Price': (result := price_asian_call_HW_all_in_one(S0, K_ATM, r, sigma, T, N, M=Ma))['price'],
            'Time': result['time'],
            'Average Number': Ma
        }
        for Ma in M_conv
    ])



    # ============================================================================
    # 7 PLOT et  AFFICHAGE DES RESULTATS
    # ============================================================================


    print("\n---- Paramètres du modèle : ----")
    print(f"S0: {S0:.2f}")
    print(f"Strike ATM: {K_ATM:.2f}")
    print(f"Strike OTM: {K_5:.2f}")
    print(f"Nombre de pas: {N}")
    print(f"Maturité en année: {T:.2f}")
    print(f"Volatilité annualisée: {sigma:.4f}")
    print(f"Taux d'intérêt sans risque: {r:.4f}")

    print("\n---- Aperçu des données ----")
    print(data.head(10))
    print("\n---- Évolution des prix de clôture d'AAPL sur les 90 derniers jours ----")
    plt.figure(figsize=(10, 6))
    plt.plot(data['Day'], data['Close'], linewidth=1.3, color='purple')
    plt.xlabel('Day')
    plt.ylabel('Price')
    plt.xlim(1, 90)
    plt.xticks(range(1, 91, 5))
    plt.grid(True, alpha=0.3)
    plt.show()


    print("\n---- Monte Carlo paths ----")
    M = 100 # nombre de chemins
    paths = simulate_monte_carlo_paths(S0, r, sigma, T, N, M)
    plt.figure(figsize=(12, 6))
    # Trace chaque chemin
    for i in range(M):
        plt.plot(range(N+1), paths[i, :], alpha=0.9, linewidth=0.5)
    # Tracer le prix initial
    plt.axhline(y=S0, color='black', linestyle='--', linewidth=1, label=f'S0 = {S0:.2f}')
    plt.title(f'{M} trajectoires')
    plt.xlabel('time steps')
    plt.ylabel('Prices')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    print("\n---- Tableau comparatif entre Monte Carlo (1000000 simulations) et Hull & White (4 moyennes) ----")
    print(MC_HW_combined_data_frame)


    print("\n--- Convergence d'Hull & White vers Monte Carlo ---")
    plt.figure(figsize=(10, 6))
    plt.plot(Convergence_HW['Average Number'], Convergence_HW['Price'], linewidth=1.1, color='blue', marker='o')
    plt.axhline(y=resultMC_ATM['price'], color='black', linewidth=0.7, linestyle='--', label='Monte Carlo price')
    plt.xlabel('Average Number', fontsize=12)
    plt.ylabel('Price', fontsize=12)
    plt.xlim(4, 128)
    plt.xticks([4, 8, 16, 32, 64, 128])
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


    print("\n---- Temps d'exécution d'Hull & White en fonction du nombre de moyennes ----")
    plt.figure(figsize=(10, 6))
    plt.plot(Convergence_HW['Time'], Convergence_HW['Average Number'], linewidth=1.1, color='orange', marker='o')
    plt.xlabel('Time (seconds)', fontsize=12)
    plt.ylabel('Average Number', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    print("\n---- Convergence d'Hull & White vers Monte Carlo en détail ----")
    print(Convergence_HW)


if __name__ == "__main__":
    main()
