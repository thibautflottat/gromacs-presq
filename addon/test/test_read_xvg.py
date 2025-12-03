#!/usr/bin/env python3
"""Test reading complex_obs.xvg and compute F(q,t)"""

import numpy as np
import matplotlib.pyplot as plt

# Lire le fichier XVG (ignorer les lignes de commentaire)
data = np.loadtxt('complex_obs.xvg', comments=['#', '@'])

print(f"Shape des données : {data.shape}")
print(f"Nombre de frames : {data.shape[0]}")
print(f"Nombre de colonnes : {data.shape[1]}")

# Extraire le temps
time = data[:, 0]
print(f"Temps : {time[0]:.3f} à {time[-1]:.3f} ps ({len(time)} frames)")

# Extraire les parties réelles et imaginaires (50 q-vecteurs)
nq = 50
real_parts = data[:, 1::2]  # Colonnes impaires : Re(q0), Re(q1), ...
imag_parts = data[:, 2::2]  # Colonnes paires : Im(q0), Im(q1), ...

print(f"\nParties réelles shape : {real_parts.shape}")
print(f"Parties imaginaires shape : {imag_parts.shape}")

# Calculer n(q,t) complexe
n_qt = real_parts + 1j * imag_parts

print(f"\nn(q,t) shape : {n_qt.shape}")
print(f"Premier n(q0,t=0) = {n_qt[0, 0]:.2f}")
print(f"Norme |n(q0,t=0)| = {np.abs(n_qt[0, 0]):.2f}")

# Calculer F(q,t) = <n*(q,0) n(q,t)> pour un q-vecteur
# Ici on fait juste n*(q,0) * n(q,t) sans moyenne temporelle
q_idx = 0  # Premier q-vecteur
n_q0 = n_qt[0, q_idx]  # n(q, t=0)
F_qt = np.conj(n_q0) * n_qt[:, q_idx]  # F(q,t) pour ce q

print(f"\nF(q={q_idx}, t) calculé")
print(f"F(q,0) = {F_qt[0]:.2f}")
print(f"|F(q,0)| = {np.abs(F_qt[0]):.2f}")

# Moyenner sur les 50 q-vecteurs
F_avg = np.zeros(len(time), dtype=complex)
for i in range(nq):
    n_q0_i = n_qt[0, i]
    F_avg += np.conj(n_q0_i) * n_qt[:, i]
F_avg /= nq

print(f"\nF(q,t) moyenné sur {nq} q-vecteurs")
print(f"F_avg(0) = {F_avg[0]:.2f}")
print(f"|F_avg(0)| = {np.abs(F_avg[0]):.2f}")

# Tracer |F(q,t)| pour vérifier la décroissance
plt.figure(figsize=(10, 6))
plt.subplot(2, 1, 1)
plt.plot(time, np.abs(F_avg), 'o-')
plt.xlabel('Time (ps)')
plt.ylabel('|F(q,t)|')
plt.title('Intermediate Scattering Function (averaged over 50 q-vectors)')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(time, F_avg.real, 'o-', label='Re[F(q,t)]')
plt.plot(time, F_avg.imag, 's-', label='Im[F(q,t)]')
plt.xlabel('Time (ps)')
plt.ylabel('F(q,t)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('F_qt.png', dpi=150)
print(f"\nGraphique sauvegardé : F_qt.png")
