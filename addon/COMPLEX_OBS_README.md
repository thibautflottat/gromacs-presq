# Complex Observable Output - Fichier XVG

## Description

GROMACS génère automatiquement un fichier `complex_obs.xvg` contenant les valeurs de **n(q,t)** pour 50 q-vecteurs à chaque frame d'énergie.

Ce fichier permet de calculer la fonction intermédiaire de diffusion :
```
F(q,t) = <n*(q,0) n(q,t)>
```

## Format du fichier

### Structure
```
# Header XVG standard
# Liste des 50 q-vecteurs générés
# q[ 0] = (  3.5565,   1.7783,  11.3224)  |q| =  12.0003
# q[ 1] = ( -3.5565,  -1.7783, -11.3224)  |q| =  12.0003
# ...
# q[49] = ( -8.8913,  -1.7783,   7.5490)  |q| =  11.7985

Time Re(q0) Im(q0) Re(q1) Im(q1) ... Re(q49) Im(q49)
0.0  68.71   2.94   68.71  -2.94  ...
0.1  43.88   31.78  43.88  -31.78 ...
```

### Colonnes
- Colonne 1 : Temps (ps)
- Colonnes 2-101 : 50 paires (partie réelle, partie imaginaire) de n(q,t)

## Utilisation avec Python

### Lecture basique
```python
import numpy as np

# Lire les données (ignorer les lignes de commentaire)
data = np.loadtxt('complex_obs.xvg', comments=['#', '@'])

time = data[:, 0]
real_parts = data[:, 1::2]  # Colonnes impaires : Re(q0), Re(q1), ...
imag_parts = data[:, 2::2]  # Colonnes paires : Im(q0), Im(q1), ...

# Construire n(q,t) complexe
n_qt = real_parts + 1j * imag_parts  # Shape: (n_frames, 50)
```

### Calculer F(q,t)
```python
# Pour chaque q-vecteur : F_i(q,t) = n*(q,0) × n(q,t)
nq = n_qt.shape[1]  # 50 q-vecteurs
F_qt_per_q = np.zeros((len(time), nq), dtype=complex)

for i in range(nq):
    n_q0_i = n_qt[0, i]  # n(q_i, t=0)
    F_qt_per_q[:, i] = np.conj(n_q0_i) * n_qt[:, i]

# Moyenner sur les 50 q-vecteurs
F_qt_avg = np.mean(F_qt_per_q, axis=1)

# Tracer |F(q,t)|
import matplotlib.pyplot as plt
plt.plot(time, np.abs(F_qt_avg))
plt.xlabel('Time (ps)')
plt.ylabel('|F(q,t)|')
plt.show()
```

## Continuation de run (Restart avec -append)

### Comportement

Lors d'un restart avec `-cpi` et `-append`, GROMACS :
1. ✅ Ouvre `complex_obs.xvg` en mode **append** (pas de réécriture)
2. ✅ N'écrit **pas** le header à nouveau
3. ⚠️ **Réécrit la frame du checkpoint** (comportement standard GROMACS)

Cela crée **un doublon de temps** au point de restart.

### Exemple
```bash
# Premier run : 0 → 0.5 ps
gmx mdrun -deffnm run -nsteps 500

# Continuation : 0.5 → 1.5 ps
gmx mdrun -deffnm run -cpi run.cpt -append -nsteps 1500
```

Le fichier `complex_obs.xvg` contiendra :
```
0.0
0.1
...
0.5   ← Écrit par le premier run
0.5   ← Réécrit par le restart (DOUBLON)
0.6
...
1.5
```

### Solution : Nettoyer les doublons

Utilisez le script fourni `clean_xvg_duplicates.py` :

```bash
# Nettoyer en créant un nouveau fichier
python clean_xvg_duplicates.py complex_obs.xvg complex_obs_clean.xvg

# Ou remplacer le fichier original
python clean_xvg_duplicates.py complex_obs.xvg
```

Le script :
- Garde la **première occurrence** de chaque temps
- Préserve le header complet
- Affiche les doublons supprimés

### Vérification manuelle

```bash
# Vérifier les doublons de temps
grep -v '^[#@]' complex_obs.xvg | awk '{print $1}' | sort | uniq -d

# Compter les frames uniques
grep -v '^[#@]' complex_obs.xvg | awk '{print $1}' | sort -u | wc -l
```

## Q-vecteurs générés

### Paramètres
- **q_value** = 12.0 rad/nm (modifiable dans `energyoutput.cpp`)
- **q_min** = 0.95 × q_value = 11.4 rad/nm
- **q_max** = 1.05 × q_value = 12.6 rad/nm
- **Nombre** : 50 q-vecteurs les plus proches de q_value

### Génération
Les q-vecteurs sont générés sur le **réseau réciproque** :
```
q = 2π × inv(box)^T × (nx, ny, nz)
```

Seuls les vecteurs avec |q| ∈ [q_min, q_max] sont gardés, puis les 50 **les plus proches de q_value** sont sélectionnés.

### Persistance
Les q-vecteurs sont :
- Générés **une seule fois** au premier appel
- Stockés en mémoire pendant toute la simulation
- **Réutilisés** en cas de restart (pas de régénération)

⚠️ **Important** : Si la boîte de simulation change significativement (NPT), les q-vecteurs restent fixes. Pour un calcul rigoureux, utilisez NVT ou régénérez les q-vecteurs.

## Fichiers fournis

- `test_read_xvg.py` : Exemple de lecture et calcul de F(q,t)
- `clean_xvg_duplicates.py` : Nettoyage des doublons après restart
- `test_continuation.sh` : Test automatique de continuation

## Références

Ce module implémente le calcul de :
- **Facteur de structure** : S(q) = |n(q)|²
- **Fonction intermédiaire** : F(q,t) = ⟨n*(q,0)n(q,t)⟩

Pour plus de détails sur la théorie, voir :
- Hansen & McDonald, "Theory of Simple Liquids"
- Frenkel & Smit, "Understanding Molecular Simulation"
