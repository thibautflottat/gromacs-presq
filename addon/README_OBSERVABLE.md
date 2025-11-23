# RÉSUMÉ : Ajouter une Observable dans les fichiers .edr GROMACS

## Question posée
> "Où et comment sont générés les .edr ? Comment rajouter une observable dans ces fichiers ?"
> "J'aimerai calculer l'observable uniquement pour l'écriture, sans accumuler, parce que c'est un calcul lourd."

## Réponse courte

Les fichiers `.edr` sont générés par la classe `EnergyOutput` dans `src/gromacs/mdlib/energyoutput.cpp`.

Pour ajouter une observable **calculée uniquement lors de l'écriture** (optimisation pour calculs lourds) :

1. **Modifier `energyoutput.h`** : Ajouter un membre `int iMyObservable_` et modifier la signature de `printStepToEnergyFile()`
2. **Modifier le constructeur** : Réserver l'espace avec `get_ebin_space()`
3. **Calculer dans `printStepToEnergyFile()`** : Juste avant `do_enx()`, calculer et appeler `add_ebin()`
4. **Mettre à jour les appels** : Dans `md.cpp`, `minimize.cpp`, `rerun.cpp`, `mimic.cpp`

## Fichiers créés pour vous

J'ai créé 4 fichiers dans `/Users/thibaut/dev/gromacs/gromacs-2024.2/` :

### 1. `CUSTOM_OBSERVABLE_GUIDE.md`
Guide complet avec :
- Architecture du système .edr
- Deux approches possibles (avec comparaison)
- Avantages/inconvénients
- Exemple du rayon de giration

### 2. `EXAMPLE_CUSTOM_OBSERVABLE.cpp`
Code exemple annoté montrant :
- Modifications dans `energyoutput.h`
- Modifications dans le constructeur
- Fonction de calcul personnalisée
- Modifications dans `printStepToEnergyFile()`
- Mises à jour des appels dans tous les fichiers mdrun

### 3. `FLUX_OBSERVABLE.md`
Diagrammes visuels montrant :
- Comparaison approche classique vs optimisée
- Diagramme de flux détaillé
- Code simplifié annoté
- Calcul de performance

### 4. `PATCH_RADIUS_GYRATION.cpp`
Patch complet et fonctionnel pour :
- Ajouter le rayon de giration comme observable
- Code prêt à copier-coller
- Tous les fichiers modifiés avec numéros de lignes
- Instructions de compilation et test

## Réponse détaillée

### Architecture des fichiers .edr

```
EnergyOutput (energyoutput.cpp)
    ├─ Constructeur : Initialise les termes d'énergie
    ├─ addDataAtEnergyStep() : Accumule les données (appelé souvent)
    └─ printStepToEnergyFile() : Écrit dans .edr (appelé rarement)
           └─ do_enx() : Écriture binaire réelle
```

### Système d'accumulation (ebin_)

GROMACS utilise un système de "bins" d'énergie (`t_ebin`) qui :
- Accumule les valeurs à chaque pas
- Calcule les moyennes
- Les écrit périodiquement dans le .edr

### Votre besoin : Calcul lourd sans accumulation

**Problème** : Si vous utilisez le flux normal (`addDataAtEnergyStep()`), votre calcul lourd sera exécuté à chaque pas d'énergie.

**Solution** : Calculer directement dans `printStepToEnergyFile()`, qui n'est appelé que lors de l'écriture.

### Code minimal

```cpp
// energyoutput.h
class EnergyOutput {
private:
    int iMyObservable_ = -1;
public:
    void printStepToEnergyFile(..., const rvec* x, const gmx_mtop_t* mtop, ...);
};

// energyoutput.cpp - Constructeur
iMyObservable_ = get_ebin_space(ebin_, 1, "My-Observable", "unit");

// energyoutput.cpp - printStepToEnergyFile, AVANT do_enx()
if (bEne && iMyObservable_ >= 0 && x != nullptr) {
    real value = calculateMyObservable(x, mtop);
    add_ebin(ebin_, iMyObservable_, 1, &value, false);
}

// md.cpp
energyOutput.printStepToEnergyFile(..., state_->x.rvec_array(), top_global, ...);
```

### Fichiers à modifier

1. **`src/gromacs/mdlib/energyoutput.h`** (2 modifications)
   - Ajouter membre privé `iMyObservable_`
   - Modifier signature de `printStepToEnergyFile()`

2. **`src/gromacs/mdlib/energyoutput.cpp`** (3 modifications)
   - Constructeur : `get_ebin_space()`
   - Fonction de calcul : `calculateMyObservable()`
   - `printStepToEnergyFile()` : Calcul et `add_ebin()`

3. **`src/gromacs/mdrun/md.cpp`** (1 modification)
   - Appel de `printStepToEnergyFile()` avec nouveaux paramètres

4. **`src/gromacs/mdrun/minimize.cpp`** (~10 modifications)
   - Tous les appels de `printStepToEnergyFile()`

5. **`src/gromacs/mdrun/rerun.cpp`** (1 modification)
   - Appel de `printStepToEnergyFile()`

6. **`src/gromacs/mdrun/mimic.cpp`** (1 modification)
   - Appel de `printStepToEnergyFile()`

### Performance attendue

Si `nstenergy = 100` (écriture tous les 100 pas) :

```
Approche classique (addDataAtEnergyStep) :
- 100 calculs pour 1 écriture
- Temps : T

Approche optimisée (printStepToEnergyFile) :
- 1 calcul pour 1 écriture
- Temps : T/100
- Gain : 100x
```

### Test de l'implémentation

```bash
# 1. Recompiler
cd build
cmake ..
make -j8
make install

# 2. Exécuter une simulation
gmx mdrun -s topol.tpr -deffnm test

# 3. Vérifier que l'observable est présente
gmx energy -f test.edr
# → Cherchez "My-Observable" ou le nom que vous avez choisi

# 4. Extraire les valeurs
echo "My-Observable" | gmx energy -f test.edr -o my_obs.xvg

# 5. Visualiser
xmgrace my_obs.xvg
```

### Alternatives si calcul TRÈS lourd

Si même avec `nstenergy=100` c'est trop lourd :

1. **Option A** : Augmenter `nstenergy` dans le .mdp
   ```mdp
   nstenergy = 1000  ; Écriture tous les 1000 pas
   ```

2. **Option B** : Ajouter un flag dans le .mdp pour activer/désactiver
   ```cpp
   if (inputrec.opts.calculate_my_observable) {
       real value = calculateMyObservable(...);
       add_ebin(..., value, ...);
   }
   ```

3. **Option C** : Post-traitement externe
   - Écrire positions dans .xtc (léger)
   - Calculer l'observable après coup avec un script Python

## Prochaines étapes recommandées

1. **Lire** `FLUX_OBSERVABLE.md` pour comprendre le flux complet
2. **Étudier** `PATCH_RADIUS_GYRATION.cpp` comme exemple concret
3. **Adapter** le code à votre observable spécifique
4. **Compiler** et tester sur un petit système d'abord
5. **Vérifier** que les valeurs sont cohérentes
6. **Mesurer** le gain de performance

## Questions fréquentes

**Q : L'observable apparaît-elle automatiquement dans gmx energy ?**
R : Oui, une fois dans le .edr, elle est accessible via `gmx energy`.

**Q : Puis-je avoir des moyennes ?**
R : Avec `add_ebin(..., false)`, ce sont des valeurs instantanées. Pour des moyennes, utilisez `gmx analyze` après coup.

**Q : Cela marche-t-il avec MPI ?**
R : Oui, mais assurez-vous que tous les rangs ont les mêmes données ou que seul le rank 0 écrit.

**Q : Puis-je ajouter plusieurs observables ?**
R : Oui, répétez le processus pour chaque observable.

**Q : Cela affecte-t-il les checkpoints (.cpt) ?**
R : Non, les valeurs instantanées ne sont pas dans les checkpoints.

## Ressources

- Code source : `src/gromacs/mdlib/energyoutput.cpp`
- Tests : `src/gromacs/mdlib/tests/energyoutput.cpp`
- Format .edr : `src/gromacs/fileio/enxio.cpp`
- Documentation : `docs/` (après compilation)

## Contact et support

Pour des questions spécifiques à votre observable :
1. Consultez d'abord les 4 fichiers créés
2. Regardez les exemples dans le code (DISRES, ORIRES, etc.)
3. Utilisez le forum GROMACS : https://gromacs.bioexcel.eu

---

**Bon courage avec votre implémentation !** 🚀
