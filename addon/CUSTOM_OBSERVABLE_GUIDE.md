# Guide : Ajouter une Observable Calculée Uniquement lors de l'Écriture .edr

Ce guide explique comment ajouter une observable lourde qui n'est calculée que lors de l'écriture dans le fichier .edr, **sans accumulation** à chaque pas de temps.

## Architecture du système

### Flux actuel
```
Boucle MD (chaque pas) :
  ├─ addDataAtEnergyStep()    ← Accumule les énergies
  └─ printStepToEnergyFile()  ← Écrit dans .edr (tous les nstenergy pas)
        └─ do_enx()
```

### Deux approches possibles

## Approche 1 : Calcul dans printStepToEnergyFile() [RECOMMANDÉ]

Cette approche calcule l'observable **uniquement** lors de l'écriture, sans passer par le système d'accumulation.

### Étape 1 : Modifier le header energyoutput.h

Ajoutez les paramètres nécessaires à la signature de `printStepToEnergyFile()` :

```cpp
// Dans src/gromacs/mdlib/energyoutput.h, ligne ~215

void printStepToEnergyFile(ener_file*        fp_ene,
                          bool              bEne,
                          bool              bDR,
                          bool              bOR,
                          FILE*             log,
                          int64_t           step,
                          double            time,
                          t_fcdata*         fcd,
                          gmx::Awh*         awh,
                          // AJOUTEZ VOS PARAMÈTRES ICI :
                          const rvec*       x,        // Positions atomiques
                          const gmx_mtop_t* mtop);    // Topologie (si nécessaire)
```

### Étape 2 : Implémenter le calcul dans printStepToEnergyFile()

```cpp
// Dans src/gromacs/mdlib/energyoutput.cpp, dans printStepToEnergyFile()
// Juste AVANT l'appel à do_enx(fp_ene, &fr); (ligne ~1273)

void EnergyOutput::printStepToEnergyFile(ener_file* fp_ene,
                                         bool       bEne,
                                         bool       bDR,
                                         bool       bOR,
                                         FILE*      log,
                                         int64_t    step,
                                         double     time,
                                         t_fcdata*  fcd,
                                         gmx::Awh*  awh,
                                         const rvec* x,
                                         const gmx_mtop_t* mtop)
{
    // ... code existant ...
    
    // AJOUTEZ VOTRE CALCUL ICI (juste avant do_enx) :
    
    // 1. Calculez votre observable lourde
    real myHeavyObservable = 0.0;
    if (bEne)  // Seulement si on écrit l'énergie
    {
        // Votre calcul lourd ici
        myHeavyObservable = calculateMyHeavyObservable(x, mtop);
        
        // 2. Ajoutez un bloc personnalisé au frame
        int customBlockIndex = fr.nblock;
        fr.nblock += 1;
        add_blocks_enxframe(&fr, fr.nblock);
        
        add_subblocks_enxblock(&(fr.block[customBlockIndex]), 1);
        fr.block[customBlockIndex].id = enxNR + 1; // ID unique
        fr.block[customBlockIndex].sub[0].nr = 1;  // 1 valeur
        
        // Allouez et stockez la valeur
        real* customValue = nullptr;
        snew(customValue, 1);
        customValue[0] = myHeavyObservable;
        
#if !GMX_DOUBLE
        fr.block[customBlockIndex].sub[0].type = XdrDataType::Float;
        fr.block[customBlockIndex].sub[0].fval = customValue;
#else
        fr.block[customBlockIndex].sub[0].type = XdrDataType::Double;
        fr.block[customBlockIndex].sub[0].dval = customValue;
#endif
    }
    
    // do the actual I/O
    do_enx(fp_ene, &fr);
    
    // ... reste du code ...
}
```

### Étape 3 : Fonction de calcul personnalisée

```cpp
// Ajoutez votre fonction de calcul (peut être dans un fichier séparé)

static real calculateMyHeavyObservable(const rvec* x, const gmx_mtop_t* mtop)
{
    real result = 0.0;
    
    // Exemple : calcul d'une propriété géométrique complexe
    // Ceci n'est exécuté QUE lors de l'écriture dans .edr
    
    for (int i = 0; i < mtop->natoms; i++)
    {
        // Votre calcul lourd ici
        // Par exemple : distances, angles, propriétés géométriques, etc.
        result += expensiveCalculation(x[i]);
    }
    
    return result;
}
```

### Étape 4 : Mettre à jour les appels dans md.cpp, minimize.cpp, etc.

```cpp
// Dans src/gromacs/mdrun/md.cpp, ligne ~1981

energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  do_ene,
                                  do_dr,
                                  do_or,
                                  do_log ? fpLog_ : nullptr,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  awh.get(),
                                  state_->x.rvec_array(),  // AJOUTÉ
                                  top_global);             // AJOUTÉ
```

## Approche 2 : Observable avec Accumulation Conditionnelle

Si vous voulez quand même utiliser le système ebin mais éviter le calcul à chaque pas :

### Dans addDataAtEnergyStep()

```cpp
void EnergyOutput::addDataAtEnergyStep(bool bDoDHDL,
                                      bool bSum,  // ← KEY: false quand on ne somme pas
                                      ...)
{
    // ... code existant ...
    
    // Calculez uniquement si c'est un pas d'écriture (bSum == true pour nstenergy)
    if (bSum && bEne)  // Seulement aux pas d'écriture
    {
        real heavyValue = calculateMyHeavyObservable(/* params */);
        add_ebin(ebin_, iMyObservable_, 1, &heavyValue, bSum);
    }
}
```

### Dans md.cpp

```cpp
// Le paramètre bCalcEnerStep contrôle si on somme ou pas
energyOutput.addDataAtEnergyStep(outputDHDL,
                                bCalcEnerStep,  // ← true seulement aux pas nstenergy
                                t,
                                md->tmass,
                                ...);
```

## Avantages et Inconvénients

### Approche 1 (Calcul direct dans printStepToEnergyFile)
✅ **Avantages** :
- Calcul SEULEMENT lors de l'écriture
- Pas d'accumulation inutile
- Performance optimale
- Contrôle total sur le format

❌ **Inconvénients** :
- Pas de moyennes automatiques
- Format de bloc personnalisé
- Plus de code à écrire

### Approche 2 (Accumulation conditionnelle)
✅ **Avantages** :
- Utilise l'infrastructure ebin existante
- Moyennes automatiques (si plusieurs frames entre écritures)
- Format standard

❌ **Inconvénients** :
- Encore appelé à chaque pas (avec une condition)
- Moins flexible

## Exemple Complet : Rayon de Giration

```cpp
// energyoutput.cpp - Dans le constructeur
if (/* condition pour activer cette observable */)
{
    iRadiusOfGyration_ = get_ebin_space(ebin_, 1, "Radius-Gyration", "nm");
}

// energyoutput.cpp - Dans printStepToEnergyFile, AVANT do_enx()
if (bEne && x != nullptr && mtop != nullptr)
{
    // Calcul du rayon de giration (lourd pour gros systèmes)
    real Rg = 0.0;
    rvec centerOfMass = {0, 0, 0};
    real totalMass = 0.0;
    
    // Centre de masse
    for (int i = 0; i < mtop->natoms; i++)
    {
        real mass = mtop->atomtypes.atomNumberFromAtomType(i); // Simplifié
        centerOfMass[XX] += mass * x[i][XX];
        centerOfMass[YY] += mass * x[i][YY];
        centerOfMass[ZZ] += mass * x[i][ZZ];
        totalMass += mass;
    }
    svmul(1.0/totalMass, centerOfMass, centerOfMass);
    
    // Rayon de giration
    for (int i = 0; i < mtop->natoms; i++)
    {
        real mass = mtop->atomtypes.atomNumberFromAtomType(i);
        rvec dr;
        rvec_sub(x[i], centerOfMass, dr);
        Rg += mass * norm2(dr);
    }
    Rg = std::sqrt(Rg / totalMass);
    
    // Ajouter au frame (via un bloc personnalisé)
    int rgBlock = fr.nblock;
    fr.nblock += 1;
    add_blocks_enxframe(&fr, fr.nblock);
    // ... configuration du bloc comme montré plus haut ...
}
```

## Résumé

Pour une observable **vraiment lourde** :
1. **Recommandé** : Approche 1 - Calcul dans `printStepToEnergyFile()`
2. Ajoutez les paramètres nécessaires (positions, topologie, etc.)
3. Calculez juste avant `do_enx()`
4. Créez un bloc personnalisé dans `t_enxframe`
5. Mettez à jour tous les appels dans les fichiers mdrun

**La clé** : `printStepToEnergyFile()` n'est appelée que tous les `nstenergy` pas, donc votre calcul lourd ne s'exécutera que lors de l'écriture réelle dans le fichier .edr.
