# Flux de Calcul : Observable sans Accumulation

## Comparaison des deux approches

### ❌ APPROCHE CLASSIQUE (avec accumulation - LOURD)

```
Pas 0    : addDataAtEnergyStep() → CALCUL LOURD → accumule
Pas 1    : addDataAtEnergyStep() → CALCUL LOURD → accumule
Pas 2    : addDataAtEnergyStep() → CALCUL LOURD → accumule
...
Pas 99   : addDataAtEnergyStep() → CALCUL LOURD → accumule
Pas 100  : addDataAtEnergyStep() → CALCUL LOURD → accumule
           printStepToEnergyFile() → écrit moyenne → .edr
           
Pas 101  : addDataAtEnergyStep() → CALCUL LOURD → accumule
...

TOTAL : 100 calculs lourds pour 1 écriture !
```

### ✅ APPROCHE OPTIMISÉE (calcul à l'écriture - LÉGER)

```
Pas 0    : addDataAtEnergyStep() → (rien pour notre observable)
Pas 1    : addDataAtEnergyStep() → (rien)
Pas 2    : addDataAtEnergyStep() → (rien)
...
Pas 99   : addDataAtEnergyStep() → (rien)
Pas 100  : addDataAtEnergyStep() → (rien)
           printStepToEnergyFile() → CALCUL LOURD 1 FOIS → .edr
           
Pas 101  : addDataAtEnergyStep() → (rien)
...
Pas 200  : addDataAtEnergyStep() → (rien)
           printStepToEnergyFile() → CALCUL LOURD 1 FOIS → .edr

TOTAL : 1 calcul lourd pour 1 écriture !
Gain : 100x plus rapide (si nstenergy=100)
```

## Diagramme de flux détaillé

```
╔═══════════════════════════════════════════════════════════════╗
║                    BOUCLE MD PRINCIPALE                        ║
╚═══════════════════════════════════════════════════════════════╝
                             │
                             ▼
        ┌────────────────────────────────────┐
        │  Intégration (positions, vitesses) │
        └────────────────────────────────────┘
                             │
                             ▼
        ┌────────────────────────────────────┐
        │   Calcul des forces                │
        └────────────────────────────────────┘
                             │
                             ▼
        ┌────────────────────────────────────┐
        │   Calcul des énergies standard     │
        │   (potentielle, cinétique, etc.)   │
        └────────────────────────────────────┘
                             │
                             ▼
                  ┌──────────────────┐
                  │ bCalcEner ?      │
                  └──────────────────┘
                    │              │
                  OUI             NON
                    │              │
                    ▼              ▼
    ┌───────────────────────┐   ┌─────────────────────┐
    │ addDataAtEnergyStep() │   │ recordNonEnergyStep │
    │                       │   └─────────────────────┘
    │ - Accumule énergies   │
    │ - PAS notre observable│   (juste incrémente compteur)
    └───────────────────────┘
                    │
                    ▼
          ┌──────────────────┐
          │ do_ene ?         │
          │ (step % nstener) │
          └──────────────────┘
            │              │
          OUI             NON
            │              │
            ▼              └─────► Continue boucle
┌───────────────────────────────┐
│  printStepToEnergyFile()      │
│                               │
│  1. Prépare le frame          │
│                               │
│  2. ★ CALCUL OBSERVABLE ★    │◄─── ICI SEULEMENT !
│     calculateCustomObservable()│
│     (1 seul appel tous les    │
│      nstenergy pas)           │
│                               │
│  3. add_ebin() instantané     │
│                               │
│  4. do_enx() → écrit .edr     │
│                               │
└───────────────────────────────┘
            │
            ▼
    Continue boucle MD
```

## Code simplifié annoté

### Fichier md.cpp (boucle principale)

```cpp
// Boucle principale de simulation
do {
    step++;
    
    // === Intégration, forces, etc. ===
    do_force(...);
    update_coords(...);
    
    // === Gestion de l'énergie ===
    bool bCalcEner = do_per_step(step, inputrec->nstenergy);
    bool do_ene    = do_per_step(step, inputrec->nstenergy);
    
    if (bCalcEner) {
        // Accumule les énergies STANDARD (rapides)
        energyOutput.addDataAtEnergyStep(
            outputDHDL,
            bCalcEner,  // true tous les nstenergy pas
            t,
            md->tmass,
            enerd_,
            ir->fepvals.get(),
            lastbox,
            ptCouplingArrays,
            state_->fep_state,
            total_vir,
            pres,
            ekind_,
            mu_tot,
            constr_
        );
        // Note: NOTRE observable N'est PAS calculée ici !
    }
    
    if (do_ene) {
        // Écriture dans .edr
        energyOutput.printStepToEnergyFile(
            mdoutf_get_fp_ene(outf),
            do_ene,
            do_dr,
            do_or,
            do_log ? fpLog_ : nullptr,
            step,
            t,
            fr_->fcdata.get(),
            awh.get(),
            // === PARAMÈTRES AJOUTÉS ===
            state_->x.rvec_array(),  // Positions pour le calcul
            state_->v.rvec_array(),  // Vitesses (si nécessaire)
            top_global,              // Topologie
            state_->box              // Boîte
        );
        // C'est ICI que notre observable est calculée !
    }
    
} while (step < nsteps);
```

### Fichier energyoutput.cpp

```cpp
void EnergyOutput::printStepToEnergyFile(
    ener_file* fp_ene,
    bool bEne,
    bool bDR,
    bool bOR,
    FILE* log,
    int64_t step,
    double time,
    t_fcdata* fcd,
    gmx::Awh* awh,
    const rvec* x,        // ← Positions
    const rvec* v,        // ← Vitesses
    const gmx_mtop_t* mtop, // ← Topologie
    const matrix box)     // ← Boîte
{
    // 1. Préparer la frame d'énergie
    t_enxframe fr;
    init_enxframe(&fr);
    fr.nre = (bEne) ? ebin_->nener : 0;
    fr.ener = ebin_->e;
    
    // 2. Traiter les blocs standards (disres, orires, etc.)
    // ... code existant ...
    
    // 3. ★★★ CALCUL DE NOTRE OBSERVABLE ★★★
    if (bEne && iCustomEnergy_ >= 0 && x != nullptr) {
        
        // CE CALCUL LOURD N'EST FAIT QUE MAINTENANT !
        // Pas à chaque pas, seulement tous les nstenergy pas
        real customValue = calculateCustomObservable(
            x, v, mtop, box, mtop->natoms
        );
        
        // Ajouter à ebin pour cette frame
        // false = pas d'accumulation, valeur instantanée
        add_ebin(ebin_, iCustomEnergy_, 1, &customValue, false);
    }
    
    // 4. Blocs free energy, AWH, etc.
    if (dhc_) {
        mde_delta_h_coll_handle_block(dhc_.get(), &fr, fr.nblock);
    }
    if (awh != nullptr) {
        awh->writeToEnergyFrame(step, &fr);
    }
    
    // 5. ★★★ ÉCRITURE DANS .EDR ★★★
    do_enx(fp_ene, &fr);
    
    // 6. Réinitialiser les sommes
    if (fr.nre) {
        reset_ebin_sums(ebin_);
    }
    
    free_enxframe(&fr);
}
```

## Paramètres .mdp pertinents

```mdp
; Contrôle la fréquence d'écriture dans .edr
; Votre calcul lourd sera fait tous les nstenergy pas
nstenergy       = 100     ; Écriture tous les 100 pas
                          ; → Calcul lourd 10x par ns (si dt=0.002)
                          ; au lieu de 500x par ns !

; Pour référence
nsteps          = 50000   ; 100 ps avec dt=0.002
dt              = 0.002   ; 2 fs

; Résultat:
; - 50000 pas de simulation
; - 500 écritures dans .edr (tous les 100 pas)
; - 500 calculs de votre observable (au lieu de 50000 !)
; - Gain de performance : 100x
```

## Vérification du résultat

```bash
# Lire le fichier .edr
gmx energy -f ener.edr

# Sélectionner votre observable
# Tapez le numéro ou le nom : "Custom-Observable"

# Visualiser
xmgrace energy.xvg
```

## Points clés à retenir

1. **addDataAtEnergyStep()** est appelé à chaque pas d'énergie
   → Ne PAS y faire de calculs lourds

2. **printStepToEnergyFile()** est appelé seulement lors de l'écriture
   → Parfait pour calculs lourds

3. **nstenergy** contrôle la fréquence
   → Plus grand = moins de calculs = plus rapide

4. **Pas de moyennes automatiques** avec cette approche
   → Valeurs instantanées dans le .edr
   → Vous pouvez moyenner après coup avec gmx analyze

5. **Accès aux données nécessaires**
   → Positions (x), vitesses (v), topologie (mtop), boîte (box)
   → Passés en paramètres depuis md.cpp

## Performance

Pour un système de 100000 atomes, nstenergy=100, calcul O(N²) :

```
AVANT (calcul à chaque pas) :
- 100000 atomes × 100000 atomes × 50000 pas = catastrophe
- Temps estimé : plusieurs jours

APRÈS (calcul à l'écriture) :
- 100000 atomes × 100000 atomes × 500 écritures = gérable
- Temps estimé : quelques minutes
- Gain : ~100x
```
