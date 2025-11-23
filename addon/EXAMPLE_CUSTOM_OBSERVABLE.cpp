/*
 * EXEMPLE : Ajouter une observable calculée uniquement lors de l'écriture .edr
 * 
 * Cet exemple montre comment ajouter une observable "CustomEnergy" qui est
 * calculée UNIQUEMENT lors de l'écriture dans le fichier .edr, sans accumulation
 * à chaque pas de temps (pour économiser du calcul).
 */

// ============================================================================
// PARTIE 1 : Modifications dans energyoutput.h
// ============================================================================

// Dans la classe EnergyOutput, ajoutez un membre privé pour stocker l'indice :
class EnergyOutput
{
private:
    // ... membres existants ...
    int iCustomEnergy_ = -1;  // Indice pour votre observable personnalisée
};

// Modifiez la signature de printStepToEnergyFile pour accepter les données nécessaires :
void printStepToEnergyFile(ener_file*              fp_ene,
                          bool                    bEne,
                          bool                    bDR,
                          bool                    bOR,
                          FILE*                   log,
                          int64_t                 step,
                          double                  time,
                          t_fcdata*               fcd,
                          gmx::Awh*               awh,
                          // NOUVEAUX PARAMÈTRES pour le calcul :
                          const rvec*             x,      // Positions atomiques
                          const rvec*             v,      // Vitesses (optionnel)
                          const gmx_mtop_t*       mtop,   // Topologie
                          const matrix            box);   // Boîte de simulation


// ============================================================================
// PARTIE 2 : Modifications dans le constructeur de EnergyOutput
// ============================================================================

EnergyOutput::EnergyOutput(/* ... paramètres existants ... */)
{
    // ... code existant d'initialisation ...
    
    // Ajoutez votre terme d'énergie personnalisé
    // Note: Ceci crée juste l'espace dans ebin_, mais ne sera rempli
    // que lors de l'écriture effective dans printStepToEnergyFile()
    if (true)  // Remplacez par votre condition (ex: inputrec.opts.custom_observable)
    {
        static const char* custom_energy_name = "Custom-Observable";
        iCustomEnergy_ = get_ebin_space(ebin_, 1, &custom_energy_name, unit_energy);
    }
    
    // ... reste du code du constructeur ...
}


// ============================================================================
// PARTIE 3 : Fonction de calcul personnalisée
// ============================================================================

namespace
{
/*! \brief Calcule votre observable personnalisée (lourd)
 *
 * Cette fonction n'est appelée QUE lors de l'écriture dans .edr,
 * pas à chaque pas de temps.
 *
 * \param[in] x     Positions atomiques
 * \param[in] v     Vitesses atomiques (peut être nullptr si non utilisé)
 * \param[in] mtop  Topologie moléculaire
 * \param[in] box   Boîte de simulation
 * \param[in] natoms Nombre d'atomes
 * \return La valeur de l'observable
 */
real calculateCustomObservable(const rvec*       x,
                              const rvec*       v,
                              const gmx_mtop_t* mtop,
                              const matrix      box,
                              int               natoms)
{
    real result = 0.0;
    
    // EXEMPLE 1 : Calcul basé sur les positions
    // (ex: moment dipolaire personnalisé, rayon de giration, etc.)
    for (int i = 0; i < natoms; i++)
    {
        // Votre calcul lourd ici
        // Exemple simple : somme des coordonnées x
        result += x[i][XX];
    }
    
    // EXEMPLE 2 : Calcul basé sur les paires d'atomes (O(N²) - très lourd!)
    // for (int i = 0; i < natoms - 1; i++)
    // {
    //     for (int j = i + 1; j < natoms; j++)
    //     {
    //         rvec dx;
    //         pbc_dx(/* pbc structure */, x[i], x[j], dx);
    //         real r2 = norm2(dx);
    //         result += 1.0 / sqrt(r2);  // Exemple : potentiel coulombien simplifié
    //     }
    // }
    
    // EXEMPLE 3 : Calcul basé sur les vitesses
    // if (v != nullptr)
    // {
    //     for (int i = 0; i < natoms; i++)
    //     {
    //         result += norm2(v[i]);  // Énergie cinétique alternative
    //     }
    // }
    
    return result;
}
} // namespace anonyme


// ============================================================================
// PARTIE 4 : Modification de printStepToEnergyFile
// ============================================================================

void EnergyOutput::printStepToEnergyFile(ener_file*        fp_ene,
                                         bool              bEne,
                                         bool              bDR,
                                         bool              bOR,
                                         FILE*             log,
                                         int64_t           step,
                                         double            time,
                                         t_fcdata*         fcd,
                                         gmx::Awh*         awh,
                                         const rvec*       x,
                                         const rvec*       v,
                                         const gmx_mtop_t* mtop,
                                         const matrix      box)
{
    t_enxframe fr;
    init_enxframe(&fr);
    fr.t       = time;
    fr.step    = step;
    fr.nsteps  = ebin_->nsteps;
    fr.dt      = delta_t_;
    fr.nsum    = ebin_->nsum;
    fr.nre     = (bEne) ? ebin_->nener : 0;
    fr.ener    = ebin_->e;
    
    // ... code existant pour disres, orires, etc. ...
    
    // ========================================================================
    // AJOUT : Calculer l'observable personnalisée UNIQUEMENT lors de l'écriture
    // ========================================================================
    
    if (bEne && iCustomEnergy_ >= 0 && x != nullptr && mtop != nullptr)
    {
        // CE CALCUL N'EST FAIT QUE TOUS LES nstenergy PAS !
        real customValue = calculateCustomObservable(x, v, mtop, box, mtop->natoms);
        
        // Stocker directement dans ebin_ pour cette frame
        // Note: on utilise false pour bSum car c'est une valeur instantanée
        add_ebin(ebin_, iCustomEnergy_, 1, &customValue, false);
    }
    
    // ... reste du code existant ...
    
    /* AWH bias blocks. */
    if (awh != nullptr)
    {
        awh->writeToEnergyFrame(step, &fr);
    }

    /* do the actual I/O */
    do_enx(fp_ene, &fr);
    
    if (fr.nre)
    {
        /* We have stored the sums, so reset the sum history */
        reset_ebin_sums(ebin_);
    }
    
    free_enxframe(&fr);
    
    // ... reste du code ...
}


// ============================================================================
// PARTIE 5 : Mise à jour des appels dans md.cpp
// ============================================================================

// Dans src/gromacs/mdrun/md.cpp, ligne ~1981
// Modifiez l'appel pour passer les nouvelles données :

if (do_log || do_ene || do_dr || do_or)
{
    energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                      do_ene,
                                      do_dr,
                                      do_or,
                                      do_log ? fpLog_ : nullptr,
                                      step,
                                      t,
                                      fr_->fcdata.get(),
                                      awh.get(),
                                      // NOUVEAUX ARGUMENTS :
                                      state_->x.rvec_array(),    // Positions
                                      state_->v.rvec_array(),    // Vitesses
                                      top_global,                 // Topologie
                                      state_->box);              // Boîte
}


// ============================================================================
// PARTIE 6 : Mise à jour des appels dans minimize.cpp
// ============================================================================

// Dans src/gromacs/mdrun/minimize.cpp, ligne ~1432
// Même modification pour les minimiseurs :

energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  do_x,
                                  do_x,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  // NOUVEAUX ARGUMENTS :
                                  state_global->x.rvec_array(),
                                  nullptr,              // Pas de vitesses en minimisation
                                  top_global,
                                  state_global->box);


// ============================================================================
// PARTIE 7 : Mise à jour des appels dans rerun.cpp
// ============================================================================

// Dans src/gromacs/mdrun/rerun.cpp, ligne ~826

energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  bDisRe,
                                  bOriRe,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  // NOUVEAUX ARGUMENTS :
                                  state.x.rvec_array(),
                                  state.v.rvec_array(),
                                  top_global,
                                  state.box);


// ============================================================================
// NOTES IMPORTANTES
// ============================================================================

/*
 * PERFORMANCE :
 * - Le calcul dans calculateCustomObservable() n'est fait que tous les
 *   nstenergy pas (défini dans le fichier .mdp)
 * - Exemple: si nstenergy=100, le calcul se fait 100x moins souvent
 *            qu'un calcul dans addDataAtEnergyStep()
 *
 * MOYENNES :
 * - Si vous utilisez add_ebin() avec bSum=false, vous obtenez la valeur
 *   instantanée à chaque écriture
 * - Si vous utilisez add_ebin() avec bSum=true, GROMACS calculera
 *   automatiquement la moyenne sur les pas intermédiaires (mais il faudrait
 *   alors appeler add_ebin() à chaque pas dans addDataAtEnergyStep)
 *
 * LECTURE DU .edr :
 * - Votre observable apparaîtra dans le fichier .edr
 * - Utilisable avec "gmx energy" : tapez le nom "Custom-Observable"
 * - Accessible via l'API Python de GROMACS
 *
 * TESTS :
 * - Ajoutez des tests dans src/gromacs/mdlib/tests/energyoutput.cpp
 * - Vérifiez que la valeur est correctement écrite et lue
 */
