/*
 * PATCH COMPLET : Ajouter le Rayon de Giration comme Observable
 * 
 * Cet exemple montre un patch complet et fonctionnel pour ajouter
 * le calcul du rayon de giration (calcul lourd pour gros systèmes)
 * uniquement lors de l'écriture dans .edr
 */

// ============================================================================
// FICHIER 1 : src/gromacs/mdlib/energyoutput.h
// ============================================================================

// MODIFICATION 1A : Ajouter le membre privé (ligne ~420)
class EnergyOutput
{
private:
    // ... membres existants ...
    
    //! Index for radius of gyration in ebin
    int iRadiusGyration_ = -1;
    
    // ... autres membres ...
};

// MODIFICATION 1B : Modifier la signature de printStepToEnergyFile (ligne ~215)
/*! \brief Writes current quantites to log and energy files.
 *
 * Prints current values of energies, pressure, temperature, restraint
 * data, etc. to energy output file and to the log file (if not nullptr).
 */
void printStepToEnergyFile(ener_file*              fp_ene,
                          bool                    bEne,
                          bool                    bDR,
                          bool                    bOR,
                          FILE*                   log,
                          int64_t                 step,
                          double                  time,
                          t_fcdata*               fcd,
                          gmx::Awh*               awh,
                          const rvec*             x,
                          const gmx_mtop_t*       mtop,
                          const matrix            box);


// ============================================================================
// FICHIER 2 : src/gromacs/mdlib/energyoutput.cpp
// ============================================================================

// MODIFICATION 2A : Dans le constructeur EnergyOutput::EnergyOutput (ligne ~545)

EnergyOutput::EnergyOutput(ener_file*                fp_ene,
                           const gmx_mtop_t&         mtop,
                           const t_inputrec&         inputrec,
                           const pull_t*             pull_work,
                           FILE*                     fp_dhdl,
                           bool                      isRerun,
                           const StartingBehavior    startingBehavior,
                           const bool                simulationsShareState,
                           const MDModulesNotifiers& mdModulesNotifiers) :
    haveFepLambdaMoves_(haveFepLambdaMoves(inputrec))
{
    // ... code existant ...
    
    // AJOUT : Juste avant la fin du constructeur (ligne ~543)
    
    // Add radius of gyration to energy output
    // This will be computed only when writing to .edr file
    if (!isRerun)  // Ne pas calculer en rerun
    {
        static const char* rg_name = "Radius-Gyration";
        iRadiusGyration_ = get_ebin_space(ebin_, 1, &rg_name, "nm");
    }
    
    // ... reste du code du constructeur (do_enxnms, etc.) ...
}


// MODIFICATION 2B : Fonction de calcul du rayon de giration (à ajouter avant printStepToEnergyFile)

namespace
{
/*! \brief Calculate radius of gyration of the system
 *
 * This is a computationally expensive operation that should only
 * be called when writing energy output.
 *
 * \param[in] x       Atomic coordinates
 * \param[in] mtop    Molecular topology
 * \param[in] box     Simulation box
 * \return Radius of gyration in nm
 */
real calculateRadiusOfGyration(const rvec* x, const gmx_mtop_t* mtop, const matrix box)
{
    if (x == nullptr || mtop == nullptr)
    {
        return 0.0;
    }
    
    const int natoms = mtop->natoms;
    
    // Calculate center of mass
    rvec centerOfMass = { 0.0, 0.0, 0.0 };
    real totalMass    = 0.0;
    
    int globalAtomIndex = 0;
    for (const gmx_molblock_t& molblock : mtop->molblock)
    {
        const gmx_moltype_t& moltype = mtop->moltype[molblock.type];
        const int            natomsPerMol = moltype.atoms.nr;
        
        for (int mol = 0; mol < molblock.nmol; mol++)
        {
            for (int a = 0; a < natomsPerMol; a++)
            {
                const real mass = moltype.atoms.atom[a].m;
                
                centerOfMass[XX] += mass * x[globalAtomIndex][XX];
                centerOfMass[YY] += mass * x[globalAtomIndex][YY];
                centerOfMass[ZZ] += mass * x[globalAtomIndex][ZZ];
                totalMass += mass;
                
                globalAtomIndex++;
            }
        }
    }
    
    if (totalMass <= 0.0)
    {
        return 0.0;
    }
    
    svmul(1.0 / totalMass, centerOfMass, centerOfMass);
    
    // Calculate radius of gyration
    real radiusSquared = 0.0;
    
    globalAtomIndex = 0;
    for (const gmx_molblock_t& molblock : mtop->molblock)
    {
        const gmx_moltype_t& moltype = mtop->moltype[molblock.type];
        const int            natomsPerMol = moltype.atoms.nr;
        
        for (int mol = 0; mol < molblock.nmol; mol++)
        {
            for (int a = 0; a < natomsPerMol; a++)
            {
                const real mass = moltype.atoms.atom[a].m;
                
                rvec dr;
                // Consider periodic boundary conditions if needed
                pbc_dx_aiuc(nullptr, x[globalAtomIndex], centerOfMass, dr);
                
                radiusSquared += mass * norm2(dr);
                
                globalAtomIndex++;
            }
        }
    }
    
    return std::sqrt(radiusSquared / totalMass);
}
} // namespace


// MODIFICATION 2C : Dans printStepToEnergyFile (ligne ~1153, AVANT do_enx())

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
    
    /* whether we are going to write anything out: */
    if (fr.nre || ndisre || nr[enxOR] || nr[enxORI])
    {
        // ... code existant de préparation des blocs ...
        
        // === AJOUT : Calcul du rayon de giration ===
        // Ce code s'exécute UNIQUEMENT lors de l'écriture dans .edr
        if (bEne && iRadiusGyration_ >= 0 && x != nullptr && mtop != nullptr)
        {
            // Calculate radius of gyration (expensive operation)
            real radiusGyration = calculateRadiusOfGyration(x, mtop, box);
            
            // Add to energy frame (false = instantaneous value, no averaging)
            add_ebin(ebin_, iRadiusGyration_, 1, &radiusGyration, false);
        }
        
        /* Free energy perturbation blocks */
        if (dhc_)
        {
            mde_delta_h_coll_handle_block(dhc_.get(), &fr, fr.nblock);
        }

        /* we can now free & reset the data in the blocks */
        if (dhc_)
        {
            mde_delta_h_coll_reset(dhc_.get());
        }

        /* AWH bias blocks. */
        if (awh != nullptr)
        {
            awh->writeToEnergyFrame(step, &fr);
        }

        /* do the actual I/O */
        do_enx(fp_ene, &fr);
        
        // ... reste du code ...
    }
    
    // ... reste de la fonction ...
}


// ============================================================================
// FICHIER 3 : src/gromacs/mdrun/md.cpp
// ============================================================================

// MODIFICATION 3 : Appel de printStepToEnergyFile (ligne ~1981)

// Dans la boucle principale do_md(), section d'écriture de l'énergie :

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
                                      // === NOUVEAUX ARGUMENTS ===
                                      state_->x.rvec_array(),  // Positions
                                      top_global,              // Topologie
                                      state_->box);            // Boîte
}


// ============================================================================
// FICHIER 4 : src/gromacs/mdrun/minimize.cpp
// ============================================================================

// MODIFICATION 4A : Premier appel (ligne ~1432)
energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  do_x,
                                  do_x,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  // === NOUVEAUX ARGUMENTS ===
                                  state_global->x.rvec_array(),
                                  top_global,
                                  state_global->box);

// MODIFICATION 4B : Deuxième appel (ligne ~1893)
energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  do_x,
                                  do_x,
                                  do_log ? fpLog_ : nullptr,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  state_global->x.rvec_array(),
                                  top_global,
                                  state_global->box);

// MODIFICATION 4C : Troisième appel (ligne ~1943)
energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  false,
                                  false,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  state_global->x.rvec_array(),
                                  top_global,
                                  state_global->box);

// ... et ainsi de suite pour tous les autres appels dans minimize.cpp


// ============================================================================
// FICHIER 5 : src/gromacs/mdrun/rerun.cpp
// ============================================================================

// MODIFICATION 5 : Appel dans do_rerun (ligne ~826)
energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  bDisRe,
                                  bOriRe,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  state.x.rvec_array(),
                                  top_global,
                                  state.box);


// ============================================================================
// FICHIER 6 : src/gromacs/mdrun/mimic.cpp
// ============================================================================

// MODIFICATION 6 : Appel dans do_mimic (ligne ~743)
energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf),
                                  true,
                                  bDisRe,
                                  bOriRe,
                                  fpLog_,
                                  step,
                                  t,
                                  fr_->fcdata.get(),
                                  nullptr,
                                  state.x.rvec_array(),
                                  top_global,
                                  state.box);


// ============================================================================
// COMPILATION ET TEST
// ============================================================================

/*
# 1. Recompiler GROMACS
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make -j8
make install

# 2. Tester avec une simulation
gmx mdrun -s topol.tpr -deffnm test

# 3. Vérifier que Radius-Gyration est dans le .edr
gmx energy -f test.edr

# 4. Sélectionner "Radius-Gyration" et générer le graphique
echo "Radius-Gyration" | gmx energy -f test.edr -o rg.xvg
xmgrace rg.xvg

# 5. Vérifier les performances
# Le calcul du Rg ne doit être fait que tous les nstenergy pas
*/

// ============================================================================
// NOTES FINALES
// ============================================================================

/*
AVANTAGES de cette implémentation :
✓ Calcul du Rg seulement lors de l'écriture (tous les nstenergy pas)
✓ Pas d'overhead à chaque pas de temps
✓ Intégration transparente dans le système ebin_
✓ Compatible avec gmx energy
✓ Valeurs instantanées (pas de moyenne sur les pas intermédiaires)

PERFORMANCE :
- Pour un système de 100k atomes avec nstenergy=100 :
  - AVANT : 100k calculs de Rg (si dans addDataAtEnergyStep)
  - APRÈS : 1k calculs de Rg (si simulation de 100k pas)
  - Gain : ~100x

PERSONNALISATION :
- Remplacez calculateRadiusOfGyration par votre propre fonction
- Adaptez les paramètres passés (x, v, box, etc.) selon vos besoins
- Ajoutez d'autres observables en suivant le même pattern

LIMITATIONS :
- Pas de moyennes automatiques (valeurs instantanées uniquement)
- Nécessite de modifier plusieurs fichiers
- Nécessite une recompilation complète de GROMACS
*/
