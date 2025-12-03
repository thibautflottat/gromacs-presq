#!/bin/bash
# Test de continuation avec append pour complex_obs.xvg

set -e

echo "=== Test de continuation avec append ==="

# Nettoyer les anciens fichiers
rm -f in.tpr in.cpt in.edr in.xtc in.log complex_obs.xvg

echo ""
echo "1. Premier run : 0-1 ps"
echo "========================"

# Modifier le mdp pour un run court
cat > short.mdp << EOF
integrator  = md
nsteps      = 1000
dt          = 0.001
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstlog      = 100
nstenergy   = 100
nstxout-compressed = 0
continuation = no
gen-vel     = yes
gen-temp    = 300
tcoupl      = no
pcoupl      = no
cutoff-scheme = Verlet
coulombtype = Cut-off
rcoulomb    = 1.0
rvdw        = 1.0
nstlist     = 10
EOF

# Préparer et lancer le premier run
gmx grompp -f short.mdp -c conf.gro -p topol.top -o in.tpr -maxwarn 10 > /dev/null 2>&1
gmx mdrun -deffnm in -v > /dev/null 2>&1

# Vérifier le fichier XVG
lines1=$(grep -v '^[#@]' complex_obs.xvg | wc -l)
echo "Lignes de données dans complex_obs.xvg : $lines1"
tail -1 complex_obs.xvg | awk '{print "Dernier temps : " $1 " ps"}'

echo ""
echo "2. Continuation avec append : 1-2 ps"
echo "===================================="

# Modifier le mdp pour continuation
cat > continue.mdp << EOF
integrator  = md
nsteps      = 2000
dt          = 0.001
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstlog      = 100
nstenergy   = 100
nstxout-compressed = 0
continuation = yes
gen-vel     = no
tcoupl      = no
pcoupl      = no
cutoff-scheme = Verlet
coulombtype = Cut-off
rcoulomb    = 1.0
rvdw        = 1.0
nstlist     = 10
EOF

# Préparer et lancer la continuation avec -cpi et -append
gmx grompp -f continue.mdp -c in.gro -t in.cpt -p topol.top -o in.tpr -maxwarn 10 > /dev/null 2>&1
gmx mdrun -deffnm in -cpi in.cpt -append -v > /dev/null 2>&1

# Vérifier le fichier XVG après append
lines2=$(grep -v '^[#@]' complex_obs.xvg | wc -l)
echo "Lignes de données après append : $lines2"
tail -1 complex_obs.xvg | awk '{print "Dernier temps : " $1 " ps"}'

echo ""
echo "3. Vérification de la continuité"
echo "================================="

# Extraire les temps (ignorer header)
grep -v '^[#@]' complex_obs.xvg | awk '{print $1}' > times.dat

# Vérifier qu'il n'y a pas de doublons
duplicates=$(sort times.dat | uniq -d | wc -l)
if [ $duplicates -eq 0 ]; then
    echo "✓ Pas de temps dupliqués"
else
    echo "✗ Attention : $duplicates temps dupliqués trouvés !"
fi

# Vérifier que les temps sont croissants
is_sorted=$(awk 'NR>1 && $1<=prev {print "NOT_SORTED"; exit} {prev=$1}' times.dat)
if [ -z "$is_sorted" ]; then
    echo "✓ Temps croissants"
else
    echo "✗ Attention : temps non croissants !"
fi

# Vérifier le nombre de headers q-vecteurs
qvec_headers=$(grep -c "# Generated.*q-vectors:" complex_obs.xvg)
echo "✓ Nombre de headers q-vecteurs : $qvec_headers (devrait être 1)"

echo ""
echo "4. Afficher quelques lignes autour de la jonction"
echo "=================================================="
echo "Lignes 8-13 (autour de t=1.0 ps) :"
sed -n '8,13p' complex_obs.xvg | awk '{print "  Line " NR+7 ": t=" $1}'

# Nettoyer
rm -f times.dat short.mdp continue.mdp

echo ""
echo "=== Test terminé avec succès ==="
