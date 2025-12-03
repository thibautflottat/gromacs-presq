#!/usr/bin/env python3
"""
Nettoyer les doublons de temps dans complex_obs.xvg après un restart avec -append.

Lors d'un restart GROMACS avec -append, la frame du checkpoint est réécrite,
créant un doublon. Ce script garde la première occurrence de chaque temps.
"""

import sys

def clean_xvg(input_file, output_file=None):
    """
    Nettoie un fichier XVG en supprimant les doublons de temps.
    
    Args:
        input_file: Chemin du fichier XVG d'entrée
        output_file: Chemin du fichier XVG de sortie (None = remplacer input)
    """
    # Lire le fichier en séparant header et données
    header_lines = []
    data_lines = []
    times_seen = set()
    
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                header_lines.append(line)
            else:
                # Ligne de données
                cols = line.split()
                if len(cols) > 0:
                    time = float(cols[0])
                    if time not in times_seen:
                        times_seen.add(time)
                        data_lines.append(line)
                    else:
                        print(f"Doublon supprimé : t={time:.6f} ps", file=sys.stderr)
    
    # Écrire le résultat
    if output_file is None:
        output_file = input_file
    
    with open(output_file, 'w') as f:
        f.writelines(header_lines)
        f.writelines(data_lines)
    
    n_removed = len([line for line in open(input_file) 
                     if not (line.startswith('#') or line.startswith('@'))]) - len(data_lines)
    
    print(f"✓ Fichier nettoyé : {len(data_lines)} frames uniques", file=sys.stderr)
    print(f"  ({n_removed} doublon(s) supprimé(s))", file=sys.stderr)
    
    return len(data_lines), n_removed

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python clean_xvg_duplicates.py complex_obs.xvg [output.xvg]")
        print("  Si output.xvg n'est pas spécifié, le fichier d'entrée est remplacé.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    clean_xvg(input_file, output_file)
