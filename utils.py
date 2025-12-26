"""
Core utility functions for protein structure analysis.
Includes disease mapping, PDB search, and structure cleaning.
"""

import requests
from typing import Optional, Tuple, Dict
from config import DISEASE_PROTEIN_MAP


def map_disease_to_protein(disease_input: str) -> Optional[str]:
    """Map a disease name or condition to its protein target."""
    disease_input = disease_input.lower().strip()
    
    for category, conditions in DISEASE_PROTEIN_MAP.items():
        if disease_input in conditions:
            return conditions[disease_input]
        if disease_input == category:
            return list(conditions.values())[0]
    
    for category, conditions in DISEASE_PROTEIN_MAP.items():
        for condition_key, protein_name in conditions.items():
            if disease_input in condition_key or condition_key in disease_input:
                return protein_name
    
    return None


def search_pdb_for_first_hit(protein_name: str) -> Optional[str]:
    """Search RCSB PDB and return the first result found."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "struct.title",
                "operator": "contains_phrase",
                "value": protein_name
            }
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}]
        }
    }
    
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    try:
        response = requests.post(url, json=query, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        result_set = data.get('result_set', [])
        if not result_set:
            return None
            
        return result_set[0]['identifier']
        
    except requests.exceptions.RequestException:
        return None


def remove_ligands_from_pdb(pdb_content: str, keep_chain: str = 'A') -> Tuple[str, Dict]:
    """Remove ligands (HETATM) and keep only a single chain."""
    lines = pdb_content.split('\n')
    cleaned_lines = []
    stats = {
        'hetatm_removed': 0,
        'atoms_kept': 0,
        'chains_removed': set(),
        'conect_removed': 0
    }

    atom_serials = set()

    for line in lines:
        if line.startswith('ATOM'):
            try:
                chain_id = line[21].strip()
                serial = int(line[6:11].strip())
                if chain_id == keep_chain:
                    atom_serials.add(serial)
            except (ValueError, IndexError):
                pass

    for line in lines:
        if line.startswith('ATOM'):
            chain_id = line[21].strip()
            if chain_id == keep_chain:
                cleaned_lines.append(line)
                stats['atoms_kept'] += 1
            else:
                stats['chains_removed'].add(chain_id)

        elif line.startswith('HETATM'):
            stats['hetatm_removed'] += 1
            continue

        elif line.startswith('CONECT'):
            try:
                parts = line.split()
                if len(parts) > 1:
                    serials = [int(x) for x in parts[1:] if x.isdigit()]
                    if all(s in atom_serials for s in serials):
                        cleaned_lines.append(line)
                    else:
                        stats['conect_removed'] += 1
            except (ValueError, IndexError):
                stats['conect_removed'] += 1
                continue

        elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS',
                             'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
                             'SEQRES', 'MODRES', 'HELIX', 'SHEET', 'CRYST1',
                             'ORIGX', 'SCALE', 'MTRIX', 'MODEL', 'ENDMDL',
                             'MASTER', 'END', 'TER')):
            cleaned_lines.append(line)

    return '\n'.join(cleaned_lines), stats
