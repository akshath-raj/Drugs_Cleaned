"""
Core utility functions for protein structure analysis.
Includes disease mapping, UniProt search, PDB filtering, and structure download.
"""

import requests
from typing import Optional, List, Dict, Tuple
from config import DISEASE_PROTEIN_MAP
import os


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


def search_uniprot_for_reviewed_human(protein_name: str) -> Optional[str]:
    """
    Search UniProt for reviewed (Swiss-Prot) human protein entries.
    Returns the UniProt accession ID.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"(protein_name:{protein_name}) AND (organism_id:9606) AND (reviewed:true)",
        "format": "json",
        "size": 1
    }
    
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        results = data.get('results', [])
        if not results:
            return None
        
        return results[0]['primaryAccession']
        
    except requests.exceptions.RequestException as e:
        print(f"UniProt search error: {e}")
        return None


def get_pdb_ids_from_uniprot(uniprot_id: str) -> List[str]:
    """
    Get all PDB IDs associated with a UniProt entry.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        pdb_ids = []
        cross_refs = data.get('uniProtKBCrossReferences', [])
        
        for ref in cross_refs:
            if ref.get('database') == 'PDB':
                pdb_id = ref.get('id')
                if pdb_id:
                    pdb_ids.append(pdb_id)
        
        return pdb_ids
        
    except requests.exceptions.RequestException as e:
        print(f"Error fetching PDB IDs from UniProt: {e}")
        return []


def get_pdb_resolution(pdb_id: str) -> Optional[float]:
    """
    Get the resolution of a PDB structure.
    Returns None if resolution is not available (e.g., NMR structures) or not X-ray.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        # Check experimental method - MUST be X-ray
        experimental_method = data.get('exptl', [{}])[0].get('method', '')
        if 'X-RAY' not in experimental_method.upper():
            return None
        
        # Try to get resolution from different possible fields
        resolution = None
        
        # Check for X-ray crystallography resolution
        if 'rcsb_entry_info' in data:
            resolution = data['rcsb_entry_info'].get('resolution_combined')
        
        # Also check refine section
        if resolution is None and 'refine' in data:
            refine_list = data['refine']
            if refine_list and len(refine_list) > 0:
                resolution = refine_list[0].get('ls_d_res_high')
        
        # Handle case where resolution might be a list
        if isinstance(resolution, list):
            resolution = resolution[0] if resolution else None
        
        return float(resolution) if resolution is not None else None
        
    except (requests.exceptions.RequestException, ValueError, KeyError, TypeError) as e:
        print(f"Error fetching resolution for {pdb_id}: {e}")
        return None


def check_mutations_in_pdb(pdb_id: str) -> bool:
    """
    Check if a PDB structure contains mutations in ANY chain.
    Returns True if mutations are found, False otherwise.
    """
    # Get the number of polymer entities first
    try:
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry_response = requests.get(entry_url, timeout=30)
        entry_response.raise_for_status()
        entry_data = entry_response.json()
        
        # Get number of polymer entities
        num_entities = entry_data.get('rcsb_entry_info', {}).get('polymer_entity_count', 1)
        
    except requests.exceptions.RequestException as e:
        print(f" [Error getting entity count: {e}]", end="")
        num_entities = 1  # Default to checking at least entity 1
    
    # Check each polymer entity (chain)
    for entity_id in range(1, num_entities + 1):
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        
        try:
            response = requests.get(url, timeout=30)
            
            # If entity doesn't exist, skip to next
            if response.status_code == 404:
                continue
                
            response.raise_for_status()
            data = response.json()
            
            # Check the rcsb_mutation_count field in entity_poly
            if 'entity_poly' in data:
                mutation_count = data['entity_poly'].get('rcsb_mutation_count', 0)
                if mutation_count > 0:
                    return True
            
            # Also check rcsb_polymer_entity for mutation_count (alternative location)
            if 'rcsb_polymer_entity' in data:
                mutation_count = data['rcsb_polymer_entity'].get('mutation_count', 0)
                if mutation_count > 0:
                    return True
            
        except requests.exceptions.RequestException as e:
            # If we can't check this entity, continue to next one
            continue
    
    # No mutations found in any chain
    return False


def download_pdb_file(pdb_id: str, output_dir: str = "proteins") -> Optional[str]:
    """
    Download a PDB file and save it to the specified directory.
    Returns the file path if successful, None otherwise.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        with open(output_path, 'w') as f:
            f.write(response.text)
        
        print(f"Downloaded {pdb_id}.pdb to {output_path}")
        return output_path
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading PDB file {pdb_id}: {e}")
        return None


def find_best_pdb_structure(protein_name: str, max_check: int = 100) -> Optional[Tuple[str, str]]:
    """
    Find the best PDB structure for a protein:
    1. Search UniProt for reviewed human protein
    2. Get associated PDB IDs
    3. Check structures iteratively for X-ray method, resolution and mutations
    4. Download the best structure without mutations
    5. Auto-stop if resolution < 1.5Å and no mutations
    
    Args:
        protein_name: Name of the protein to search for
        max_check: Maximum number of structures to check (default: 100)
    
    Returns: Tuple of (pdb_id, file_path) or None
    """
    print(f"Searching UniProt for: {protein_name}")
    uniprot_id = search_uniprot_for_reviewed_human(protein_name)
    
    if not uniprot_id:
        print(f"No reviewed human protein found for: {protein_name}")
        return None
    
    print(f"Found UniProt ID: {uniprot_id}")
    
    pdb_ids = get_pdb_ids_from_uniprot(uniprot_id)
    
    if not pdb_ids:
        print(f"No PDB structures found for UniProt ID: {uniprot_id}")
        return None
    
    print(f"Found {len(pdb_ids)} PDB structures")
    print(f"Checking up to {min(max_check, len(pdb_ids))} X-ray structures for best resolution without mutations...\n")
    
    # Process structures iteratively to find the best one
    best_candidate = None
    checked_count = 0
    
    for i, pdb_id in enumerate(pdb_ids):
        if checked_count >= max_check:
            print(f"\nReached maximum check limit ({max_check} structures)")
            break
        
        # Get resolution (also checks if X-ray method)
        resolution = get_pdb_resolution(pdb_id)
        if resolution is None:
            continue
        
        checked_count += 1
        print(f"[{checked_count}/{min(max_check, len(pdb_ids))}] {pdb_id}: {resolution}Å (X-ray)", end="")
        
        # ALWAYS check for mutations - this is critical
        has_mutations = check_mutations_in_pdb(pdb_id)
        
        if has_mutations:
            print(" - has mutations, skipping")
            continue
        
        print(" - no mutations ✓")
        
        # Update best candidate if this is better
        if best_candidate is None or resolution < best_candidate[1]:
            best_candidate = (pdb_id, resolution)
            print(f"  → New best candidate!")
            
            # Auto-stop if we found excellent resolution without mutations
            if resolution < 1.5:
                print(f"\n✓ Found excellent resolution (< 1.5Å) without mutations! Stopping search.")
                break
    
    if best_candidate is None:
        print("\nNo suitable structure found (all contain mutations, lack resolution data, or not X-ray)")
        return None
    
    pdb_id, resolution = best_candidate
    print(f"\n✓ Best structure: {pdb_id} (resolution: {resolution}Å, X-ray, no mutations)")
    print(f"Downloading...")
    
    file_path = download_pdb_file(pdb_id)
    
    if file_path:
        return (pdb_id, file_path)
    
    print("Failed to download the selected structure")
    return None


# Legacy function for backwards compatibility
def search_pdb_for_first_hit(protein_name: str) -> Optional[str]:
    """Search RCSB PDB and return the first result found (legacy function)."""
    result = find_best_pdb_structure(protein_name)
    if result:
        return result[0]  # Return just the PDB ID
    return None