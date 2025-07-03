"""IO-related functions for MultiModulon class."""

from __future__ import annotations

import json
import pickle
import gzip
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from .core import MultiModulon
    from .species_data import SpeciesData

logger = logging.getLogger(__name__)

# COG category color mapping
COG_COLORS = {
    # INFORMATION STORAGE AND PROCESSING
    'Translation, ribosomal structure and biogenesis': 'black',
    'RNA processing and modification': 'lightblue',
    'Transcription': 'sandybrown',
    'Replication, recombination and repair': 'fuchsia',
    'Chromatin structure and dynamics': 'mediumpurple',
    
    # CELLULAR PROCESSES AND SIGNALING
    'Cell cycle control, cell division, chromosome partitioning': 'y',
    'Nuclear structure': 'goldenrod',
    'Defense mechanisms': 'lightgray',
    'Signal transduction mechanisms': 'lime',
    'Cell wall/membrane/envelope biogenesis': 'mediumvioletred',
    'Cell motility': 'orchid',
    'Cytoskeleton': 'darkviolet',
    'Extracellular structures': 'darkgoldenrod',
    'Intracellular trafficking, secretion, and vesicular transport': 'saddlebrown',
    'Posttranslational modification, protein turnover, chaperones': 'skyblue',
    'Post-translational modification, protein turnover, and chaperones': 'skyblue',  # Alternative spelling
    
    # METABOLISM
    'Energy production and conversion': 'lightgreen',
    'Carbohydrate transport and metabolism': 'pink',
    'Amino acid transport and metabolism': 'red',
    'Nucleotide transport and metabolism': 'c',
    'Coenzyme transport and metabolism': 'green',
    'Lipid transport and metabolism': 'turquoise',
    'Inorganic ion transport and metabolism': 'blue',
    'Secondary metabolites biosynthesis, transport and catabolism': 'dodgerblue',
    'Secondary metabolites biosynthesis, transport, and catabolism': 'dodgerblue',  # Alternative spelling
    
    # POORLY CHARACTERIZED
    'Function unknown': 'slategray',
    
    # No annotation
    'No COG annotation': 'lightskyblue'
}

# COG category letter codes
COG_LETTER_CODES = {
    'J': 'Translation, ribosomal structure and biogenesis',
    'A': 'RNA processing and modification',
    'K': 'Transcription',
    'L': 'Replication, recombination and repair',
    'B': 'Chromatin structure and dynamics',
    'D': 'Cell cycle control, cell division, chromosome partitioning',
    'Y': 'Nuclear structure',
    'V': 'Defense mechanisms',
    'T': 'Signal transduction mechanisms',
    'M': 'Cell wall/membrane/envelope biogenesis',
    'N': 'Cell motility',
    'Z': 'Cytoskeleton',
    'W': 'Extracellular structures',
    'U': 'Intracellular trafficking, secretion, and vesicular transport',
    'O': 'Posttranslational modification, protein turnover, chaperones',
    'C': 'Energy production and conversion',
    'G': 'Carbohydrate transport and metabolism',
    'E': 'Amino acid transport and metabolism',
    'F': 'Nucleotide transport and metabolism',
    'H': 'Coenzyme transport and metabolism',
    'I': 'Lipid transport and metabolism',
    'P': 'Inorganic ion transport and metabolism',
    'Q': 'Secondary metabolites biosynthesis, transport and catabolism',
    'S': 'Function unknown'
}


def save_bbh(multi_modulon: 'MultiModulon', output_path: str) -> None:
    """
    Save BBH results to file.
    
    Parameters
    ----------
    multi_modulon : MultiModulon
        The MultiModulon instance containing BBH results
    output_path : str
        Path to save BBH results
    """
    if multi_modulon._bbh is None:
        logger.warning("No BBH results to save")
        return
    
    output_path = Path(output_path)
    
    # Convert to serializable format
    bbh_data = {}
    for (sp1, sp2), df in multi_modulon._bbh.items():
        key = f"{sp1}___{sp2}"
        bbh_data[key] = df.to_dict('records')
    
    # Save as JSON
    with open(output_path.with_suffix('.json'), 'w') as f:
        json.dump(bbh_data, f, indent=2)
    
    # Also save as pickle for easier loading
    with open(output_path.with_suffix('.pkl'), 'wb') as f:
        pickle.dump(multi_modulon._bbh, f)
    
    logger.info(f"BBH results saved to {output_path}")


def load_bbh(multi_modulon: 'MultiModulon', input_path: str) -> None:
    """
    Load BBH results from file.
    
    Parameters
    ----------
    multi_modulon : MultiModulon
        The MultiModulon instance to load BBH results into
    input_path : str
        Path to load BBH results from
    """
    input_path = Path(input_path)
    
    if input_path.with_suffix('.pkl').exists():
        # Load from pickle
        with open(input_path.with_suffix('.pkl'), 'rb') as f:
            multi_modulon._bbh = pickle.load(f)
    elif input_path.with_suffix('.json').exists():
        # Load from JSON
        with open(input_path.with_suffix('.json'), 'r') as f:
            bbh_data = json.load(f)
        
        multi_modulon._bbh = {}
        for key, records in bbh_data.items():
            sp1, sp2 = key.split('___')
            multi_modulon._bbh[(sp1, sp2)] = pd.DataFrame(records)
    else:
        raise FileNotFoundError(f"BBH file not found: {input_path}")
    
    logger.info(f"BBH results loaded from {input_path}")


def get_orthologs(multi_modulon: 'MultiModulon', species1: str, species2: str) -> pd.DataFrame:
    """
    Get ortholog pairs between two species.
    
    Parameters
    ----------
    multi_modulon : MultiModulon
        The MultiModulon instance containing BBH results
    species1, species2 : str
        Species names
        
    Returns
    -------
    pd.DataFrame
        Ortholog pairs with columns: query_id, subject_id, evalue, bitscore
    """
    if multi_modulon._bbh is None:
        raise ValueError("BBH analysis not run. Initialize with run_bbh=True or load BBH results.")
    
    key = (species1, species2)
    if key not in multi_modulon._bbh:
        raise ValueError(f"No BBH results found for {species1} vs {species2}")
    
    return multi_modulon._bbh[key]


def save_to_json_multimodulon(multi_modulon: 'MultiModulon', save_path: str) -> None:
    """
    Save the entire MultiModulon object to a JSON file.
    
    This function serializes all data including:
    - All species data (log_tpm, log_tpm_norm, X, M, A, gene_table, sample_sheet, M_thresholds, presence_matrix)
    - Combined gene database
    - Input folder path
    
    Parameters
    ----------
    multi_modulon : MultiModulon
        The MultiModulon instance to save
    save_path : str
        Path to save the JSON file. If the path ends with '.json.gz', 
        the file will be gzip compressed.
        
    Examples
    --------
    >>> save_to_json_multimodulon(multiModulon, "multimodulon_data.json")
    >>> save_to_json_multimodulon(multiModulon, "multimodulon_data.json.gz")  # Compressed
    """
    save_path = Path(save_path)
    
    # Create the data structure to save
    data = {
        'input_folder_path': str(multi_modulon.input_folder_path),
        'species_data': {},
        'combined_gene_db': None
    }
    
    # Save combined gene database if available
    if hasattr(multi_modulon, '_combined_gene_db') and multi_modulon._combined_gene_db is not None:
        data['combined_gene_db'] = {
            'data': multi_modulon._combined_gene_db.to_dict('records'),
            'columns': list(multi_modulon._combined_gene_db.columns)
        }
    
    
    # Save each species data
    for species_name, species_data_obj in multi_modulon._species_data.items():
        species_dict = {
            'species_name': species_data_obj.species_name,
            'data_path': str(species_data_obj.data_path),
            'log_tpm': None,
            'log_tpm_norm': None,
            'X': None,
            'M': None,
            'A': None,
            'sample_sheet': None,
            'gene_table': None,
            'M_thresholds': None,
            'presence_matrix': None
        }
        
        # Save DataFrames if they exist
        if species_data_obj._log_tpm is not None:
            species_dict['log_tpm'] = {
                'data': species_data_obj._log_tpm.to_dict('records'),
                'index': list(species_data_obj._log_tpm.index),
                'columns': list(species_data_obj._log_tpm.columns)
            }
        
        if species_data_obj._log_tpm_norm is not None:
            species_dict['log_tpm_norm'] = {
                'data': species_data_obj._log_tpm_norm.to_dict('records'),
                'index': list(species_data_obj._log_tpm_norm.index),
                'columns': list(species_data_obj._log_tpm_norm.columns)
            }
        
        if species_data_obj._X is not None:
            species_dict['X'] = {
                'data': species_data_obj._X.to_dict('records'),
                'index': list(species_data_obj._X.index),
                'columns': list(species_data_obj._X.columns)
            }
        
        if species_data_obj._M is not None:
            species_dict['M'] = {
                'data': species_data_obj._M.to_dict('records'),
                'index': list(species_data_obj._M.index),
                'columns': list(species_data_obj._M.columns)
            }
        
        if species_data_obj._A is not None:
            species_dict['A'] = {
                'data': species_data_obj._A.to_dict('records'),
                'index': list(species_data_obj._A.index),
                'columns': list(species_data_obj._A.columns)
            }
        
        if species_data_obj._sample_sheet is not None:
            species_dict['sample_sheet'] = {
                'data': species_data_obj._sample_sheet.to_dict('records'),
                'index': list(species_data_obj._sample_sheet.index),
                'columns': list(species_data_obj._sample_sheet.columns)
            }
        
        if species_data_obj._gene_table is not None:
            species_dict['gene_table'] = {
                'data': species_data_obj._gene_table.to_dict('records'),
                'index': list(species_data_obj._gene_table.index),
                'columns': list(species_data_obj._gene_table.columns)
            }
        
        if species_data_obj._M_thresholds is not None:
            species_dict['M_thresholds'] = {
                'data': species_data_obj._M_thresholds.to_dict('records'),
                'index': list(species_data_obj._M_thresholds.index),
                'columns': list(species_data_obj._M_thresholds.columns)
            }
        
        if species_data_obj._presence_matrix is not None:
            species_dict['presence_matrix'] = {
                'data': species_data_obj._presence_matrix.to_dict('records'),
                'index': list(species_data_obj._presence_matrix.index),
                'columns': list(species_data_obj._presence_matrix.columns)
            }
        
        data['species_data'][species_name] = species_dict
    
    # Save to JSON file with optional compression
    if str(save_path).endswith('.json.gz'):
        # Save as compressed gzip file
        with gzip.open(save_path, 'wt', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        logger.info(f"MultiModulon object saved to {save_path} (compressed)")
    else:
        # Save as regular JSON file
        with open(save_path, 'w') as f:
            json.dump(data, f, indent=2)
        logger.info(f"MultiModulon object saved to {save_path}")


def load_json_multimodulon(load_path: str) -> 'MultiModulon':
    """
    Load a MultiModulon object from a JSON file.
    
    This function creates a new MultiModulon object and populates it
    with all the data saved in the JSON file.
    
    Parameters
    ----------
    load_path : str
        Path to the JSON file to load. If the path ends with '.json.gz',
        the file will be treated as gzip compressed.
        
    Returns
    -------
    MultiModulon
        A new MultiModulon object with all loaded data
        
    Examples
    --------
    >>> multiModulon = load_json_multimodulon("multimodulon_data.json")
    >>> multiModulon = load_json_multimodulon("multimodulon_data.json.gz")  # Compressed
    """
    from .core import MultiModulon
    from .species_data import SpeciesData
    
    load_path = Path(load_path)
    
    # Load JSON data with optional decompression
    if str(load_path).endswith('.json.gz'):
        # Load from compressed gzip file
        with gzip.open(load_path, 'rt', encoding='utf-8') as f:
            data = json.load(f)
    else:
        # Load from regular JSON file
        with open(load_path, 'r') as f:
            data = json.load(f)
    
    # Create a new MultiModulon object
    # We need to handle the case where the input folder might not exist
    input_folder = data['input_folder_path']
    
    # Create a minimal MultiModulon object
    multi_modulon = MultiModulon.__new__(MultiModulon)
    multi_modulon.input_folder_path = Path(input_folder)
    multi_modulon._species_data = {}
    multi_modulon._bbh = None
    
    # Load combined gene database
    if data.get('combined_gene_db') is not None:
        df_data = data['combined_gene_db']
        multi_modulon._combined_gene_db = pd.DataFrame(df_data['data'])
        if df_data['data']:  # Only set columns if there's data
            multi_modulon._combined_gene_db = multi_modulon._combined_gene_db[df_data['columns']]
    else:
        multi_modulon._combined_gene_db = None
    
    
    # Load species data
    for species_name, species_dict in data['species_data'].items():
        # Create SpeciesData object
        species_data = SpeciesData.__new__(SpeciesData)
        species_data.species_name = species_dict['species_name']
        species_data.data_path = Path(species_dict['data_path'])
        
        # Load DataFrames
        def load_dataframe(df_data):
            if df_data is None:
                return None
            df = pd.DataFrame(df_data['data'])
            if df_data['data']:  # Only set index/columns if there's data
                df.index = df_data['index']
                df.columns = df_data['columns']
            return df
        
        species_data._log_tpm = load_dataframe(species_dict.get('log_tpm'))
        species_data._log_tpm_norm = load_dataframe(species_dict.get('log_tpm_norm'))
        species_data._X = load_dataframe(species_dict.get('X'))
        species_data._M = load_dataframe(species_dict.get('M'))
        species_data._A = load_dataframe(species_dict.get('A'))
        species_data._sample_sheet = load_dataframe(species_dict.get('sample_sheet'))
        species_data._gene_table = load_dataframe(species_dict.get('gene_table'))
        species_data._M_thresholds = load_dataframe(species_dict.get('M_thresholds'))
        species_data._presence_matrix = load_dataframe(species_dict.get('presence_matrix'))
        
        # Add to MultiModulon
        multi_modulon._species_data[species_name] = species_data
    
    logger.info(f"MultiModulon object loaded from {load_path}")
    logger.info(f"Loaded {len(multi_modulon._species_data)} species")
    
    return multi_modulon