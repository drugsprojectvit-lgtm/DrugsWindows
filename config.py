"""
Configuration file for protein structure analysis application.
Contains disease-protein mappings, API tokens, and global state.
"""

# Disease to Protein mapping
DISEASE_PROTEIN_MAP = {
    "inflammation": {
        "anti-inflammatory": "Cyclooxygenase-2",
        "cox-2": "Cyclooxygenase-2",
        "pro-inflammatory": "Tumor necrosis factor-alpha",
        "rheumatoid arthritis": "Tumor necrosis factor-alpha",
        "tnf-alpha": "Tumor necrosis factor-alpha"
    },
    "oncology": {
        "neuro tumor": "Vascular endothelial growth factor receptor 2",
        "glioblastoma": "Vascular endothelial growth factor receptor 2",
        "vegfr-2": "Vascular endothelial growth factor receptor 2"
    },
    "metabolic": {
        "pre-diabetes": "Dipeptidyl peptidase 4",
        "diabetes": "Dipeptidyl peptidase 4",
        "dpp-4": "Dipeptidyl peptidase 4",
        "sglt2": "Sodium/glucose cotransporter 2",
        "obesity": "Glucagon-like peptide 1 receptor",
        "glp-1r": "Glucagon-like peptide 1 receptor"
    },
    "neurodegenerative": {
        "parkinson's disease": "Leucine-rich repeat kinase 2",
        "parkinsons disease": "Leucine-rich repeat kinase 2",
        "lrrk2": "Leucine-rich repeat kinase 2",
        "alzheimer's disease": "Beta-secretase 1",
        "alzheimers disease": "Beta-secretase 1",
        "bace1": "Beta-secretase 1"
    }
}

# SWISS-MODEL Configuration
API_TOKEN = "9e8b3ac03b851bb3834cdb311045c78021087d1d"
BASE_URL = "https://swissmodel.expasy.org"
HEADERS = {"Authorization": f"Token {API_TOKEN}"}

# Global variable to store current PDB info
current_pdb_info = {
    "pdb_id": None, 
    "pdb_path": None, 
    "prepared_pdbqt": None, 
    "docking_results": None,
    "prankweb_csv": None
}

# Directory paths
PROTEINS_DIR = "proteins"
RAMPLOT_OUTPUT_DIR = "my_analysis_folder"
PRANKWEB_OUTPUT_DIR = "prankweb_results"
PREPARED_PROTEIN_DIR = "prepared_protein_meeko"
DOCKING_RESULTS_DIR = "docking_results"
LIGAND_DIR = "pdbqt"
