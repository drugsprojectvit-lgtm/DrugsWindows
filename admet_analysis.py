import os
import glob
import pandas as pd
import numpy as np
from rdkit import Chem
from adme_py import ADME
from admet_ai import ADMETModel

# ==========================================
# 1. CONFIGURATION (UPDATED)
# ==========================================

# Reduced column list based on your UI requirement image
DISPLAY_COLUMNS = [
    "Filename",                # Identifier
    "Docking Score",           # 1. Docking score
    "Final Decision",          # 2. Final decision
    "SA Score",                # 3. SA score
    "QED",                     # 4. QED
    "hERG",                    # 5a. hERG
    "Ames",                    # 5b. Ames
    "CYP3A4 Inhibition",       # 5c. CYP Flags (Part 1)
    "CYP2D6 Inhibition",       # 5d. CYP Flags (Part 2)
    "Developability Score"     # 6. Composite score
]

# ==========================================
# 2. HELPER FUNCTIONS
# ==========================================

def pdb_to_smiles(pdb_path):
    """Convert PDB file to SMILES string using RDKit."""
    try:
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    except Exception:
        return None

def extract_docking_score(pdb_file):
    """Extract docking score from PDB file REMARK lines."""
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    if len(parts) >= 4:
                        return float(parts[3])
    except:
        pass
    return -7.0  # Default mock value if not found

# ==========================================
# 3. FILTERING LOGIC
# ==========================================

def apply_primary_filters(row):
    """
    Stage 1: Primary Filters (Hard Stops - Auto Reject)
    Returns: ('PASS' | 'REVIEW' | 'REJECT', reason_string)
    """
    # --- HARD REJECTS ---
    if row.get('Lipinski') == 'Fail': 
        return 'REJECT', 'Lipinski Fail'
    
    if row.get('PAINS') == 'Yes': 
        return 'REJECT', 'PAINS Alert'
    
    if row.get('Brenk') == 'Yes': 
        return 'REJECT', 'Brenk Alert'
    
    if row.get('Ames') in ['Positive', 1]: 
        return 'REJECT', 'Ames Positive'
    
    # hERG Logic
    herg_val = row.get('hERG', '')
    if herg_val == 'High' or (isinstance(herg_val, (int, float)) and herg_val >= 0.7):
        return 'REJECT', 'hERG High Risk'
    elif herg_val == 'Medium' or (isinstance(herg_val, (int, float)) and 0.3 <= herg_val < 0.7):
        return 'REVIEW', 'hERG Medium Risk'

    # Carcinogenicity
    if row.get('Carcinogenicity') in ['Yes', 1]: 
        return 'REJECT', 'Carcinogenic'
    
    # DILI Logic
    dili_val = row.get('DILI', '')
    if dili_val == 'High' or (isinstance(dili_val, (int, float)) and dili_val >= 0.7):
        return 'REJECT', 'DILI High Risk'
    elif dili_val == 'Medium' or (isinstance(dili_val, (int, float)) and 0.3 <= dili_val < 0.7):
        return 'REVIEW', 'DILI Medium Risk'
    
    return 'PASS', None

def apply_developability_filters(row):
    """
    Stage 2: Developability Filters (Soft Constraints)
    Returns: ('ACCEPT' | 'REVIEW' | 'REJECT', reason_string)
    """
    reject_reasons = []
    review_reasons = []

    # Helper to safely get floats
    def get_float(key, default):
        try:
            return float(row.get(key, default))
        except:
            return default

    sa = get_float('SA Score', 5.0)
    qed = get_float('QED', 0.6)
    mw = get_float('MW', 400)
    tpsa = get_float('TPSA', 100)
    wlogp = get_float('WLogP', 3)

    # SA Score
    if sa > 6.0: reject_reasons.append('SA > 6.0')
    elif 5.0 < sa <= 6.0: review_reasons.append('SA 5.0-6.0')
    
    # QED
    if qed < 0.4: reject_reasons.append('QED < 0.4')
    elif 0.4 <= qed < 0.6: review_reasons.append('QED 0.4-0.6')
    
    # MW
    if mw > 550: reject_reasons.append('MW > 550')
    elif 500 <= mw <= 550: review_reasons.append('MW 500-550')
    
    # TPSA
    if tpsa > 160: reject_reasons.append('TPSA > 160')
    elif 140 <= tpsa <= 160: review_reasons.append('TPSA 140-160')
    
    # WLogP
    if wlogp > 6: reject_reasons.append('WLogP > 6')
    elif 5 <= wlogp <= 6: review_reasons.append('WLogP 5-6')
    
    if reject_reasons:
        return 'REJECT', '; '.join(reject_reasons)
    elif review_reasons:
        return 'REVIEW', '; '.join(review_reasons)
    else:
        return 'ACCEPT', None

def calculate_adme_penalties(row):
    """
    Stage 3: ADME Risk Flags (Penalty Based)
    Returns: Penalty points to subtract.
    """
    penalty = 0
    
    # 1. GI Absorption
    gi = row.get('GI Absorption', 'High')
    if gi in ['Medium', 'Moderate']:
        penalty += 5
    
    # 2. Caco-2
    caco2 = row.get('Caco-2 (Wang)', '')
    if caco2 == 'Moderate':
        penalty += 5
    
    # 3. BBB
    bbb = row.get('BBB (Martins)', '')
    if bbb == 'Borderline':
        penalty += 5

    # 4. PPB
    ppb = row.get('PPB (AZ)', 90)
    try: ppb = float(ppb)
    except: ppb = 90
    if ppb >= 95:
        penalty += 5
    
    # 5. CYP Inhibitions
    if row.get('CYP3A4 Inhibition') in ['Yes', 'Weak', 1]:
        penalty += 5
    if row.get('CYP2D6 Inhibition') in ['Yes', 'Weak', 1]:
        penalty += 5
    
    # 6. hERG Medium Risk
    herg = row.get('hERG', '')
    if herg == 'Medium' or (isinstance(herg, (int, float)) and 0.3 <= herg < 0.7):
        penalty += 10
    
    return penalty

def calculate_developability_score(row, docking_score):
    """
    Stage 4: Composite Score Calculation
    """
    base_score = 100
    penalty = 0
    
    # SA Score Penalty
    sa = row.get('SA Score', 5.0)
    try: sa = float(sa)
    except: sa = 5.0
    
    if 5.0 <= sa <= 6.0:
        penalty += 10
    elif sa > 6.0:
        return 0 
    
    # QED Penalty
    qed = row.get('QED', 0.6)
    try: qed = float(qed)
    except: qed = 0.6
    
    if qed < 0.6:
        penalty += 10
    
    # Add ADME penalties
    penalty += calculate_adme_penalties(row)
    
    final_score = base_score - penalty
    return max(0, final_score)

def make_final_decision(developability_score):
    """
    Stage 5: Score-based Decision
    """
    if developability_score >= 75:
        return 'ACCEPT'
    elif 60 <= developability_score < 75:
        return 'REVIEW'
    else:
        return 'REJECT'

# ==========================================
# 4. MAIN EXECUTION PIPELINE
# ==========================================

def run_admet_prediction():
    print("--- Starting ADMET Prediction Pipeline ---")
    
    # 1. FIND FILES
    search_path = os.path.join("docking_results", "pdb", "*ligand_poses.pdb")
    pdb_files = glob.glob(search_path)
    
    if not pdb_files:
        print("No ligand PDB files found.")
        return None

    # 2. PREPARE DATA
    compounds = []
    smiles_list = []
    
    print(f"Found {len(pdb_files)} files. Extracting SMILES...")
    for pdb_file in pdb_files:
        smiles = pdb_to_smiles(pdb_file)
        if smiles:
            compounds.append({
                "Filename": os.path.basename(pdb_file),
                "SMILES": smiles,
                "Docking Score": extract_docking_score(pdb_file)
            })
            smiles_list.append(smiles)
            
    if not compounds:
        print("Failed to extract SMILES.")
        return None

    df_base = pd.DataFrame(compounds)

    # 3. RUN ADME-PY (Physicochemical)
    print("Running ADME-Py...")
    adme_data = []
    for cmp in compounds:
        row = {}
        try:
            res = ADME(cmp["SMILES"]).calculate()
            p = res.get("physiochemical", {})
            m = res.get("medicinal", {})
            pk = res.get("pharmacokinetics", {})
            
            row["MW"] = p.get("molecular_weight")
            row["TPSA"] = p.get("tpsa")
            row["HBD"] = p.get("num_h_donors")
            row["HBA"] = p.get("num_h_acceptors")
            row["Rotatable Bonds"] = p.get("num_rotatable_bonds")
            row["Fsp3"] = p.get("sp3_carbon_ratio")
            row["WLogP"] = res.get("lipophilicity", {}).get("wlogp")
            row["GI Absorption"] = pk.get("gastrointestinal_absorption")
            row["Lipinski"] = "Pass" if res.get("druglikeness", {}).get("lipinski") else "Fail"
            row["PAINS"] = "Yes" if m.get("pains") else "No"
            row["Brenk"] = "Yes" if m.get("brenk") else "No"
            row["SA Score"] = m.get("synthetic_accessibility")
        except:
            pass
        adme_data.append(row)
    df_adme = pd.DataFrame(adme_data)

    # 4. RUN ADMET-AI (Toxicity/Metabolism)
    print("Running ADMET-AI...")
    try:
        model = ADMETModel()
        preds_df = model.predict(smiles_list)
        
        rename_map = {
            "Caco2_Wang": "Caco-2 (Wang)",
            "BBB_Martins": "BBB (Martins)",
            "PPBR_AZ": "PPB (AZ)",
            "CYP3A4_Veith": "CYP3A4 Inhibition",
            "CYP2D6_Veith": "CYP2D6 Inhibition",
            "hERG": "hERG",
            "AMES": "Ames",
            "DILI": "DILI",
            "Carcinogens_Lagunin": "Carcinogenicity",
            "LD50_Zhu": "LD50 (Zhu)",
            "QED": "QED"
        }
        available_cols = [c for c in rename_map.keys() if c in preds_df.columns]
        admet_mapped = preds_df[available_cols].rename(columns=rename_map)
    except Exception as e:
        print(f"ADMET-AI Failed: {e}")
        # Create empty DataFrame with same length as input if prediction fails
        admet_mapped = pd.DataFrame(index=range(len(compounds)))

    # 5. MERGE (FIXED: Reset index to avoid conflicts)
    # Ensure all DataFrames have standard RangeIndex before concat
    df_base.reset_index(drop=True, inplace=True)
    df_adme.reset_index(drop=True, inplace=True)
    admet_mapped.reset_index(drop=True, inplace=True)

    # Force strict concatenation by position (ignoring index labels)
    final_df = pd.concat([df_base, df_adme, admet_mapped], axis=1, ignore_index=False)

    # 6. APPLY PIPELINE LOGIC
    print("Applying Decision Filters...")
    decisions = []
    dev_scores = []

    for idx, row in final_df.iterrows():
        # A. Primary Filters
        p_res, p_reason = apply_primary_filters(row)
        if p_res == 'REJECT':
            decisions.append(f"REJECT ({p_reason})")
            dev_scores.append(0)
            continue
        
        # B. Developability Filters
        d_res, d_reason = apply_developability_filters(row)
        if d_res == 'REJECT':
            decisions.append(f"REJECT ({d_reason})")
            dev_scores.append(0)
            continue
        
        # C. Calculate Score
        score = calculate_developability_score(row, row.get('Docking Score'))
        dev_scores.append(score)
        
        # D. Final Decision & OVERRIDE LOGIC
        final_dec = make_final_decision(score)
        
        warnings = []
        if p_res == 'REVIEW': warnings.append(p_reason)
        if d_res == 'REVIEW': warnings.append(d_reason)
        
        if warnings:
            warning_text = "; ".join(warnings)
            if final_dec == 'ACCEPT':
                decisions.append(f"REVIEW ({warning_text})")
            else:
                decisions.append(f"{final_dec} ({warning_text})")
        else:
            decisions.append(final_dec)

    final_df['Developability Score'] = dev_scores
    final_df['Final Decision'] = decisions

    # 7. CLEANUP & SAVE
    # Filter columns to only those present in DISPLAY_COLUMNS
    cols_to_keep = [c for c in DISPLAY_COLUMNS if c in final_df.columns]
    final_df_clean = final_df[cols_to_keep]

    # Round numeric columns safely
    numeric_cols = final_df_clean.select_dtypes(include=['float64', 'float32']).columns
    final_df_clean.loc[:, numeric_cols] = final_df_clean[numeric_cols].round(2)

    os.makedirs("results", exist_ok=True)
    csv_path = os.path.join("results", "final_admet_report.csv")
    final_df_clean.to_csv(csv_path, index=False)
    
    msg = f"Analysis Complete. Processed {len(final_df_clean)} ligands."
    return msg, final_df_clean, csv_path