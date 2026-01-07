from Bio import PDB

input_pdb = "batch_results\KRAS\\01_structure_search\\4LPK.pdb"
protein_out = "protein.pdb"
ligand_out = "ligands.pdb"

parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure("complex", input_pdb)

class ProteinSelect(PDB.Select):
    def accept_residue(self, residue):
        # Standard protein residues
        return residue.id[0] == " "

class LigandSelect(PDB.Select):
    def accept_residue(self, residue):
        # Docked ligands, ions, waters, cofactors
        return residue.id[0] != " "

io = PDB.PDBIO()
io.set_structure(structure)

# Write protein
io.save(protein_out, ProteinSelect())

# Write ligands
io.save(ligand_out, LigandSelect())

print("Split completed:")
print(" - protein.pdb")
print(" - ligands.pdb")
