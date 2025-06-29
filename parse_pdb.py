from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import sys
from collections import namedtuple

Missing = namedtuple('Missing', ["ifmissing", "missing_ress", "num_per_chain"])
Residue = namedtuple('Residue', ["chain", "name", "num"])
Ligand = namedtuple('Ligand', ["name", "numatom", "chain", "syn"])

def sep_lig_rec(pdb):
    with open(pdb, 'r') as f:
        lines = f.readlines()
        title = None
        res = None
        chain = []
        missing_res = Missing(ifmissing=False, missing_ress=[], num_per_chain=[])
        ligands = []

        for line in lines:
            if line.startswith('TITLE'):
                title = ' '.join(line.split()[1:])
            if line.startswith('REMARK   2 RESOLUTION'):
                res = line.split()[3]
            if line.startswith('COMPND   3 CHAIN'):
                chain = [x.rstrip(';|,') for x in line.split()[3:]]
            if line.startswith('REMARK 465') and 'MISSING RESIDUES' in line:
                missing_res = missing_res._replace(ifmissing=True)
            if line.startswith('REMARK 465') and missing_res.ifmissing:
                if any(keyword in line for keyword in ['MISSING RESIDUES', 'RES C SSSEQI', '---', 'RESIDUES', 'EXPERIMENT', 'IDENTIFIER']):
                    continue
                miss_line = line.split()
                if len(miss_line) >= 5 and miss_line[0] == 'REMARK' and miss_line[1] == '465':
                    missing_residue = Residue(chain=miss_line[3], name=miss_line[2], num=miss_line[4])
                    missing_res.missing_ress.append(missing_residue)
            if line.startswith('HET   '):
                parts = line.split()
                ligands.append(Ligand(name=parts[1], numatom=parts[4], chain=parts[2], syn=[]))
            if line.startswith('HETNAM') or line.startswith('HETSYN'):
                parts = line.split()
                for ligand in ligands:
                    if ligand.name in parts:
                        for item in parts[1:]:
                            ligand.syn.append(item)
                        ligand.syn.sort()
        
        for chain_ in chain:
            missing_in_chain_ = []
            for missing in missing_res.missing_ress:
                if missing.chain == chain_:
                    missing_in_chain_.append((missing.name, missing.num))
            missing_res.num_per_chain.append((chain_, len(missing_in_chain_)))

        print(title)
        print(res)
        print(chain)
        if missing_res.ifmissing:
            print("There are missing residues")
            for chain__ in missing_res.num_per_chain:
                print(f"for chain {chain__[0]}: {chain__[1]}")

        chosen_chain = 'A'
        if len(chain) > 1:
            chosen_chain = input(f"Choose the chain you want {chain}: ")
            print("you have chosen the chain, ", chosen_chain)
            # Create cleaned PDB with headers preserved
            with open(f"{pdb.split('.')[0]}_cleaned.pdb", 'w') as f:
                # Write all header lines first
                for line in lines:
                    if line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 
                                      'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK', 
                                      'DBREF', 'SEQADV', 'SEQRES', 'HET   ', 'HETNAM', 
                                      'HETSYN', 'FORMUL', 'CRYST1', 'ORIGX', 'SCALE')):
                        f.write(line)
                
                # Write ATOM records for chosen chain
                for line in lines:
                    if line.startswith('ATOM'):
                        if line.split()[4] == chosen_chain:
                            f.write(line)
                
                # Write END record
                f.write('END\n')
                f.close()

        for ligand in ligands:
            if len(chain) > 1:
                if ligand.chain == chosen_chain:
                    print(f"{ligand.name}, number of atoms: {ligand.numatom}, synonyms: {ligand.syn}")
            else:
                print(f"{ligand.name}, number of atoms: {ligand.numatom}, synonyms: {ligand.syn}")
        
        chosen_lig = None
        if ligands and len(ligands) > 1: 
            chosen_lig = input("chose the ligand you want: ")

            ligand_lines = []
            for line in lines:
                if line.startswith(('HETATM', 'ATOM')) and chosen_lig in line:
                    if chosen_chain:
                        line_chain = line[21] if len(line) > 21 else ' '
                        if line_chain == chosen_chain:
                            ligand_lines.append(line)
                    else:
                        ligand_lines.append(line)
        
            ligand_pdb = "".join(ligand_lines) + 'END\n'
            with open(f"{chosen_lig}.pdb", 'w') as f:
                f.write(ligand_pdb)
        
        return chosen_chain, chosen_lig


def get_ideal_lig_rcsb(ligand):
    res = requests.get(f"https://files.rcsb.org/ligands/download/{ligand}_ideal.sdf")
    with open(f"{ligand}_ideal.sdf", 'w') as f:
        f.write(res.text)


def fix_lig(lig):
    get_ideal_lig_rcsb(lig)
    mole_i = Chem.SDMolSupplier(f"{lig}_ideal.sdf")[0]
    mole_ = Chem.MolFromPDBFile(f"{lig}.pdb")
    AllChem.AssignBondOrdersFromTemplate(mole_i, mole_)
    writer = Chem.SDWriter(f"{lig}.sdf")
    writer.write(mole_)
    writer.close()

def ligand_to_pdbqt(sdf):
    pass

def ligand_to_mol2():
    pass

def fix_missing_res(pdb):
    fixer = PDBFixer(filename=f"{pdb}.pdb")
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    PDBFile.writeFile(fixer.topology, fixer.positions, file=open(f'{pdb}_fixed.pdb', 'w'))



def clean_pdb_for_meeko(pdb_file):
    """Clean PDB file to avoid RDKit valence errors"""
    import re
    
    cleaned_file = pdb_file.replace('.pdb', '_cleaned_for_meeko.pdb')
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    cleaned_lines = []
    for line in lines:
        if line.startswith('ATOM'):
            # Fix common issues that cause valence errors
            # Remove alternative location indicators
            if len(line) > 16:
                line = line[:16] + ' ' + line[17:]
            
            # Fix occupancy and B-factor if they're problematic
            if len(line) > 54:
                # Set occupancy to 1.00 and reasonable B-factor
                line = line[:54] + '  1.00 20.00' + line[66:]
        
        cleaned_lines.append(line)
    
    with open(cleaned_file, 'w') as f:
        f.writelines(cleaned_lines)
    
    return cleaned_file


print("""

  How can I help you?
  1) seprate ligand and receptor
  2) pdb ligand to sdf file
  3) fix missing res and atoms
  4) clean for meeko

""")

option = input("Enter the number: ")

if int(option) == 1:
    pdb = input("Enter the pdb with .pdb ext: ")
    sep_lig_rec(pdb)
elif int(option) == 2:
    print("Attention! This should only be the ligand saved from pervieus step")
    lig = input("Enter the ligand name: ")
    fix_lig(lig)
elif int(option) == 3:
    print("Enter the pdb id without ext")
    pdb = input("Ehter the pdb id: ")
    fix_missing_res(pdb)
elif int(option) == 4:
    pdb = input("Enter the pdb file name: ")
    clean_pdb_for_meeko(pdb)

