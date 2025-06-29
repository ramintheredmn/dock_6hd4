import sys
from collections import namedtuple

Missing = namedtuple('Missing', ["ifmissing", "missing_ress", "num_per_chain"])
Residue = namedtuple('Residue', ["chain", "name", "num"])
Ligand = namedtuple('Ligand', ["name", "numatom", "chain", "syn"])

pdb = sys.argv[1]


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
        if len(chain)>1:
            chosen_chain = input(f"Choose the chain you want {chain}: ")
            print("you have chosen the chain, ", chosen_chain)

            with open(f"{sys.argv[1].split('.')[0]}_cleaned.pdb", 'w') as f:
                for line in lines:
                    if line.startswith('ATOM'):
                        if line.split()[4] == chosen_chain:
                            f.writelines(line)

        for ligand in ligands:
            if len(chain)>1:
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
