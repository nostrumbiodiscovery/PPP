import os
import PPP.main as mn


TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def test_smiles_ppp(pdb="PR_1A28_xray_-_minimized.pdb", outpdb="PR_1A28_xray_-_minimized_processed.pdb",
    smiles="C2CCC1CCCCC1C2", ligand_pdb="STR.pdb"):
    if os.path.exists(os.path.join(TEST_DIR, outpdb)):
        os.remove(os.path.join(TEST_DIR, outpdb))
    mn.main(os.path.join(TEST_DIR, pdb), TEST_DIR, constrain_smiles=smiles, 
    ligand_pdb=os.path.join(TEST_DIR, ligand_pdb), output_pdb=[os.path.join(TEST_DIR, outpdb)])
    assert os.path.exists(os.path.join(TEST_DIR, outpdb))

def test_ppp_fix(pdb="1l63_proc.pdb", out="1l63_proc_processed.pdb"):
    if os.path.exists(os.path.join(TEST_DIR, out)):
        os.remove(os.path.join(TEST_DIR, out))
    mn.main(os.path.join(TEST_DIR, pdb), TEST_DIR, output_pdb=[os.path.join(TEST_DIR, out)])
    assert os.path.exists(os.path.join(TEST_DIR,out))


