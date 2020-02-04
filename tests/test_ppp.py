import os
import PPP.main as mn


TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def test_ppp(pdb="1l63_proc.pdb", out="1l63_proc_processed.pdb"):
    mn.main(os.path.join(TEST_DIR,pdb), TEST_DIR)
    assert os.path.exists(os.path.join(TEST_DIR,out))
