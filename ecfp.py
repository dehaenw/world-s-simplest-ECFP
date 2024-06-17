from rdkit import Chem, DataStructs
from numpy.linalg import matrix_power
import sys

def get_ecfp(mol,radius=3,fpsize=2048):
    invariants = []
    ecfp = DataStructs.ExplicitBitVect(fpsize)
    am = Chem.GetAdjacencyMatrix(mol,useBO=True)
    for i,a in enumerate(mol.GetAtoms()):
        am[i][i] = -1*a.GetAtomicNum()
    for r in range(radius+1):
        pm = matrix_power(am,r+1)
        invariants += [pm[i][i] for i in range(am.shape[0])]
    ecfp.SetBitsFromList([hash(inv)%fpsize for inv in invariants])
    return ecfp
    
if __name__ == "__main__":
    if len(sys.argv) > 1:
        mol = Chem.MolFromSmiles(sys.argv[1])
    else:
        mol = mol = Chem.MolFromSmiles("c1ncccc1C")
    print(list(get_ecfp(mol).GetOnBits()))
