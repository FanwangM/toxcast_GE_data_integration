from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import pickle
from numpy import zeros,array_str,set_printoptions,int
set_printoptions(threshold="nan")

# Path to folder containing the data
CALC_DATA_DIR='/scratch/az338/ucc-fileserver/ucc_az///toxCast_lincs_integration///data/'

# Import SMILES 
f = open(CALC_DATA_DIR+"chemicals.tab","r")
colnames = f.readline().strip("\n").split("\t")
chemicals = []
for l in f :
	chemicals.append(l.strip("\n").split("\t"))
	
f.close()

# Compute Morgan fingerprint
for c in chemicals :
	info = {}
	mol = Chem.MolFromSmiles(c[1])
	fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, bitInfo=info)
	arr = zeros((1,),dtype=int)
	DataStructs.ConvertToNumpyArray(fp,arr)
	c.append(arr)
	c.append(info)
	
# Write fingerprint and lincs_id to file to join to other datasets
f2 = open(CALC_DATA_DIR+"datasets/fingerprints.tab","w")
for c in chemicals :
	f2.write(c[0])
        for bit in c[-2] :
                f2.write("\t"+str(bit)) 
        f2.write("\n")
	
f2.close()

# Also save data to pickle object for bit retrieval/explanation
with open(CALC_DATA_DIR+'datasets/fingerprints.dat', 'wb') as f3:
      pickle.dump(chemicals, f3)


