import rdkit
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import combinations


#define a function to check whether the smiles str is correct
def valid_smls(smls):
    mol = Chem.MolFromSmiles(smls)
    if mol is not None:
        return True
    else:
        return False

#define a function to calculate TanimotoSimilarity between two smiles strs. 
def smilarity_cal(sml1,sml2):
    fp1 = Chem.RDKFingerprint(Chem.MolFromSmiles(sml1))
    fp2 = Chem.RDKFingerprint(Chem.MolFromSmiles(sml2))
    simiarity = rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)
    return(simiarity)

#main body

#read .smi file and delete wrong smiles (if any)
test_smls = []
with open('test_smiles.smi', 'r') as file:
    for line in file:
        parts = line.split()
        if valid_smls(parts[0]):
            test_smls.append((parts[1],parts[0]))
        else:
            print(f"No.{parts[1]} smiles str is wrong")        

#calculate the similarity score between every two compounds
#and pick out the two compounds which have the highest simlilarity score  
sim_list = []
for pair in combinations(test_smls,2):
    sim_score = smilarity_cal(pair[0][1],pair[1][1])
    sim_list.append((pair,sim_score))
most_sim = max(sim_list,key=lambda x:x[1])
print(f"""The two most similar compounds are: 
      No.{most_sim[0][0][0]}: {most_sim[0][0][1]}
      No.{most_sim[0][1][0]}: {most_sim[0][1][1]}
      The similarity score (TanimotoSimilarity) between them: {most_sim[1]}
      """)

#delete two most similar compounds and write the remaining compounds to a new .smi file
to_print_list = []
for items in test_smls:
    if items != most_sim[0][0] and items != most_sim[0][1]:
        to_print_list.append(str(items[1])+"  "+str(items[0]))
    else:
        print(f"No.{items[0]}: {items[1]} has been deleted")
base_filename = "filtered_smiles.smi"
index = 0
while True:
    if index == 0:
        new_filename = base_filename
    else:
        new_filename = f"filtered_smiles_{index}.smi"
    if not os.path.exists(new_filename):
        break
    index += 1
with open(new_filename, 'w') as file:
    for line in to_print_list:
        file.write(line + '\n')
print(f"Filtered smiles have been written to: {new_filename}")

