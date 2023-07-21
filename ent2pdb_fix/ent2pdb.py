from Bio.PDB import PDBParser
import Bio.PDB.PDBIO as PDBIO
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
import os
import shutil

def ent2pdb(ents,outdir):
    p = PDBParser()
    path_out = outdir + '/' + os.path.split(ents)[0]
    n1 = os.path.basename(ents)
    n2 = os.path.splitext(n1)
    s = p.get_structure(n1,ents)
    model = s[0]
    ppb=PPBuilder()
    if len(ppb.build_peptides(model)) != 0:
        
        list_seq = []
        for i in range(len(ppb.build_peptides(model))):
            try:
                if len(ppb.build_peptides(model.get_list()[i])) != 0:
                    chaini = model.get_list()[i]
                    pepi = ppb.build_peptides(chaini)[0]
                    seq = pepi.get_sequence()
                    
                    if seq not in list_seq:
                        list_seq.append(seq)
                       
                        io = PDBIO()
                        io.set_structure(chaini)
                        filename = n2[0] + '_' + str(chaini.get_id()) + '.pdb'
                        io.save(filename)
                        print('pathout:' + path_out)
                        if not os.path.exists(path_out):
                            os.makedirs(path_out)
                        source = filename
                        destination = path_out
                        shutil.move(source, destination)
            except:
                pass
                continue


rootdir =  'pdb'
for parent, dirnames, filenames in os.walk(rootdir):
    for filename in filenames:
        paths = str(os.path.join(parent,filename))
        print('current_file:' + paths)
        if paths.endswith('.ent'):
            ent2pdb(paths, 'pdb_fixed')