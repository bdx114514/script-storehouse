from pdbfixer import PDBFixer
from openmm.app import PDBFile
def fixpdb(inputpdb,outputpdb):
    fixer = PDBFixer(filename=inputpdb)
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    keys = fixer.missingResidues.keys()
    for key in keys:
        chain = chains[key[0]]
        if key[1] == 0 or key[1] == len(list(chain.residues())):
            del fixer.missingResidues[key]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(outputpdb, 'w'))

from collections import defaultdict

import os
from multiprocessing import Pool

def findAllFile(base):
    for root, ds, fs in os.walk(base):
        for f in fs:
            yield root+'/'+f


pdb_dir = findAllFile("/home/chaos/L_database/bio_tmp/rcsb_simple_pdb/")
print(pdb_dir)

def run(_input):
    outf = _input[-8:-4]
    outf1 = 'fix_out_1/'+outf[1:3]+'/'+outf+'.pdb'
    os.system('mkdir -p '+'fix_out_1/'+outf[1:3])
    fixpdb(_input, outf1)
    return outf

with Pool(30) as p:
    try:
        print(p.map(run,list(pdb_dir)))
    except:
        print('error')






