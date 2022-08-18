#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:21:49 2022

@author: u1044662
"""

import os
import gzip
import sys
import shutil
import pickle
import tempfile

def BuscoCode():
    if len(sys.argv) == 1:
        raise SystemExit('\nError: Please provide the assemblies directory\n')
    if len(sys.argv) == 2:
        if sys.argv[1] == '-h':
            print('\nSyntax:\n')
            print('BuscoCode.py <ASSEMBLIES_DIR>\n')
            print('<ASSEMBLIES_DIR> The path to all the assemblies.\n.')
            print('\nNotes:\n\t\t *The name of each assembly should be the taxon name that you want to appear in the final alignment.\n')
            print('\t\t *Assemblies should be in fasta format.\n')
            print('\t\t *BuscoCode accepts gzip-compressed files.\n')
            return
        else:
            if os.path.exists(sys.argv[1]):
                species = []
                for DirName, subdirList, FileList in os.walk(sys.argv[1]):
                    if FileList:
                        for file in FileList:
                            species.append(os.path.basename(file).split('.')[0])
                            outDir = runBusco(DirName,file)
                with open('Alignments/species.pickle', 'wb') as fp:
                    pickle.dump(species, fp)
                logf = open(os.path.join(outDir, 'logs', 'busco.log'), 'r')
                if 'BUSCO analysis done.' in logf:
                    print('\n{} Busco folders finished succefuly for taxa {}\n'.format(len(species), species))
                    print('All taxa names can be found in ./Alignments/species.pickle \n')
            else:
                raise SystemExit('Error: No such directory {}'.format(sys.argv[1]))
    elif len(sys.argv) > 2:
        raise SystemExit('Error: too many arguments')

def runBusco(Dirname,file):
    #This code calls Busco vis os.system
    '''
    Parameters
    ----------
    file: str
         Assembly files in fasta format.
    Dirname : str
        the directory where the assembly files are stored.
    '''
    wd = os.getcwd()
    out = r'Alignments/' + os.path.basename(file).split('.')[0]
    if not os.path.exists(os.path.join(wd,out)):
        if is_gz_file(os.path.join(Dirname,file)):
            with gzip.open(os.path.join(Dirname,file), 'rb') as inf:
                with tempfile.NamedTemporaryFile() as tmp:
                    shutil.copyfileobj(inf, tmp)
                    bus = " ".join(['busco -i',tmp.name,'-o',out,'-m genome --auto-lineage-euk'])
                    os.system(bus)
            tmp = tempfile.NamedTemporaryFile(delete=True)
        else:
            bus = " ".join(['busco -i',os.path.join(Dirname,file),'-o',out,'-m genome --auto-lineage-euk'])
            os.system(bus)
    else:
        force = input("\nDirectory {} already exist. Do you want to rewrite existing files? Y/N \t".format(os.path.join(wd,out)))
        if force == 'Y' or force == 'y':
            if is_gz_file(os.path.join(Dirname,file)):
                with gzip.open(os.path.join(Dirname,file), 'rb') as inf:
                    with tempfile.NamedTemporaryFile() as tmp:
                        shutil.copyfileobj(inf, tmp)
                        bus = " ".join(['busco -i',tmp.name,'-o',out,'-m genome --auto-lineage-euk -f'])
                        os.system(bus)
                tmp = tempfile.NamedTemporaryFile(delete=True)
            else:
                bus = " ".join(['busco -i',os.path.join(Dirname,file),'-o',out,'-m genome --auto-lineage-euk -f'])
                os.system(bus)
    return os.path.join(wd,out)

def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'
    
if __name__ == '__main__':
    BuscoCode()
