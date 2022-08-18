#!/usr/bin/env python
#!/usr/bin/python

"""
Created on Wed Jul 27 11:43:48 2022

@author: Suha Naser-Khdour
"""

import os
import re
import sys
import pickle
from tqdm import tqdm
from io import StringIO
from Bio import AlignIO
from Bio.Nexus import Nexus
from Bio.Align.Applications import MafftCommandline

def buscoAln():
    if len(sys.argv) == 1:
        print('Please provide the Busco directory and the sequence type.\nfor help use -h')
        arg0 = input("do you want to proceed? Y/N \t")
        if arg0 == 'Y' or arg0 == 'y':
            arg1 = input("Please provide the Busco directory\t")
            arg2 = input("Please provide the sequence type (DNA or AA)\t")
            arg3 = input("do you want to provide the output file name? Y/N. If No, the output will be BUSCO_DIRECTORY\extracted\MSA.nex\t")
            if arg3 == 'Y' or arg3 == 'y':
                arg4 = input("Please provide the output path\t")
                superAln(arg1,arg2,arg4)
            elif arg3 == 'N' or arg3 == 'n':
                superAln(arg1,arg2,os.path.join(arg1,'extracted','MSA.nex'))
            else:
                raise SystemExit('Error: unknown argument {}\n'.format(arg3))
        else:
            sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == '-h':
            print('\n')
            print('buscoAln.py creates a multi-sequence alignment for the busco genes of the different taxa\n')
            print('All the extracted busco multi sequences will be saved in BUSCO_DIRECTORY\extracted\ \n')
            print('Syntax:\n.')
            print('buscoAln.py <BUSCO_DIRECTORY> <SEQ_TYPE> <OUTPUT_ALN>\n')
            print('<BUSCO_DIRECTORY>  \t is the directory that contains all the taxa subdirectories.\n.')
            print('\t\t\t Make sure that the only subdirectories in BUSCO_DIRECTORY are the ones created by BUSCO.\n.')
            print('<SEQ_TYPE>  \t\t is the sequence type. AA or DNA.\n.')
            print ('<OUTPUT_ALN> \t\t Path where the final MSA will be saved. By default it will be saved in BUSCO_DIRECTORY\extracted\MSA.nex\n')
            print ('Note: All intermediate sequences of this code will be saved in BUSCO_DIRECTORY\extracted\ \n')
            return
        else:
            print('You need to provide the Busco directory and the sequence type.\nfor help use -h')
            arg0 = input("do you want to proceed? Y/N \t")
            if arg0 == 'Y' or arg0 == 'y':
                arg1 = input("Please provide the Busco directory\n")
                arg2 = input("Please provide the sequence type (DNA or AA)\n")
                arg3 = input("do you want to provide the output file name? Y/N. If No, the output will be BUSCO_DIRECTORY\extracted\MSA.nex\n")
                if arg3 == 'Y' or arg3 == 'y':
                    arg4 = input("Please provide the output path\n")
                    superAln(arg1,arg2,arg4)
                elif arg3 == 'N' or arg3 == 'n':
                    superAln(arg1,arg2,os.path.join(arg1,'extracted','MSA.nex'))
                else:
                    raise SystemExit('Error: unknown argument {}\n'.format(arg3))
            else:
                sys.exit(0)
    elif len(sys.argv) == 3:
        if (sys.argv[2] != 'AA') and (sys.argv[2] != 'DNA'):
            raise SystemExit('Error: unknown sequence type {}. Please choose AA or DNA.\n'.format(sys.argv[2]))
        else:
            superAln(sys.argv[1],sys.argv[2],os.path.join(sys.argv[1],'extracted','MSA.nex'))
    elif len(sys.argv) == 4:
        superAln(sys.argv[1],sys.argv[2],sys.argv[3])
    elif len(sys.argv) > 4:
        raise SystemExit('Error: too many arguments')

def superAln(wd, mode, outf):
    '''
    Parameters
    ----------
    wd : path
        the directory that contains all the species subdirectories. 
        This the same directory where all the Busco runs were saved.
    mode : str
        the mode, AA or DNA, for amino acid and nucleotide sequences, respectively.
    outf : path
        output file of the final multi-sequense alignment
    '''  
    os.chdir(wd)
    output = os.path.join(wd,'extracted') # Check if the output directory exists
    if not os.path.exists(output): #if output directory does not exist, create a directory
        os.mkdir(output)
    else: #if output directory existe, first clean the old output files
        for file in os.listdir(output):
            os.remove(os.path.join(output, file))
    with open('species.txt', 'rb') as f:
        species = pickle.load(f)
    extractBuscoPhylo(species,wd,mode)
    file_list = []
    for DirName, subdirList, FileList in os.walk(os.path.join(wd,'extracted')):
        for f in FileList:
            if 'trimmed' in f:
                file_list.append(os.path.join(wd,'extracted',f))
    msa = concat_without_missing_taxa(file_list, same_taxa=True)
    msa.write_nexus_data(filename=open(outf, 'w'))

def extractBuscoPhylo(species, wd, mode):
    '''
    Parameters
    ----------
    species : list
        list of species.
    wd : path
        the directory that contains all the species subdirectories. 
        This the same directory where all the Busco runs were saved.
    mode : str
        the mode, AA or DNA, for amino acid and nucleotide sequences, respectively.
    '''  
    # Keep only Busco genes that exist in all species
    complete = set([])
    for a_species in species:
        try:
            full_table = open('%s/%s/run_arthropoda_odb10/full_table.tsv' % (wd, a_species), 'r')
        except FileNotFoundError as e:
            continue
        complete_tmp = []
        for full_table_line in full_table:
            try:
                if 'Complete' in full_table_line.split()[1]: #  !!! >  if you don't want the best duplicated
                #if 'Complete' in full_table_line.split()[1] or 'Duplicated' in full_table_line.split()[1]: # if you want to include the best scoring duplicated. In our experience, it is less reliable, but we have not really benchmarked that.
                    complete_tmp.append(full_table_line.split()[0])
            except IndexError:
                pass
        if complete:
            complete = complete.intersection(set(complete_tmp))
        else:
            complete = set(complete_tmp)
    # Now using this list, retrieve all sequences, one file for each BUSCO. When duplicate, keep the best score
    outDir = os.path.join(wd,'extracted')
    for busco in tqdm(complete):
        if mode == 'AA':
            output_file = os.path.join(outDir,busco+'.faa')
        elif mode == 'DNA':
            output_file = os.path.join(outDir,busco+'.fna')
        else:
            print('Wrong mode specified. Please choose "AA" for amino acid sequences or "DNA" for DNA sequences')
        for a_species in species:
            seq = fetch(busco, '%s/%s/run_arthropoda_odb10' % (wd, a_species), mode)
            search_string = r'^>.*\n'
            seq = re.sub(search_string, '>%s\n' % a_species, seq)
            with open(output_file, 'a') as inf:
                inf.write('%s\n' % seq)
            trimAl(output_file, mode)

def fetch(busco_id, folder, mode):
    '''
    Parameters
    ----------
    busco_id : str
        the busco id
    folder : path
        the run folder
    mode : str
        the mode, AA or DNA, for amino acid and nucleotide sequences, respectively.
    Returns
    ----------
    sequence : the sequence
    '''
    if folder[-1] != '/':
        folder += '/'
    # load the full table
    full_table = {}
    full_table_file = open('%sfull_table.tsv' % (folder), 'r')
    for line in full_table_file:
        if not line.startswith('#'):
            if line.strip().split()[0] not in full_table:
                full_table.update({line.strip().split()[0]: [line.strip().split()[1:]]})
            else:
                full_table[line.strip().split()[0]].append(line.strip().split()[1:])
    if mode == 'AA':
        return _fetch_prot(busco_id, folder, full_table)
    elif mode == 'DNA':
        return _fetch_nucl(busco_id, folder, full_table)
    else:
        print('Wrong mode specified. Please choose "AA" for amino acid sequences or "DNA" for DNA sequences')
        
def _fetch_prot(busco_id, folder, full_table):
    '''
    Parameters
    ----------
    busco_id : str
        the busco id
    folder : path
        the run folder
    full_table : dict
        dictionary with all the unique busco genes (include the best scoring duplicate only)
    Returns
    ----------
    sequence : the best sequence for a BUSCO id, extracted from a protein BUSCO run     
    '''
    id_to_return = _fetch_best(busco_id, folder, full_table)
    # Now return the correct sequence
    sequences = open('%smetaeuk_output/combined_pred_proteins.fas' % folder, 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip().split(' ')[0] == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip().upper()
    return sequence

def _fetch_nucl(busco_id, folder, full_table):
    '''
    Parameters
    ----------
    busco_id : str
        the busco id
    folder : path
        the run folder
    full_table : dict
        dictionary with all the unique busco genes (include the best scoring duplicate only)
    Returns
    ----------
    sequence : the best sequence for a BUSCO id, extracted from a protein BUSCO run     
    '''
    id_to_return = _fetch_best(busco_id, folder, full_table)
    # Now return the correct sequence
    sequences = open('%smetaeuk_output/combined_nucl_seqs.fas' % folder, 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip().split(' ')[0] == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip().upper()
    return sequence

def _fetch_best(busco_id, folder, full_table):
    '''
    Parameters
    ----------
    busco_id : str
        the busco id
    folder : path
        the run folder
    full_table : dict
        dictionary with all the unique busco genes (include the best scoring duplicate only)  
    Returns
    ----------
    id_to_return : str
        the sequence id to retrieve, with the best score
    '''
    result_scores = []
    good_results = []
    for entry in full_table[busco_id]:
        result_scores.append(float(entry[5]))
    # open all hmm result file for this busco id
    for file in os.listdir('%shmmer_output/initial_run_results/' % folder):
        if file.startswith(busco_id):
            for line in open('%shmmer_output/initial_run_results/%s' % (folder, file), 'r'):
                if not line.startswith('#'):
                    if float(line.split()[7]) in result_scores:
                        good_results.append(line)
    if not good_results:
        for file in os.listdir('%shmmer_output/rerun_results/' % folder):
            if file.startswith(busco_id):
                for line in open('%shmmer_output/rerun_results/%s' % (folder, file), 'r'):
                    if not line.startswith('#'):
                        if float(line.split()[7]) in result_scores:
                            good_results.append(line)
    # Recheck that the best result has the longest protein, should be the best score as well. 
    # Edit: in fact no, for transcriptomes it is not since the 6 frames translation are evaluated and only one is correct, so the warning below are not really useful, don't focus too much on it.
    id_to_return = None
    best_score = 0
    for result in good_results:
        if float(result.split()[7]) > best_score:
            id_to_return = result.split()[0]
    return id_to_return

def trimAl(infile, mode):
    '''
    Parameters
    ----------
    infile : path
        the multi sequence file
    mode : str
        the mode, AA or DNA, for amino acid and nucleotide sequences, respectively.
    '''
    aln = AlignMafft(infile, mode)
    if mode == 'AA':
        outf = '%s_trimmed.faa' % os.path.splitext(infile)[0]
    if mode == 'DNA':
        outf = '%s_trimmed.fna' % os.path.splitext(infile)[0]
    runTrimal = " ".join(['trimal -in', aln, '-out', outf, '-nexus -gappyout -keepheader'])
    os.system(runTrimal)
    return

def AlignMafft(infile, mode):
    '''
    Parameters
    ----------
    infile : path
        the multi sequence file
    mode : str
        the mode, AA or DNA, for amino acid and nucleotide sequences, respectively.
    Returns
    ----------
    outf : path
        the multi sequence alignment
    '''
    if mode == 'AA':
        outf = '%s_aligned.faa' % os.path.splitext(infile)[0]
    elif mode == 'DNA':
        outf = '%s_aligned.fna' % os.path.splitext(infile)[0]
    mafft_cline = MafftCommandline(input=infile)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(StringIO(stdout), "fasta")
    AlignIO.write(align,outf, "fasta")
    return outf

def check_taxa(matrices):
    """Verify Nexus instances have the same taxa information.

    Checks that nexus instances in a list [(name, instance)...] have
    the same taxa, provides useful error if not and returns None if
    everything matches
    """
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]
        if first_only:
            missing = ', '.join(first_only)
            msg = '%s taxa %s not in martix %s' % (matrices[0][0], missing, name)
            raise Nexus.NexusError(msg)
        elif new_only:
            missing = ', '.join(new_only)
            msg = '%s taxa %s not in all matrices'  % (name, missing)
            raise Nexus.NexusError(msg)
    return None # will only get here if it hasn't thrown an exception

def concat_without_missing_taxa(file_list, same_taxa=True):
    """Combine multiple nexus data matrices in one partitioned file.

    By default this will only work if the same taxa are present in each file
    use same_taxa=False if you are not concerned by this
    """
    nexi = [(os.path.basename(fname), Nexus.Nexus(fname)) for fname in file_list]
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)

if __name__ == '__main__':
    buscoAln()
