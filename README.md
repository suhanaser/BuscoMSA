# BuscoMSA
A repository that contains python scripts to create a multi-sequense alignment from busco genes .
## Scripts
### 1. BuscoCode.py
This script takes a directory of all the assemblies and creates a Busco directory for each assembly. It requires python 3.10.x or higher and dependencies as in the header of the script.

Note: This script calls Busco software via `os.system()`. Therefore, first you need to install Busco (see instructions [here](https://busco.ezlab.org/busco_userguide.html)) and if neseccary adjust the script to call busco from the relevant path.

#### Running the script
You can either download just the `BuscoCode.py` script from this repo, or you can get the whole repo. It's proabably easiest to do the latter, like this:
```
git clone git clone https://github.com/suhanaser/BuscoMSA.git
cd BuscoMSA
chmod +x BuscoCode.py
```
Now you can run the script as follows
```
python BuscoMSA/BuscoCode.py
```
To get the help information, just do this:
```
BuscoMSA/BuscoCode.py -h
```
That should show you the following help:
```
Syntax:
BuscoCode.py <ASSEMBLIES_DIR>
<ASSEMBLIES_DIR>  The path to all the assemblies 

Notes: 
  * The name of each assembly should be the taxon name that you want to appear in the final alignment.
  * Assemblies should be in fasta format.
  * BuscoCode accepts gzip-compressed files.
```
This script should produce a directory called `Alignments` which contains subdirectories for each taxon and `species.pickle` file that contains the list of all taxa. 

### 2. buscoAln.py
This script takes as input the directory that contains the `BuscoCode.py` output. 
if you didn't modify the `BuscoCode.py`script, the name of this directory should be `/Alignments`.

 In addition, you should specify the type of sequences that you want to align (AA or DNA).
 
 Optional: you can give the sctipt the path to the output MSA.
 
 #### Running the script
 If you cloned the whole BuscoMSA repo, you can run the script as follows
```
python BuscoMSA/buscoAln.py
```
To get the help information, just do this:
```
python BuscoMSA/BuscoCode.py -h
```
The syntax for this script is:
```
buscoAln.py <BUSCOS_DIR> <SEQ_TYP> <OUTPUT>
<BUSCOS_DIR>  is the directory where all the busco sequenses are saved. This is the  `/Alignments` directory output from running `BuscoCode.py` script.
<SEQ_TYP>     AA or DNA
<OUTPUT>      the output multi-sequence alignment (optional). If not provided the output will be `BUSCO_DIR\extracted\MSA.nex`.
```
