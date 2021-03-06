# ENSP_crystal_structure_match  
###  Before processing and downloading the sequences from Ensembl:
###  Remember to convert 'X' , '*' and 'U' to 'C'  
```  
$ sed -i 's/X/C/g' *
$ sed -i 's/*/C/g' *
$ sed -i 's/U/C/g' *
```  
###  Library packages needed:  
``` 
pandas, biopython, python2.7 
# import pandas,  numpy,  Bio,  re  
```  
## local PC version  
### Instructions:  
  ```
 1. the 'SNPs_with_crystal_structure' contains all proteins need to be process  
 2. 'pdb_only/ENSPxxxxxxxxxx.out' contains a list of crystal structure, which will compare with ENSPxxxxxxxxxx.seq which can be found in '/test/'  
 3. all the sequence of crystal structure are located in 'fasta/'  
 4. run file : multi_mapping.py    
 ```
## HPC version
**_hpc/_ folder contains gnuparallel method processing on SupreMike**   

```
split SNPs_with_crystal_structure into several part(every little part has titile), feed them to cores coorespondently.  
0. load module python, gcc, gnuparallel all you needed in ~/.modules
    ~/.modules : module load python/2.7.14-anaconda
                 module load gnuparallel/20180222/INTEL-18.0.0
    $ conda create -n your_env python=2.7
    $ source activate your_env
    $ conda install pandas
    $ conda install biopython
   
   ***** all the path needed in scripts should be changed to our own path ******
   
1. 'python gen_input_lst.py -numLst N' sepreate SNPs_with_crystal_structure into N parts  
    (N % 16 ==0, hpc each node has 16 cores )  
    and generate 'input.lst' which will feed to 'pbs.script' 
2. 'qsub pbs.script' to hpc, the number of nodes applied should be N/16. 
   (eg. N=16 ~ 1node, N=48 ~ 3nodes, N=128 ~ 8nodes)  
   (8 nodes: 45 mins on supermike)
3. 'python combine_folder.py' can combine N part into one file.
  
```
