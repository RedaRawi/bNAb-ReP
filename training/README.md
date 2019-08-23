## Generate Training data

### Please install all dependencies (see main GitHub page)

### 3 Input arguments required
1. Training sequences in FASTA alignment file format.
2. Neutralization file: Each line contains either a 0 (resistant) or 1 (sensitive) for the corresponding sequences in the input sequences. 
3. Path to installed MAFFT software (e.g. "/usr/local/bin/mafft")

## Perform GBM model training
1. Number of cores (CPUs); Required for multiprocessing
2. Memory (in GB); Please provide number of GB in memory you want the training to allocate/use
3. (Full) Path to the current directory where the scripts are executed (output will be saved there) (e.g. "/home/rawir/bNAb-ReP")
4. (Full) Path to the training data, which was preprocessed in the previous step (e.g. "/home/rawir/bNAb-ReP/Training.txt")
