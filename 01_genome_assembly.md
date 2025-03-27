# Workflow for Genome Assembly 

1) Sequencing and Raw Read Information (PacBio Revio)
2) Whole Genome Assembly with Hifiasm
3) Assembly Quality Control with BUSCO, assemblystats.py & QUAST
4) Haplotig duplicate purging with PurgeDups
5) Contamination filtering with BlobToolKit


## 1) Sequencing and Raw Read Information
Raw reads were acquired from Dr. Ed Wilcox and the Brigham Young University sequencing center in Spring of 2024. 

To generate reads, a single female _Oligotoma nigra_ from Dr. Janice-Edgerly Rooks was flash frozen and overnight shipped with dry ice to the BYU sequencing center.
Here, they extracted DNA with a mortar and pestle in liquid nitrogen, and column cleaned using Qiagen Genomic-tip columns with the recommended buffers purchased from the manufacturer.

Extracted DNA was DNA was prepared for circular consensus sequencing (CCS) using the SMRTbell Prep Kit 3.0. After library preparation, the concentration of the library DNA was measured using the 
Qubit DNA HS assay, and the DNA size was estimated using the Fragment Analyzer. These measurements were then entered into SMRT Link (PacBio) for each library. 

The target loading concentration was set at 500 pM to ensure most Zero-Mode Waveguide cells are loaded, with a 24-hour running time. Sequencing was carried out on a PacBio Revio system running SMRT Link Version 12.

```
ZMWs input                    : 16,326,794         

ZMWs pass filters             : 5,482,200 (33.58%)
ZMWs fail filters             : 10,844,594 (66.42%)
ZMWs shortcut filters         : 0 (0.000%)

ZMWs with tandem repeats      : 213,503 (1.969%)

Exclusive failed counts
Below SNR threshold           : 250,647 (2.311%)
Median length filter          : 1,297,272 (11.96%)
Lacking full passes           : 7,013,114 (64.67%)
Heteroduplex insertions       : 0 (0.000%)
Coverage drops                : 0 (0.000%)
Insufficient draft cov        : 0 (0.000%)
Draft too different           : 0 (0.000%)
Draft generation error        : 0 (0.000%)
Draft above --max-length      : 0 (0.000%)
Draft below --min-length      : 0 (0.000%)
Reads failed polishing        : 6,169 (0.057%)
Empty coverage windows        : 12,896 (0.119%)
CCS did not converge          : 0 (0.000%)
CCS adapter concatenation     : 2 (0.000%)
CCS adapter palindrome        : 7 (0.000%)
CCS adapter residue           : 61 (0.001%)
ZMW with full-length subread  : 2,263,392 (20.87%)
ZMW with control failure      : 930 (0.009%)
ZMW with control success      : 104 (0.001%)
CCS below minimum RQ          : 0 (0.000%)
Unknown error                 : 0 (0.000%)

Additional passing metrics
ZMWs missing adapters         : 484,090 (8.830%)

- - - - - - - - - - - - - - - : - - - - -

HiFi Reads                    : 5,227,102
HiFi Yield (bp)               : 85,702,860,819
HiFi Read Length (mean, bp)   : 16,395
HiFi Read Length (median, bp) : 15,594
HiFi Read Length N50 (bp)     : 15,594
HiFi Read Quality (median)    : 33
HiFi Number of Passes (mean)  : 8

<Q20 Reads                    : 255,098
<Q20 Yield (bp)               : 4,257,843,608
<Q20 Read Length (mean, bp)   : 16,691
<Q20 Read Length (median, bp) : 15,905
<Q20 Read Quality (median)    : 17

>=Q30 Reads                   : 3,398,649
>=Q30 Yield (bp)              : 52,469,028,161
>=Q30 Read Length (mean, bp)  : 15,438
>=Q30 Read Length (median, bp): 14,773
>=Q30 Read Quality (median)   : 36
```


## 2) Whole Genome Assembly with Hifiasm
To assemble HiFi reads, I ran the following script on BYU's Fulton HPC. Note, I switch between multiple HPC clusters during this project, and will note in each script where the analysis was run.

```
# Raw Reads in Directory on BYU Fulton
/home/fslcollab384/groups/fslg_nanopore/nobackup/archive/genomics_workshop_byu_may_24/amanda/o_nigra/m84100_240417_025302_s2_fastq.zip

# Directory for HiFiasm Results on BYU Fulton
/home/fslcollab384/nobackup/archive/genomics_workshop/onigra/hifiasm
```

```
# HiFiasm script used on BYU Fulton. Save as "hifiasm.job" 

#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=14336M   # memory per CPU core
#SBATCH -J "hifiasm_onigra"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load miniconda3/4.12-pws-472
conda activate hifiasm

hifiasm -o $1.asm -l 3 -t $SLURM_NPROCS $2

awk '$1 ~/S/ {print ">"$2"\n"$3}' $1.asm.bp.p_ctg.gfa > $1.asm.bp.p_ctg.fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' $1.asm.bp.hap1.p_ctg.gfa > $1.asm.bp.hap1.p_ctg.fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' $1.asm.bp.hap2.p_ctg.gfa > $1.asm.bp.hap2.p_ctg.fasta
```

Now, I run hifiasm to assemble the _O. nigra_ genome. The script requires ordered arguments in the following format:

```
sbatch hifiasm.job onigra m84100_240417_025302_s2_fastq.zip
```
