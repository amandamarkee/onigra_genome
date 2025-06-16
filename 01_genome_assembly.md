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

## 3) Assembly Quality Control with BUSCO, assemblystats.py & QUAST
To assess genome quality, I use QUAST (or assemblystats.py) to obtain genome length statistics, and BUSCO for single-copy ortholog recovery (completeness). After assembly with hifiasm, we can assess assembly quality using the [assemblystats.py script](https://github.com/MikeTrizna/assembly_stats/tree/0.1.4) created by Mike Trizna.
- The version of assemblystats.py used here was modified by Paul Frandsen (Brigham Young University).

First, I copied this script into my working directory, and called it assemblystats.py

```
#!/usr/bin/env python

import numpy as np
from itertools import groupby
import json
import sys


def fasta_iter(fasta_file):
    """Takes a FASTA file, and produces a generator of Header and Sequences.
    This is a memory-efficient way of analyzing a FASTA files -- without
    reading the entire file into memory.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    header: str
        The string contained in the header portion of the sequence record
        (everything after the '>')
    seq: str
        The sequence portion of the sequence record
    """

    fh = open(fasta_file)
    fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in fa_iter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.upper().strip() for s in next(fa_iter))
        yield header, seq


def read_genome(fasta_file):
    """Takes a FASTA file, and produces 2 lists of sequence lengths. It also
    calculates the GC Content, since this is the only statistic that is not
    calculated based on sequence lengths.

    Parameters
    ----------
    fasta_file : str
        The file location of the FASTA file

    Returns
    -------
    contig_lens: list
        A list of lengths of all contigs in the genome.
    scaffold_lens: list
        A list of lengths of all scaffolds in the genome.
    gc_cont: float
        The percentage of total basepairs in the genome that are either G or C.
    """

    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []
    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        if "NN" in seq:
            contig_list = seq.split("NN")
        else:
            contig_list = [seq]
        for contig in contig_list:
            if len(contig):
                gc += contig.count('G') + contig.count('C')
                total_len += len(contig)
                contig_lens.append(len(contig))
    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    seq_array = np.array(seq_lens)
    stats['sequence_count'] = seq_array.size
    stats['gc_content'] = gc_cont
    sorted_lens = seq_array[np.argsort(-seq_array)]
    stats['longest'] = int(sorted_lens[0])
    stats['shortest'] = int(sorted_lens[-1])
    stats['median'] = np.median(sorted_lens)
    stats['mean'] = np.mean(sorted_lens)
    stats['total_bps'] = int(np.sum(sorted_lens))
    csum = np.cumsum(sorted_lens)
    for level in [10, 20, 30, 40, 50]:
        nx = int(stats['total_bps'] * (level / 100))
        csumn = min(csum[csum >= nx])
        l_level = int(np.where(csum == csumn)[0])
        n_level = int(sorted_lens[l_level])

        stats['L' + str(level)] = l_level
        stats['N' + str(level)] = n_level
    return stats


if __name__ == "__main__":
    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    contig_stats = calculate_stats(contig_lens, gc_cont)
    scaffold_stats = calculate_stats(scaffold_lens, gc_cont)
    stat_output = {'Contig Stats': contig_stats,
                   'Scaffold Stats': scaffold_stats}
    print(json.dumps(stat_output, indent=2, sort_keys=True))
```

Next, I changed permissions as follows to allow execution permissions.
```
chmod +x assemblystats.py
```

Lastly, I ran the assemblystats.py script on the newly generated fasta file of the fga in the format of scriptfilepath/scirptname.py assembly.fa
```
assemblystats.py assembly.fa
```

Save results as a text file as shown.
```
./assemblystats.py assembly.fa >> species-asmstats.txt
```

Alternatively, if you have the program QUAST installed on your computer or conda environment, you can use the following code instead, which will produce a directory
of more detailed plots on quality statistics:
```
module load quast 
quast assembly.fa 
```

<img width="515" alt="Screenshot 2025-05-25 at 1 31 48 PM" src="https://github.com/user-attachments/assets/7fdc2adb-86fe-45c2-a897-611b0ecacbc4" />

Note: These results were poor, and indicated that we needed to further run Purge Dups to decrease the duplicates and rid of short contigs. 
This is the first of three QUAST runs for QC: 1) Raw Assembly, 2) Post-Purge Dups, and 3) Post Decontamination 

For assessing single-copy ortholog recovery as a proxy for assembly completeness, I use the following BUSCO script with the orthodb10 Insecta dataset:

```
#!/bin/bash
#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10240M   # memory per CPU core
#SBATCH -J "onigra_busco"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# BUSCO_CONFIG_FILE=
source ~/.bashrc
source ~/.bash_profile
conda activate busco
busco -o busco_insecta-25 -i genome_assembly.fasta -l insecta_odb10 -c 24 -m genome #--offline
```

## 4) Duplicate Purging with PurgeDups
PurgeDups is a tool designed to remove haplotigs and contig overlaps in a de novo assembly based on read depth.
Below are the modified scripts I used from Dr. Paul Frandsen for execution of PurgeDups on the BYU cluster. For further information on PurgeDups, please see their [GitHub Repository](https://github.com/dfguan/purge_dups)

### Step 0: purge_dup_minimap.sh
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=16G  # memory per CPU core
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3/4.12-pws-472
conda activate minimap2

minimap2 -I6G -xasm20 -t 16 /home/fslcollab384/compute/genomics_workshop/onigra/onigra_hifiasm.asm.bp.p_ctg.fasta \
../m84100_240417_025302_s2.fastq.gz > ../m84100_240417_025302_s2.fastq.paf.gz
```

### Step 1: purge_dups_split_assembly.sh
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=32   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=16G  # memory per CPU core
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3/4.12-pws-472
conda activate purge_dups

assembly=../onigra_hifiasm.asm.bp.p_ctg.fasta

split_fa $assembly > $assembly.split && minimap2 -I6G -xasm5 -DP -t 32 $assembly.split $assembly.split | gzip -c - > $assembly.split.self.paf.gz
```

### Step 2: purge_dup_stats.sh
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=8G  # memory per CPU core
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3/4.12-pws-472
conda activate purge_dups

pbcstat ../*fastq.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs 2>calcults.log
```

### Step 3: purge_dups.sh
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=28G  # memory per CPU core
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3/4.12-pws-472
conda activate purge_dups

purge_dups -2 -T cutoffs -c PB.base.cov ../onigra_hifiasm.asm.bp.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
```

### Step 4: purge_dup_get_seqs.sh 
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --mem-per-cpu=8G  # memory per CPU core
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load miniconda3/4.12-pws-472
conda activate purge_dups

get_seqs dups.bed ../onigra_hifiasm.asm.bp.p_ctg.fasta
```

### Results from PurgeDups
<img width="419" alt="Screenshot 2025-05-25 at 1 30 27 PM" src="https://github.com/user-attachments/assets/d4232347-a8d0-4cdf-b2f8-1fa2af147481" />

## 5) Decontamination Filtering with BlobToolKit
BlobTools is a bioinformatics toolkit designed to help assess the quality and composition of genome assemblies. It visualizes metrics like GC content, read coverage, and taxonomic classification in interactive "blobplots" to identify potential contamination, symbionts, or misassemblies. This makes it especially useful for genome projects involving non-model organisms or environmental samples. For more information about the package, please see their [website](https://blobtools.readme.io/docs/the-blobtools-workflows#section-summary)

![image](https://github.com/user-attachments/assets/3bef6907-69fc-4c6a-9e27-3b6138852724)

Input Files for building BlobTools Database:
- blast hits file
- assembly fasta
- .bam of raw reads against assembly fasta

Below are the scripts I used to generate each of the input files:

### Step 1: blastn/megablast for whole genome

Note, because this is blasting the entire assembly against the NR database, one blast requires too much memory. Instead,
I run a blast array by first splitting the assembly as follows:

```
awk 'BEGIN {n_seq=0;} /^>/ \
{if(n_seq%50==0){file=sprintf("prefix_%d.fa",n_seq);} \
print >> file; n_seq++; next;} { print >> file; }' < [O_nigra_assembly].fasta
```

Then, I list the assembly split files (beginning with "prefix_" into a new file called blast.list:
```
ls prefix_* >> blast.list
```

Lastly, I use this blast.list file as input into the array job below:
```
#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=30   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=21240M   # memory per CPU core
#SBATCH -J "blob_blast"   # job name
#SBATCH --array=1-19


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
source ~/.bashrc
conda activate spidroins
module load blast-plus

i=$SLURM_ARRAY_TASK_ID
P=`awk "NR==$i" blast.list`

blastn \
-task megablast \
-db /apps/blast/databases/nt \
-query input/${P} \
-outfmt "6 qseqid staxids bitscore std" \
-max_target_seqs 20 \
-max_hsps 1 \
-evalue 1e-20 \
-num_threads $SLURM_NTASKS \
-out blast_out/${P}.blast.out
```

### Step 2: Asembly Fasta

I moved the assembly produced in the PurgeDups step prior, to my working directory for BlobTools.

### Step 3: Mapping (.bam) file with minimap2

Below is the script I use to generate the sorted .bam file for mapping raw reads back to the assembly.

```
#!/bin/bash
#SBATCH --time=25:00:00   # walltime
#SBATCH --ntasks=20   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=20144M   # memory per CPU core
#SBATCH -J "BlobTools minimap2"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# Set the max number of threads to use for programs using OpenMP. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source ~/.bashrc
conda activate spidroins
module load minimap
module load samtools

./minimap2-2.28_x64-linux/minimap2 \
-ax map-hifi \
-t $SLURM_NTASKS \
/home/fslcollab384/nobackup/archive/genomics_workshop/onigra/blobtools/onigra_purged_asm.fa \
/home/fslcollab384/nobackup/archive/genomics_workshop/onigra/RAW/m84100_240417_025302_s2.fastq.gz \
| samtools sort -@$SLURM_NTASKS -O BAM -o onigra_hifiasm_sorted.bam

samtools index onigra_hifiasm_sorted.bam
```

### Step 4: Running BlobTools to Identify Contaminants
```
#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=12288M   # memory per CPU core
#SBATCH -J "blobplots"   # job name
#SBATCH --mail-user=amarkee@amnh.org   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
source ~/.bashrc
conda activate spidroins

blobtools create -i /nobackup/archive/usr/fslcollab384/genomics_workshop/onigra/blobtools/onigra_purged_asm.fa \
-b /nobackup/archive/usr/fslcollab384/genomics_workshop/onigra/blobtools/onigra_hifiasm_sorted.bam \
-t /nobackup/archive/usr/fslcollab384/genomics_workshop/onigra/blobtools/onigra_blast-active.out \
-o onigra_blobplot-ACTIVE
```
The resulting file here will be in .json format, and used to filter contaminants later.

### Step 5: Making the Plots
```
conda activate blobtools
mkdir plots
blobtools plot -i [genome name]_blobplot.blobDB.json -o plots/
```

Output from this step will provide images of the contaminant plots, like this:
![onigra_blobplot-ACTIVE blobDB json bestsum phylum p8 span 100 blobplot bam0](https://github.com/user-attachments/assets/43fb8861-deaa-4aba-bdab-00f88b314cf6)

### Step 6: Remove Contaminants 
You can either use the built-in function in BlobTools called "seqfilter", or the function called "subseq". 
I was unable to get seqfilter to work again, so I created a list of contaminant contigs to use in seqtk as follows:

```
# create table from .json output. lists all contaminants and phyla level
blobtools view -i [blobDB in json format] -r [desired taxonomic level]
```

I isolated all the contigs I wanted to keep (i.e. contigs assigned to Arthropoda and no-hits) in Excel, because I was lazy today. Then, saved these contig names to a file called "contaminant_contigs-keep.txt". Lastly, I used the following
code to filter the assembly to only keep the listed contigs:

```
seqtk subseq original.fasta contigs-to-keep.txt > onigra_purged_asm.filtered.fasta
```
