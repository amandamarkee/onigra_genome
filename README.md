# _De novo_ genome assembly and comparative genomics of black webspinner (_Oligotoma nigra_) silk genes

The black webspinner (_Oligotoma nigra_) is a non-holometabolous insect in the order Embioptera. Webspinners are known for their ample silk production throughout development, which they primarily use to protect sub-social or gregarious colonies within silken galleries. 
Here, we describe methods used to assemble, annotate and compare the primary silk genes of embiids (_e-fibroin_, a contraction for Embioptera fibroin). This analysis is in-part an effort to characterize the first whole genome for members of this order, and to 
assess interspecific variation in webspinner silk genes by comparing genomic archetecture and ensemble repeat composition to that of _Aposthonia ceylonica_, another member of the Embioptera (Family: Oligotomidae)

</br>
<img width="895" alt="Screenshot 2025-03-12 at 1 53 00 PM" src="https://github.com/user-attachments/assets/37df8e47-224c-4349-be02-04ed5c33079d" />
<br/><br/>

## Workflow for Genome Assembly
1) Whole Genome Assembly with Hifiasm
2) Assembly Quality Control with BUSCO, assemblystats.py & QUAST
3) Haplotig duplicate purging with PurgeDups
4) Contamination filtering with BlobToolKit

## Workflow for Whole Genome Annotation (Feature)
1) Repeat Element Modeling and Masking with EarlGrey
2) RNAseq Mapping with MiniMap2
3) Protein Database Aquisition with OrthoDB
4) Whole Genome Feature Annotation with BRAKER3

## Workflow for Silk Gene Manual Annotation
1) Putitive _e-fibroin_ Identification with BLASTx, BioPython and Excel
2) Fibroin Extraction with Samtools
3) Annotation of Intron/Exon Boundaries with IGV and Sequencher
