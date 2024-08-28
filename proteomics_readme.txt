Each file contains association summary statistics for one of the 4,907 aptamers. 
The files are named using the following convention: SeqId_GeneName_ProteinName.txt.gz where SeqId uniquely identifies the aptamer.

These summary statisics files have the following columns:
Chrom: Chromosome
Pos: Position (hg38)
Name: Unique variant name
rsids: rs-name, if it exists
effectAllele: Effect allele
otherAllele: Non-effect allele
Beta: Effect (in standard deviations)
Pval: P-value
min_log10_pval: -log10 of P-value
SE: Standard error
N: Sample size
ImpMAF: Minor allele frequency.

Note that in the summary statistics files, the effectAllele is not always the minor allele, and therefore the ImpMAF is not always the frequency of the effect allele (it may be the frequency of the other allele). 
We now provide an Extra annotation file (assocvariants.annotated.txt.gz) that also includes the effect allele frequency and has the following columns:
Chrom: Chromosome
Pos: Position (hg38)
Name: Unique variant name
rsids: rs-name, if it exists
effectAllele: Effect allele
otherAllele: Non-effect allele
effectAlleleFreq: Effect allele frequency.
This file can be mapped to the other files using the "Name" column. 
Also, the summary statistics files sometimes incorrectly have effectAllele=otherAllele for multiallelic variants. 
In these cases the effectAllele is correct, but the otherAllele should be '!', meaning that the effectAllele is tested against the other (two or more) alleles (using the '!' sign as shorthand for 'not effectAllele'). 
This has been corrected in the file Extra annotation file (assocvariants.annotated.txt.gz). 

Finally, a subset of the variants in the summary statistics files should be excluded due to quality issues.
These variants are listed in a separate Excluded variants file (assocvariants.excluded.txt.gz). 
The Extra annotation file (assocvariants.annotated.txt.gz) does not include these variants.
