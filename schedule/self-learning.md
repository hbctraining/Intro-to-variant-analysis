# Introduction to Variant Analysis

## Learning Objectives

- Evaluate QC metrics for variant calling
- Call variants using GATK
- Filter variants to retain only high-quality variant calls
- Annotate variants using SnpEff and dbSNP
- Prioritize variants by their impact
- Visualize variants in IGV

## Installations

### On your desktop

1. [FileZilla Client](https://filezilla-project.org/download.php?type=client) (make sure you get ‘FileZilla Client')
2. [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)

### On your HPCC (if not using Harvard's O2 cluster)

#### Required
1. [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) version 0.11.9
2. [`bwa`](https://bio-bwa.sourceforge.net) version 0.7.17
3. [`Picard`](https://broadinstitute.github.io/picard/) version 2.27.5
4. [`MultiQC`](https://multiqc.info) version 1.12
5. [`GATK`](https://gatk.broadinstitute.org/hc/en-us) version 4.1.9.0
6. [`SnpEff and SnpSift suite`](http://pcingola.github.io/SnpEff/) version 4.3g
7. [`bcftools`](http://www.htslib.org) version 1.14

#### Optional
1. [`samtools`](http://www.htslib.org) version 1.15.1
2. [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) version 2.30.0

> ***NOTE:*** If you are not working on the O2 cluster and are using different versions of these software programs, these packages may still work with the provided commands. However,this workshop was designed on these versions specifically, so you may need to tweak some of the commands if you use different versions of this software.

## Lessons

1. [ICGC-TCGA DREAM Mutation Calling Challenge Synthetic Dataset](../lessons/01_syn3_dataset.md)
2. [Project Organization](../lessons/02_project_organization.md)
3. [File Formats](../lessons/03_file_formats.md)
4. [Evaluating Read Quality with `FastQC`](../lessons/04_fastqc.md)
5. [Sequence Read Alignment](../lessons/05_sequence_alignment_theory.md)
6. [Alignment File Processing ](../lessons/06_alignment_file_processing.md)
7. [Alignment File Quality Control](../lessons/07_alignment_QC.md)
8. [Evaluating Quality Control Metrics](../lessons/08_evaluate_QC.md)
9. [Variant Calling](../lessons/09_variant_calling.md)
10. [Variant Filtering](../lessons/10_variant_filtering.md)
11. [Variant Annotation with SnpEff](../lessons/11_variant_annotation.md)
12. [Automation of Variant Calling Pipeline](../lessons/12_automation_of_variant_calling.md)
13. [Variant Prioritization with SnpSift](../lessons/13_variant_prioritization.md)
14. [Visualization in IGV](../lessons/14_IGV.md)

> ***NOTE:*** If you aren't working on Harvard's O2 cluster the directory structure for the HPCC that you are using is likely different and you will need to modify paths to work within your HPCC's directory structure.

## Answer key
* [Exercises Answer Key](../answer_key/Answer_key.md)

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
