# Introduction to Variant Analysis

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | Shell for Bioinformatics | 4 - 2.5 hour sessions |

### Description

This ***Introduction to Variant Analysis*** workshop is aimed at providing best practices for calling variants for paired normal/tumor datasets. Importantly, while this workshop series does focus on calling variants in the context of paired tumor/normal samples, much of this workshop's pipeline and discussion is adaptable to other types of variant calling applications. This workshop series should provide participants with an ability to take raw sequence reads and process them into a VCF file with annotated variants. Furthermore, the workshop series ends with a tutorial for visualizing called variants within the Integrative Genomics Viewer (IGV). 

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.

#### Pre-Requisite

This workshop series is designed for users with a background in the fundamentals of working in the shell environment on an high-performance computing cluster (HPCC). If you don't have this background or are unsure if you are ready for this workshop series, you should first complete our [***Shell for Bioinformatics***](https://hbctraining.github.io/Shell-for-bioinformatics/) workshop series.

### Learning Objectives

- Evaluate QC metrics for variant calling
- Call variants using GATK
- Filter variants to retain only high-quality variant calls
- Annotate variants using SnpEff and dbSNP
- Prioritize variants by their impact
- Visualize variants in IGV
- Implement cBioPortal to analyze variants

### Lessons
* [Workshop schedule (trainer-led learning)](schedule/README.md)
* [Self-learning](schedule/self-learning.md)

### Installation Requirements
***All:***
- [FileZilla Client](https://filezilla-project.org/download.php?type=client) (make sure you get 'FileZilla Client')
- [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)

***Mac users:***
No additional installation requirements.  

***Windows users:***
[GitBash](https://git-scm.com/download/win)  


****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
