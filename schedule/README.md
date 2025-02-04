# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Shell for Bioinformatics](https://hbctraining.github.io/Shell-for-bioinformatics/schedule/) workshop is a prerequisite.


## Pre-reading:

* Please **study the contents** within the following lessons:
  * [Shell basics review](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/shell_review.html)
  * [Best Practices in Research Data Management (RDM)](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon-flipped/lessons/04a_data_organization.html)
  * [Introduction to Variant Calling](../lessons/00_intro_to_variant_calling.md) 


## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:10 | [Workshop Introduction](../lectures/workshop_intro_slides.pdf) | Will |
| 10:00 - 11:30 | [Introduction to Variant Calling](../lectures/Variant_calling_intro.pdf) | Elizabeth |
| 11:30 - 11:50 | [Project Organization](../lessons/01_data_organization.md) | Elizabeth |
| 11:50 - 12:00 | Overview of self-learning materials and homework submission | Will |

### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:

  1. [Evaluating Read Quality with `FastQC`](../lessons/02_fastqc.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>The first step in many NGS studies is first to evaluate the read qualites that you received from the sequencing facility. A common tool used for handling this analysis is <code>FastQC</code>. 
         <br><br>This lesson will:<br>
          <ul><li>Implement <code>FastQC</code> to evaluate read qualities</li>
          <li>Evaluate FASTQC quality metrics</li></ul>
          <hr />
        </details>

  2. [Sequence Read Alignment](../lessons/03_sequence_alignment_theory.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Once we have completed our QC on sequence reads we will be aligning the reads to a reference sequence. This alignment step places each read in genomic space and creates the bedrock for calling variants.
         <br><br>This lesson will:<br>
             <ul><li>Enumerate difficulties with alignment</li>
             <li>Create an <code>sbatch</code> script to align reads</li></ul>
             <hr />
        </details>

  3. [Alignment File Processing ](../lessons/04_alignment_file_processing.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Before we can call variants from our alignment files, we need to do some processing to clean up the alignment files. The two major concerns here are organizing (sorting) our alignment files for our analyses and removing duplicates.
         <br><br>This lesson will:<br>
             <ul><li>Differentiate between query-sorted and coordinate-sorted alignment files</li>
             <li>Describe and remove duplicate reads</li>
             <li>Process a raw SAM file for input into a BAM for <code>GATK</code></li></ul>
             <hr />
        </details>


> **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it).
> 1. Log in using `ssh rc_trainingXX@o2.hms.harvard.edu` and enter your password (replace the "XX" in the username with the number you were [assigned in class](https://docs.google.com/spreadsheets/d/1kBlYowhjjHJC9ZovmbBULmbqozKpprM17vZ2wPlhNg0/edit?usp=sharing)). 
> 2. Once you are on the login node, use `srun --pty -p interactive -t 0-2:30 --mem 1G /bin/bash` to get on a compute node.
> 3. Proceed once your command prompt has the word `compute` in it.
> 4. If you log out between lessons (using the `exit` command twice), please follow points 1. and 2. above to log back in and get on a compute node when you restart with the self learning.

II. **Complete the exercises**:
   * Each lesson above contains exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form](https://forms.gle/N6eRvvQYhGxLkU837) the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

***

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:00 | Self-learning lessons review | All |
| 10:00 - 10:30 | [Alignment File Quality Control](../lessons/05_alignment_QC.md) | Elizabeth |
| 10:30 - 10:40 | Break |  |
| 10:40 - 11:15 | [Aggregating QC metrics using MultiQC](../lessons/06_aggregate_multiqc.md) | Elizabeth |
| 11:15 - 12:00 | [Variant Calling](../lessons/07_variant_calling.md) | Will |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:

1. [Variant Filtering](../lessons/08_variant_filtering.md) 

      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Now that we have called our raw variants, we will need to filter our data for only high-quality variant calls. Low-quality variant calls can occur for a variety of reasons that we will explore and we will implement steps to exclude them.
         <br><br>This lesson will:<br>
          <ul><li>Filter raw variant calls using <code>FilterMutectCells</code> to reduce errors</li>
          <li>Remove Low-Complexity Regions from the called variants using <code>SnpSift</code> to further reduce errors</li></ul>
        </details>

2. [Variant Annotation with SnpEff](../lessons/09_variant_annotation.md) 

      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>With our high-quality variant calls, we would like to know more information about these variants. For example, we might like to know which genes our they are in or how they alter the protein-coding sequence for the genes they are in. In order to do this, we will need to provide annotations for our genes.
         <br><br>This lesson will:<br>
          <ul><li>Annotate a VCF file for functional impacts with `SnpEff`</li>
          <li>Differentiate between an unannotated and annotated VCF file</li></ul>
          <hr />
        </details>

> **NOTE:** To run through the code above, you will need to be **logged into O2** and **working on a compute node** (i.e. your command prompt should have the word `compute` in it). For login instructions, please see above.

II. **Complete the exercises**:
   * Each lesson above contains exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form](https://forms.gle/WMUnZHsJT3dBiBHi7) the **day before the next class**.
   
### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:00 | Self-learning lessons review | Elizabeth |
| 10:00 - 10:30 | [Variant Prioritization with SnpSift](../lessons/10_variant_prioritization.md) | Elizabeth |
| 10:30 - 11:00 | [Exercise](../lessons/11_SnpSift_exercise.md) ([Key](../lessons/SnpSift_exercise_KEY.md)) | Will |
| 11:00 - 11:30 | [Visualization in IGV](../lessons/12_IGV.md) | Will |
| 11:30 - 12:00 | Q & A (review of Automation) | All |

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

***

## Day 4

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 9:30 - 10:30 | Introduction to cBioPortal | Dr. Tali Mazor |
| 10:30 - 11:30 | cBioPortal Practical | Dr. Tali Mazor |
| 11:30 - 11:45 | [Oncoprint Integration](../lessons/13_oncoprint_creation.md) | Will |
| 11:45 - 12:00 | [Wrap up](../lectures/workshop_wrapup_slides.pdf) | Elizabeth |

## File Format Reference
* [File Formats](../lessons/file_formats_reference.md)

## Automation Reference
* [Automation of Variant Calling Pipeline](../lessons/automation_of_variant_calling.md)

## Answer key
* [Exercises Answer Key](../answer_key/Answer_key.md)

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Some materials used in these lessons were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
