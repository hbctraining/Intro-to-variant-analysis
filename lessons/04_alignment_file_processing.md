---
title: "Alignment File Processing"
author: "Will Gammerdinger, Meeta Mistry"
date: "December 11, 2023"
---

Approximate time: 45 minutes

## Learning objectives

- Differentiate between query-sorted and coordinate-sorted alignment files
- Describe and remove duplicate reads
- Process a raw SAM file for input into a BAM for GATK

## Alignment file processing

The alignment files that come from `bwa` are raw alignment and need some processing before they can be used for variant calling. While all of the data is present, we need to format it in a way that helps the variant calling algorithm process it. This process is very similar to flipping all of the pieces of a puzzle over to the correct side and grouping edge pieces or pieces by a similar color or pattern when starting a jigsaw puzzle. 

<p align="center">
<img src="../img/Process_alignment_workflow.png" width="800">
</p>


## Pipeline for processing alignment files with GATK

[Picard](https://broadinstitute.github.io/picard/) is a set of command line tools for **processing high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF**. It is maintained by the Broad Institute, and is open-source under the MIT license and free for all uses. The Broad Institute also maintains the tool [GATK](https://gatk.broadinstitute.org/hc/en-us), which we will be using for variant calling and discuss more later on. The all of the functions of Picard have been ported over to GATK, so instead of using Picard, we will just use Picard's functions which are now a part of GATK. To see all of the function within GATK you can use:

```
gatk --list
```

You will notice that some packages have `(Picard)` after them and these represent tools that have been brought over from Picard and incorporated into GATK now. As a result of this history and subsequent merging of the tools, you may see people reference GATK or Picard for some of these functions interchangeably. Picard tools were written without multithreading and thus when they were brought over to GATK, many, if not all, of them retained their inability to multithread.

> ### Why not use `samtools`?
> The processing of the alignment files (SAM/BAM files) can also be done with [`samtools`](https://github.com/samtools/samtools). While there are some advantages to using samtools (i.e. more user-friendly, multi-threading capability), there are slight formatting differences which we may want to take advantage of downstream. Since we will be using GATK later in this workshop (also from the Broad Institute), Picard seemed like a more suitable fit.
>
> **Near the end of the lesson, there will be a dropdown, if you would like to know how to do the processing steps in `samtools`.**

> **NOTE:** You may encounter a situation where your reads from a single sample were sequenced across different lanes/machines. As a result, each alignment will have different read group IDs, but the same sample ID (the SM tag in your SAM/BAM file). You will need to merge these files here before continuing. The dropdown menu below will detail how to do this.
>
><details>
> <summary><b>Click here if you need to merge alignment files from the same sample</b></summary>
>   You can merge alignment files with different read group IDs from the same sample in both <code>Picard</code> and <code>samtools</code>. In the dropdowns below we will outline each method:
> <details>
>   <summary><b>Click to see how to merge SAM/BAM files in <code>GATK</code></b></summary>
>   First, we need to load the <code>gatk</code> module:
>   <pre>
>   module load gatk/4.6.1.0</pre>
>   We can define our variables as:
>   <pre>
>   INPUT_BAM_1=Read_group_1.bam
>   INPUT_BAM_2=Read_group_2.bam
>   MERGED_OUTPUT=Merged_output.bam</pre>
>   Here is the command we would need to run to merge the SAM/BAM files:
>   <pre>
>   gatk MergeSamFiles \
>     --INPUT $INPUT_BAM_1 \
>     --INPUT $INPUT_BAM_2 \
>     --OUTPUT $MERGED_OUTPUT</pre>
>   We can breakdown this command:
>   <ul><li><code>gatk MergeSamFiles</code> This calls the <code>MergeSamFiles</code> from within <code>gatk</code></li>
>   <li><code>--INPUT $INPUT_BAM_1</code> This is the first SAM/BAM file that we would like to merge.</li>
>   <li><code>--INPUT $INPUT_BAM_2</code> This is the second SAM/BAM file that we would like to merge. We can continue to add <code>--INPUT</code> lines as needed.</li>
>   <li><code>--OUTPUT $MERGED_OUTPUT</code> This is the output merged SAM/BAM file</li></ul>
>   <hr />
> </details>
> <details>
>   <summary><b>Click to see how to merge SAM/BAM files in <code>samtools</code></b></summary>
>   First, we need to load the <code>samtools</code> module, which also requires <code>gcc</code> to be loaded:
>   <pre>
>   module load gcc/14.2.0
>   module load samtools/1.21</pre>
>   We can define our variables as:
>   <pre>
>   INPUT_BAM_1=Read_group_1.bam
>   INPUT_BAM_2=Read_group_2.bam
>   MERGED_OUTPUT=Merged_output.bam
>   THREADS=8</pre>
>   Here is the command we would need to run to merge the SAM/BAM files:
>   <pre>
>   samtools merge \
>     -o $MERGED_OUTPUT \
>     $INPUT_BAM_1 \
>     $INPUT_BAM_2 \
>     --output-fmt BAM \
>     --threads $THREADS</pre>
>   We can break down this command:
>   <ul><li><code>samtools merge</code> This calls the <code>merge</code> package within <code>samtools</code>.</li>
>   <li><code>-o $MERGED_OUTPUT</code> This is the merged output file.</li>
>   <li><code>$INPUT_BAM_1</code> This is the first SAM/BAM file that we would like to merge.</li>
>   <li><code>$INPUT_BAM_2</code> This is the second SAM/BAM file that we would like to merge. We can continue to add additional input SAM/BAM files to this list as needed.</li>
>   <li><code>--output-fmt BAM</code> This specifies the output format as <code>BAM</code>. If for some reason you wanted a <code>SAM</code> output file then you would use <code>--output-fmt SAM</code> instead.</li>
>   <li><code>--threads $THREADS</code> This specifies the number of threads we want to use for this process. We are using 8 threads in this example, but this could be different depending on the parameters that you would like to use.</li></ul>
> </details>
></details>

Before we start processing our alignment SAM file with `GATK`/`Picard`, let's take a quick look at the steps involved in this pipeline. 

<p align="center">
<img src="../img/Process_alignment_workflow_zoom_in.png" width="600">
</p>

Several goals need to be accomplished, and we will go through each in detail in this lesson!

Let's begin by creating a script for alignment processing. Make a new `sbatch` script within `vim`:

```
$ cd ~/variant_calling/scripts/
$ vim gatk_alignment_processing_normal.sbatch
```

As always, we start the `sbatch` script with our shebang line, description of the script and our `sbatch` directives to request appropriate resources from the O2 cluster. 

```
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 1
#SBATCH --mem 32G
#SBATCH -o gatk_alignment_processing_normal_%j.out
#SBATCH -e gatk_alignment_processing_normal_%j.err
```

Next we load the `GATK` module: 

```
# Load module
module load gatk/4.6.1.0
```

**Note: `GATK` is software that does NOT require gcc/14.2.0 to also be loaded** 

Next, let's define some variables that we will be using:

```
# Assign file paths to variables
SAMPLE_NAME=syn3_normal
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_GRCh38.p7.sam
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/gatk/${SAMPLE_NAME}/
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`
```

Finally, we can also make a directory to hold the `GATK`/`Picard` reports:

```
# Make reports directory
mkdir -p $REPORTS_DIRECTORY
```

### 1. Compress SAM file to BAM file

As you might suspect, because SAM files hold alignment information for all of the reads in a sequencing run and there are oftentimes millions of sequence reads, SAM files are very large and cumbersome to store. As a result, **SAM files are often stored in a binary compressed version called a BAM file**. Most software packages are agnostic to this difference and will accept both SAM and BAM files, despite BAM files not being human readable. It is generally considered best practice to keep your data in BAM format if it is being stored for long periods of time, unless you specifically need the SAM version of the alignment, in order to reduce unnecessary storage on a shared computing cluster.

**This step will happen in the same line code used for query-sorting (below).**

### 2. Query-sort alignment file

Alignment files are initally ordered by the order of the reads in the FASTQ file, which is not particularly useful. `GATK`/`Picard` can more exhaustively look for duplicates if the file is sorted by read-name (**query-sorted**). Oftentimes, when people discuss sorted BAM/SAM files, they are refering to **coordinate-sorted** BAM/SAM files. 

- **Query**-sorted BAM/SAM files are sorted based upon their read names and ordered lexiographically
- **Coordinate**-sorted BAM/SAM files are sorted by their aligned sequence name (chromosome/linkage group/scaffold) and position 

<p align="center">
<img src="../img/SAM_sorting.png" width="800">
</p>

GATK/Picard can mark and remove duplicates in either coordinate-sorted or query-sorted BAM/SAM files, however, if the alignments are query-sorted it can test secondary alignments for duplicates. A brief discussion of this nuance is discussed in the [`MarkDuplicates` manual of `GATK`/`Picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-). As a result, we will first **query**-sort our SAM file and convert it to a BAM file.

**While we query-sort the reads, we are also going to convert our SAM file to a BAM file**. We don't need to specify this conversion explicitly, because `GATK`/`Picard` will make this change by interpreting the file extensions that we provide in the `INPUT` and `OUTPUT` file options.

Add the following command to our script:

```
# Query-sort alginment file and convert to BAM
gatk SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname
```

The components of this command are:

* `gatk SortSam ` Calls `GATK`/`Picard`'s `SortSam` software package
* `--INPUT $SAM_FILE` This is where we provide the SAM input file
* `--OUTPUT $QUERY_SORTED_BAM_FILE` This is the BAM output file. 
* `--SORT_ORDER queryname` The options here are either `queryname` or `coordinate`.

> #### Why does this command look different from the Picard documentation?
> The **syntax that GATK/Picard uses** is quite particular and the syntax shown in the documentation is **not always consistent**. There are two main ways for providing input for GATK/Picard: Traditional and New (Barcalay) Syntax. Commands written in either syntax are **equally valid and produce the same output**. To better understand the different syntax, we recommend you take a look [at this short lesson](picard_syntax.md).

### 3. Mark and Remove Duplicates

An important step in processing a BAM file is to mark and remove PCR duplicates. These PCR duplicates can introduce artifacts because **regions that have preferential PCR amplification could be over-represented**. These reads are flagged by having identical mapping locations in the BAM file. Importantly, it is impossible to distinguish between PCR duplicates and identical fragments. However, one can reduce the latter by doing paired-end sequencing and providing appropriate amounts of input material. 

<p align="center">
<img src="../img/Duplicate_reads.png" width="800">
</p>

Now we will add the command to our script that allows us to mark and remove duplicates in Picard:

```
# Mark and remove duplicates
gatk MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true
```

The components of this command are:

* `gatk MarkDuplicates` Calls `GATK`/`Picard`'s `MarkDuplicates` program
* `--INPUT $QUERY_SORTED_BAM_FILE` Uses our query-sorted BAM file as input
* `--OUTPUT $REMOVE_DUPLICATES_BAM_FILE` Writes the output to a BAM file
* `--METRICS_FILE $METRICS_FILE` Creates a metrics file (required by `Picard MarkDuplicates`)
* `--REMOVE_DUPLICATES true` Not only are we going to mark/flag our duplicates, we can also remove them. By setting the `REMOVE_DUPLICATES` parameter equal to `true`, we can remove the duplicates.

### 4. Coordinate-sort the Alignment File

For most downstream processes, coordinate-sorted alignment files are required. As a result, we will need to **change our alignment file from being query-sorted to being coordinate-sorted** and we will once again use the `SortSam` command within `GATK` to accomplish this. Since this BAM file will be the final BAM file that we make and will use for downstream analyses, **we will need to create an index for it at the same time**. The command we will be using for coordinate-sorting and indexing our BAM file is:

```
# Coordinate-sort BAM file and create BAM index file
gatk SortSam \
  --INPUT $REMOVE_DUPLICATES_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true
```

The components of this command are:

* `gatk SortSam` Calls `GATK`'s `SortSam` program
* `--INPUT $REMOVE_DUPLICATES_BAM_FILE` Our BAM file once we have removed the duplicate reads. **NOTE: Even though the software is called `SortSam`, it can use BAM or SAM files as input and also BAM or SAM files as output.**
* `--OUTPUT $COORDINATE_SORTED_BAM_FILE` Our BAM output file sorted by coordinates.
* `--SORT_ORDER coordinate` Sort the output file by **coordinates**
* `--CREATE_INDEX true` Setting the `CREATE_INDEX` equal to `true` will create an index of the final BAM output. The index creation can also be accomplished by using the `BuildBamIndex` command within `Picard`, but this `CREATE_INDEX` functionality is built into many `Picard` functions, so you can often use it at the last stage of processing your BAM file to save having to run `BuildBamIndex` after.

Go ahead and save and quit. **Don't run it just yet!**

<details>
  <summary><b>Click here to see what our final <code>sbatch</code>code script for the normal sample should look like</b></summary> 
  <pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 1
#SBATCH --mem 32G
#SBATCH -o gatk_alignment_processing_normal_%j.out
#SBATCH -e gatk_alignment_processing_normal_%j.err<br>
# Load module
module load gatk/4.6.1.0<br>
# Assign file paths to variables
SAMPLE_NAME=syn3_normal
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_GRCh38.p7.sam
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/gatk/${SAMPLE_NAME}/
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`<br>
# Make reports directory
mkdir -p $REPORTS_DIRECTORY<br>
# Query-sort alginment file and convert to BAM
gatk SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname<br>
# Mark and remove duplicates
gatk MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true<br>
# Coordinate-sort BAM file and create BAM index file
gatk SortSam \
  --INPUT $REMOVE_DUPLICATES_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true<br>
</pre>
</details>


> #### Do I need to add read groups?
> Some pipelines will have you add read groups while procressing your alignment files. It is usually not necessary because you can typically do it during alignment. **If you are needing to add read groups, we recommend doing it first (before all the processing steps outlined above)**. You can use Picard `AddOrReplaceReadGroups`, which has the added benefit of allowing you to also sort your alignment file (our first step anyways) in the same step as adding the read group information. The dropdown below discusses how to add or replace read groups within `GATK`/`Picard`.
>
><details>
>  <summary><b>Click here if you need to add or replace read groups using <code>GATK</code>/<code>Picard</code></b></summary>
>    In order to add or replace read groups, we are going to use <code>GATK</code>/<code>Picard</code>'s <code>AddOrReplaceReadGroups</code> tool. First we would need to load the <code>GATK</code> module:
>  <pre>
>  # Load module
>  module load gatk/4.6.1.0</pre>
>  
>  The general syntax for <code>AddOrReplaceReadGroups</code> is:
>  <pre>
>  # Add or replace read group information
>  gatk AddOrReplaceReadGroups \
>    --INPUT $SAM_FILE \
>    --OUTPUT $BAM_FILE \
>    --RGID $READ_GROUP_ID \
>    --RGLB $READ_GROUP_LIBRARY \
>    --RGPL $READ_GROUP_PLATFORM \
>    --RGPU $READ_GROUP_PLATFORM_UNIT \
>    --RGSM $READ_GROUP_SAMPLE</pre>
>  
>  <ul><li><code>gatk AddOrReplaceReadGroups</code> This calls the <code>AddOrReplaceReadGroups</code> package within <code>GATK</code>/<code>Picard</code></li>
>    <li><code>--INPUT $SAM_FILE</code>This is your input file. It could be a BAM/SAM alignment file, but because we recommend doing this first if you need to do it, this would be a SAM file. You don't need to specifiy that it is a BAM/SAM file, <code>GATK</code>/<code>Picard</code> with figure that out from the provided extension.</li>
>    <li><code>--OUTPUT $BAM_FILE</code>This would be your output file. It could be BAM/SAM, but you would mostly likely pick BAM because you'd like to save space on the cluster. You don't need to specifiy that it is a BAM/SAM file, <code>GATK</code>/<code>Picard</code> with figure that out from the provided extension.</li>
>    <li><code>--RGID $READ_GROUP_ID</code>This is your read group ID and must be unique</li>
>    <li><code>--RGLB $READ_GROUP_LIBRARY</code>This is your read group library</li>
>    <li><code>--RGPL $READ_GROUP_PLATFORM</code>This is the platform used for the sequencing</li>
>    <li><code>--RGPU $READ_GROUP_PLATFORM_UNIT</code>This is the unit used to do the sequencing</li>
>    <li><code>--RGSM $READ_GROUP_SAMPLE</code>This is the sample name that the sequencing was done on</li></ul>
>    
>    We discussed the Read Group tags previously in the <a href="https://hbctraining.github.io/variant_analysis/lessons/03_sequence_alignment_theory.html#short-read-alignment">Sequence Alignment Theory</a> and more information on them can be found <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups">here</a>.
>  
>
>  If you would also sort your SAM/BAM file at the same time, you just need to add the <code>--SORT_ORDER</code> option to your command. If you don't add it, it will leave your reads in the same order as they were provided. The main two sort orders to be aware of are query-sorted and coordinate-sorted. A full discussion  of them can be found above. If you wanted the output to be query-sorted, then you could use:
>    
>  <pre>--SORT_ORDER queryname</pre>
>  
>  Or if you wanted them to be coordinate-sorted than you could use:
>  
>  <pre>--SORT_ORDER coordinate</pre>
>  <hr />
></details>

---

<details>
  <summary><b>Click here for alignment file processing using <code>Samtools</code></b></summary>
<br><code>Samtools</code> is another popular tool used for processing BAM/SAM files. The output from <code>Samtools</code> compared to <code>Picard</code> is largely the same. Below is the pipeline and explanation for how you would carry out the similar SAM/BAM processing steps within <code>Samtools</code>.<br>
<p align="center">
<img src="../img/Process_alignment_workflow_zoom_in_samtools.png" width="800">
</p><br>
<ol><li><details>
    <summary><b>Click here for setting up a <code>sbatch</code> script BAM/SAM Processing for the <code>Samtools</code> pipeline</b></summary>
<h2>Setting up <code>sbatch</code> Script</h2>
First, we are going to navigate to our scirpts folder and open a new <code>sbatch</code> submission script in <code>vim</code>:
<pre>
cd ~/variant_calling/scripts/
vim samtools_processing_normal.sbatch
</pre>
Next, we are going to need to set-up our <code>sbatch</code> submission script with our shebang line, description, <code>sbatch</code> directives, modules to load and file variables.
<pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Samtools<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o samtools_processing_normal_%j.out
#SBATCH -e samtools_processing_normal_%j.err<br>
# Load modules
module load gcc/14.2.0
module load samtools/1.21<br>
# Assign file paths to variables
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
REMOVED_DUPLICATES_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}removed_duplicates.bam`<br>
</pre>
<hr />
</details></li>
<li><details>
    <summary><b>Click here for <b>query</b>-sorting a SAM file and converting it to BAM for the <code>Samtools</code> pipeline</b></summary>
Similarly to <code>Picard</code>, we are going to need to initally <b>query</b>-sort our alignment. We are also going to be converting the SAM file into a BAM file at this step. Also similarly to <code>Picard</code>, we don't need to specify that our input files are BAM or SAM files. <code>Samtools</code> will use the extensions you provide it in your file names as guidance for whether you are providing it a BAM/SAM file. Below is the code we will use to <b>query</b>-sort our SAM file and convert it into a BAM file:<br>
<pre>
# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
  --threads 8 \
  -n \
  -o $QUERY_SORTED_BAM_FILE \
  $SAM_FILE
</pre>
    
The components of this line of code are:
    
<ul><li><code>samtools sort</code> This calls the sort function within <code>samtools</code>.</li>

<li><code>--threads 8</code> This tells <code>samtools</code> to use 8 threads when it multithreads this task. Since we requested 8 cores for this <code>sbatch</code> submission, let's go ahead and use them all.</li>

<li><code>-n</code> This argument tells <code>samtools sort</code> to sort by read name as opposed to the default sorting which is done by coordinate.</li>

<li><code>-O bam</code> This is declaring the output format of <code>.bam</code>.</li>

<li><code>-o $QUERY_SORTED_BAM_FILE</code> This is a <code>bash</code> variable that holds the path to the output file of the <code>samtools sort</code> command.</li>

<li><code>$SAM_FILE</code> This is a <code>bash</code> variable holding the path to the input SAM file.</li></ul>
<hr />
</details></li>

<li><details>    
<summary><b>Click here for fixing mate information for the <code>Samtools</code> pipeline</b></summary>
Next, we are going to add more mate-pair information to the alignments including the insert size and mate pair coordinates. It is important to note that for this command <code>samtools</code> relies on positional parameters for assigning the input and output BAM files. In this case the input BAM file (<code>$QUERY_SORTED_BAM_FILE</code>) needs to come before the output file (<code>$FIXMATE_BAM_FILE</code>):
    
<pre>
# Score mates
samtools fixmate \
  --threads 8 \
  -m \
  $QUERY_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE
</pre>

The parts of this command are:

<ul><li><code>samtools fixmate</code> This calls the <code>fixmate</code> command in <code>samtools</code></li>
  
<li><code>--threads 8</code> This tells <code>samtools</code> to use 8 threads when it multithreads this task.</li>

<li><code>-m</code> This will add the mate score tag that will be critically important later for <code>samtools markdup</code></li>

<li><code>$QUERY_SORTED_BAM_FILE</code> This is our input BAM file</li>

<li><code>$FIXMATE_BAM_FILE</code> This is our output BAM file</li></ul>
<hr />
</details></li>

<li><details>
<summary><b>Click here for coordinate-sorting a BAM file for the <code>Samtools</code> pipeline</b></summary>
    
Now that we have added the <code>fixmate</code> information, we need to <b>coordinate</b>-sort the BAM file. We can do that by: 

<pre>
# Sort BAM file by coordinate   
samtools sort \
  --threads 8 \
  -o $COORDINATE_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE
</pre>

We have gone through all of the these paramters already in the previous <code>samtools sort</code> command. The only difference in this command is that we are not using the <code>-n</code> option, which tells <code>samtools</code> to sort by read name. Now, we are excluding this and thus sorting by coordinates, the default setting.
<hr />
</details></li>
    
<li><details>
<summary><b>Click here for marking and removing duplicates for the <code>Samtools</code> pipeline</b></summary>

Now we are going to mark and remove the duplicate reads:
  
<pre>
# Mark and remove duplicates and then index the output file
samtools markdup \
  -r \
  --write-index \
  --threads 8 \
  $COORDINATE_SORTED_BAM_FILE \
  ${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai
</pre> 

The components of this command are:    
    
<ul><li><code>samtools markdup</code> This calls the mark duplicates software in <code>samtools</code></li>
    
<li><code>-r</code> This removes the duplicate reads</li>
    
<li><code>--write-index</code> This writes an index file of the output BAM file</li>
    
<li><code>--threads 8</code> This sets that we will be using 8 threads</li>
    
<li><code>$BAM_FILE</code> This is our input BAM file</li>
    
<li><code>${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai</code>This has two parts:
<ol><li>The first part (<code>${REMOVED_DUPLICATES_BAM_FILE}</code>) is our BAM output file with the duplicates removed from it</li>
<li>The second part (<code>##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai</code>) is a shortcut to creating a <code>.bai</code> index of the BAM file. If we use the <code>--write-index</code> option without this second part, it will create a <code>.csi</code> index file. <code>.bai</code> index files are a specific type of <code>.csi</code> files, so we need to specify it with the second part of this command to ensure that a <code>.bai</code> index file is created rather than a <code>.csi</code> index file.</li></ol></li></ul>
<hr />
</details></li>

<li><details>
<summary><b>Click here for the final normal sample <code>sbatch</code> script to do the BAM/SAM processing for the <code>Samtools</code> pipeline</b></summary>

The final script should look like:
    
<pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Samtools<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o samtools_processing_normal_%j.out
#SBATCH -e samtools_processing_normal_%j.err<br>
# Load modules
module load gcc/14.2.0
module load samtools/1.21<br>
# Assign file paths to variables
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_normal_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
REMOVED_DUPLICATES_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}removed_duplicates.bam`<br>
# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
  --threads 8 \
  -n \
  -o $QUERY_SORTED_BAM_FILE \
  $SAM_FILE<br>
# Score mates
samtools fixmate \
  --threads 8 \
  -m \
  $QUERY_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE<br>
# Sort BAM file by coordinate   
samtools sort \
  --threads 8 \
  -o $COORDINATE_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE<br>
# Mark and remove duplicates and then index the output file
samtools markdup \
  -r \
  --write-index \
  --threads 8 \
  $COORDINATE_SORTED_BAM_FILE \
  ${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai<br>
</pre>
<hr />
</details></li>
<li><details>
<summary><b>Click here for the final tumor sample <code>sbatch</code> script to do the BAM/SAM processing for the <code>samtools</code> pipeline</b></summary>
In order to create the tumor <code>sbatch</code> submission script to process the BAM/SAM file using <code>samtools</code>, we will once again use <code>sed</code>:<br>
<pre>
sed &#39;s/normal/tumor/g&#39; samtools_processing_normal.sbatch &gt; samtools_processing_tumor.sbatch  
</pre>
The final <code>sbatch</code> submission script for the tumor sample should look like:
<pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use in GATK using Samtools<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 8
#SBATCH --mem 16G
#SBATCH -o samtools_processing_tumor_%j.out
#SBATCH -e samtools_processing_tumor_%j.err<br>
# Load modules
module load gcc/14.2.0
module load samtools/1.15.1<br>
# Assign file paths to variables
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/syn3_tumor_GRCh38.p7.sam
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
FIXMATE_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}fixmates.bam`
COORDINATE_SORTED_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}coordinate_sorted.bam`
REMOVED_DUPLICATES_BAM_FILE=`echo ${QUERY_SORTED_BAM_FILE%query_sorted.bam}removed_duplicates.bam`<br>
# Sort SAM file and convert it to a query name sorted BAM file
samtools sort \
  --threads 8 \
  -n \
  -o $QUERY_SORTED_BAM_FILE \
  $SAM_FILE<br>
# Score mates
samtools fixmate \
  --threads 8 \
  -m \
  $QUERY_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE<br>
# Sort BAM file by coordinate   
samtools sort \
  --threads 8 \
  -o $COORDINATE_SORTED_BAM_FILE \
  $FIXMATE_BAM_FILE<br>
# Mark and remove duplicates and then index the output file
samtools markdup \
  -r \
  --write-index \
  --threads 8 \
  $COORDINATE_SORTED_BAM_FILE \
  ${REMOVED_DUPLICATES_BAM_FILE}##idx##${REMOVED_DUPLICATES_BAM_FILE}.bai<br>
</pre>
</details></li></ol>
<hr />
</details>

---

## Exercises

**1.** When inspecting a SAM file you see the following order:

<p align="center">
<img src="../img/Sort_order_question.png" width="500">
</p>

Is this SAM file's sort order likely: unsorted, query-sorted, coordinate-sorted or is it ambiguous?

---

## Creating the Tumor SAM/BAM processing script
    
Similar to the `bwa` script, we will now need use `sed` to create a `sbatch` script that will be used for processing the tumor SAM file into a BAM file that can be used as input to GATK. The `sed` command to do this would be:
  
```
sed 's/normal/tumor/g' gatk_alignment_processing_normal.sbatch > gatk_alignment_processing_tumor.sbatch  
```

_As a result your tumor `GATK`/`Picard` alignment processing script should look almost identical but `normal` has been replaced by `tumor`._


<details>
  <summary><b>Click here to see what our final <code>sbatch</code>code script for the tumor sample should look like </b></summary> 
  <pre>
#!/bin/bash
# This sbatch script is for processing the alignment output from bwa and preparing it for use<br> 
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-04:00:00
#SBATCH -c 1
#SBATCH --mem 32G
#SBATCH -o gatk_alignment_processing_tumor_%j.out
#SBATCH -e gatk_alignment_processing_tumor_%j.err<br>
# Load module
module load gatk/4.6.1.0<br>
# Assign file paths to variables
SAMPLE_NAME=syn3_tumor
SAM_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/alignments/${SAMPLE_NAME}_GRCh38.p7.sam
REPORTS_DIRECTORY=/home/${USER}/variant_calling/reports/picard/${SAMPLE_NAME}/
QUERY_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}query_sorted.bam`
REMOVE_DUPLICATES_BAM_FILE=`echo ${SAM_FILE%sam}remove_duplicates.bam`
METRICS_FILE=${REPORTS_DIRECTORY}/${SAMPLE_NAME}.remove_duplicates_metrics.txt
COORDINATE_SORTED_BAM_FILE=`echo ${SAM_FILE%sam}coordinate_sorted.bam`<br>
# Make reports directory
mkdir -p $REPORTS_DIRECTORY<br>
# Query-sort alginment file and convert to BAM
gatk SortSam \
  --INPUT $SAM_FILE \
  --OUTPUT $QUERY_SORTED_BAM_FILE \
  --SORT_ORDER queryname<br>
# Mark and remove duplicates
gatk MarkDuplicates \
  --INPUT $QUERY_SORTED_BAM_FILE \
  --OUTPUT $REMOVE_DUPLICATES_BAM_FILE \
  --METRICS_FILE $METRICS_FILE \
  --REMOVE_DUPLICATES true<br>
# Coordinate-sort BAM file and create BAM index file
gatk SortSam \
  --INPUT $REMOVE_DUPLICATES_BAM_FILE \
  --OUTPUT $COORDINATE_SORTED_BAM_FILE \
  --SORT_ORDER coordinate \
  --CREATE_INDEX true
</pre>
</details>

## Submitting scripts for processing
  
Now we are ready to submit our normal and tumor `GATK`/`Picard` processing scripts to the O2 cluster. However, we might have a problem. If you managed to go quickly into this lesson from the previous lesson, **your `bwa` alignment scripts may still be running and your SAM files are not complete yet!**
  
First, we need to check the status of our `bwa` scripts and we can do this with the command:
  
```
squeue --me
```
  
* **If you have `bwa` jobs still running,** then wait for them to complete (less than 2 hours) before continuing. There are ways to queue jobs together in SLURM using the `--dependency` option in `sbatch`. This is outside the scope of this workshop, but we cover this in the [automation lesson](automation_of_variant_calling.md), if you are interested in learning more. For now, just hang tight until the alignment is complete.
  
* **If the only job running is your interactive job,** then it should be time to start your `Picard` processing scripts. You can go ahead and submit your `sbatch` scripts for `Picard` processing with:
  
```
sbatch gatk_alignment_processing_normal.sbatch
sbatch gatk_alignment_processing_tumor.sbatch
```

> **NOTE:** These scripts will take ~2 hours for each to run!
  
[Next Lesson >>](05_alignment_QC.md)

[Back to Schedule](../schedule/README.md)
  
***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
