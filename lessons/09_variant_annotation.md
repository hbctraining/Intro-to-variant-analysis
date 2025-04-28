---
title: "Variant Annotation"
author: "Will Gammerdinger, Meeta Mistry"
date: "May 27, 2024"
---

Approximate time: 45 minutes


## Learning Objectives

- Annotate a VCF file for functional impacts with `SnpEff`
- Differentiate between an unannotated and annotated VCF file

# Variant Annotation

Now that we have a filtered VCF file for our samples, we want to annotate those variants to figure out what impacts they could have on our samples. Information we might be interested in knowing:

- Which genes are the variants we see are located in
- Where in the gene model do these variants reside (5' UTR, coding sequence, etc.)
- How do these variants impact the gene model (creating a premature stop codon, silent mutation, etc.)

In order to do these types of analyses, we have to merge our variants with annotated transcript information and we will use `SnpEff` to do this. `SnpEff` is extremely fast and also is bundled with `SnpSift`, which we will use downstream to prioritize our variants. 

<p align="center">
<img src="../img/Annotate_variants_workflow.png" width="400">
</p>

## SnpEff

[SnpEff](http://pcingola.github.io/SnpEff/) uses transcriptome annotations to predict the functional impacts of mutations in a VCF file and will modify the INFO field of the VCF file to carry the predicted functional impacts.

### Searching for SnpEff Databases

The first step in annotating your VCF file is finding the appropriate SnpEff database to use so that the annotations are consistent with your version of the reference genome. These databases hold the gene model information that is critical for annotating variants. There are two paths for getting a SnpEff database for your analysis:

- `SnpEff` has tens of thousands of genome databases pre-built that you can use and they represent many publicly availible genomes. **We will be using a pre-built database.**
- In the rare case that your genome is not represented, you can build one to use by following the instructions [here](http://pcingola.github.io/SnpEff/se_build_db/). Unless you are working with a rare organism or an unpublished reference genome, it is highly unlikely that you will need to build the database yourself. 

To see if your genome of interest is in the `SnpEff` database, we first need to load the `SnpEff` module:

```
module load snpEff/5.2f
```

With the `SnpEff` module loaded, let's use the following command to browse all of the currently available genomes:

```
java -jar $SNPEFF/snpEff.jar databases | less
```

The first column is the database name and the second column in the `Genus_species` for the organism. There is also a database download link where the database can be downloaded at but this can be ignored as SnpEff will automatically download the database if needed. As you can see there are tens of thousands of these pre-built databases. So let's exit the `less` buffer page and **see which hg38 databases are available**:

```
java -jar $SNPEFF/snpEff.jar databases | grep "hg38" 
```

We can see that this build of SnpEff has three possible hg38 databases that we can use for annotation. We will use the one labelled hg38. Now that we have found the database that we would like to use for our analysis, we can run `SnpEff`.

### Running SnpEff

Before we get into `SnpEff` we need to discuss cancer-mode in `SnpEff`. In order to run `SnpEff` in cancer-mode, we need to do two things:

**1)** Understand what cancer-mode is

**2)** Format our VCF file accordingly to accommodate for cancer-mode

#### Understanding Cancer-mode

There can be loci in the genome where, due to variation in the population, an individual's germline alleles differs from the reference sequence. Furthermore, if a somatic mutation happens at one of these loci, then from SnpEff's perspective the mutational event will look like a change from the reference allele to the somatic mutant rather than a change from the normal sample allele to the the somatic mutant. SnpEff's cancer-mode is designed to address this. The illustration below shows an example of this case. 

<p align="center">
<img src="../img/Cancer_mode.png" width="600">
</p>

A more detailed explaination of `SnpEff`'s cancer-mode can be found on their [website](https://pcingola.github.io/SnpEff/snpeff/cansersamples/).

#### Setting up Cancer-mode

In order to run cancer-mode, we will need to append to our VCF file an additional header line that contains information that `SnpEff` can use when determining, which sample is the normal sample and which sample is the tumor sample. To append the header lines we will be using a package within [`bcftools`](https://samtools.github.io/bcftools/bcftools.html). First, move to your `scripts` directory and create a new file named `syn3_normal_syn3_tumor_pedigree_header.txt` using `vim`:

```
cd ~/variant_calling/scripts/
vim syn3_normal_syn3_tumor_pedigree_header.txt
```

Once inside this file, we just need to add this text to the file:

```
##PEDIGREE=<Derived=syn3_tumor,Original=syn3_normal>
```

This text will act as the additional header information that we need to add. Then save and exit the file.

Next, we can start writing our `sbatch` submission script named `variant_annotation_normal_tumor.sbatch` using `vim`:

```
vim variant_annotation_normal_tumor.sbatch
```

Once inside insert mode, can can enter the shebang line, description and `SBATCH` directives:

```
#!/bin/bash
# Using SnpEff to annotate our variants

# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_annotation_syn3_normal_syn3_tumor_%j.out
#SBATCH -e variant_annotation_syn3_normal_syn3_tumor_%j.err
```

Next, we will add the line to load the modules that we will need: 

```
# Load modules
module load gcc/9.2.0
module load bcftools/1.14
module load snpEff/5.2f
```

Also, we will add our variables:

```
# Assign variables
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/snpeff/
SAMPLE_NAME=mutect2_syn3_normal_syn3_tumor
REFERENCE_SEQUENCE_NAME=GRCh38.p7
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=GRCh38.p7.RefSeq
DATADIR=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
FILTERED_VCF_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-pass-filt-LCR.vcf
PEDIGREE_HEADER_FILE=/home/$USER/variant_calling/scripts/syn3_normal_syn3_tumor_pedigree_header.txt
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.dbSNP.vcf.gz
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf
```

Next, we need to create a directory to hold our reports:

```
# Create reports directory
mkdir -p $REPORTS_DIRECTORY
```

Now, we can add our `bcftools` command to append the header line required for cancer-mode in `SnpEff` to our VCF file:

```
# Append Header
bcftools annotate \
  --header-lines $PEDIGREE_HEADER_FILE \
  $FILTERED_VCF_FILE \
  > $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER
```

Here is what each part of that command does:

  - `bcftools annotate` This calls the `annotate` function within `bcftools`
  - `--header-lines $PEDIGREE_HEADER_FILE` This is the header information that we want to append to the end of the VCF file's header lines
  - `$FILTERED_VCF_FILE` This is our input VCF file without the appended header
  - `> $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER` This is the output file with the appended header

Next, we need to add our `SnpEff` command:

```
# Run SnpEff
java -jar -Xmx4g $SNPEFF/snpEff.jar  eff \
  -dataDir $DATADIR \
  -cancer \
  -noLog \
  -csvStats $CSV_STATS \
  -s $HTML_REPORT \
  $REFERENCE_DATABASE \
  $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER \
  > $SNPEFF_ANNOTATED_VCF_FILE
```

Let's breakdown this command and discuss each argument:

- `java -jar $SNPEFF/snpEff.jar  eff` `Snpeff` This calls the `eff` package within `SnpEff` to annotate variants.
- `-dataDir $DATADIR` This is the path to the data directory that holds the `SnpEff` annotations
- `-cancer` Performs cancer-mode 
- `-noLog` This does not report command usage to <code>SnpEff</code>'s server
- `-csvStats $CSV_STATS` This produces a flat-text file with summary statistics regarding the variants annotated. (Optional)
- `-s $HTML_REPORT` This creates an HTML file with summary statistics regarding the variants annotated. This HTML file is mostly just an HTML stylized version of the CSV file above. (Optional)
- `$REFERENCE_DATABASE` This is the `SnpEff` database we are going to use for the annotation.
- `$FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER` This is the input VCF file to be annotated
- `> $SNPEFF_ANNOTATED_VCF_FILE` The output of `SnpEff` will be redirected into this VCF file.

### Annotate SNPs with dbSNP

dbSNP is a public-domain archive maintained by NCBI for simple genetic polymorphisms in humans. It is designed as a repository for researchers to be able to submit novel variants to as well as query for variants submitted by others. The variants within dbSNP include:

- SNPs
- Small-scale insertions/deletions
- Microsatellite repeat variations
- Retrotransposable element insertions

<p align="center">
<img src="../img/dbSNP_logo.png" width="800">
</p>

In addition to adding the annotations that `SnpEff` provides regarding types of mutations, we can also add annotations from dbSNP regarding our variants. These annotations will add lots of information to the `INFO` field, but it will also populate the `ID` field in the VCF file with the dbSNP ID for each variant, if it exists. In order to do these annotations, we need to have access to a dbSNP VCF file with the annotations along with an index of the VCF file (a `.tbi` file). **Note, that we have already created these files for you.** 

<details>
  <summary><b>Click here to see how to obtain the dbSNP VCF file and how to index it</b></summary>
  The first step is that we need to download out dbSNP VCF file from the NCBI FTP site. We can do this with a <code>curl</code> command like: 

<pre>
# YOU DO NOT NEED TO DO THIS
curl -o GRCh38.p7.dbSNP.vcf.gz -L https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
</pre>

We used the <code>curl</code> command earlier when demonstrating how you could download the LCR BED file. But as a refresher, we can breakdown this command:
    
<ul><li><code>curl</code> This calls the <code>curl</code> option</li>
  <li><code>-o GRCh38.p7.dbSNP.vcf.gz</code> This tells <code>curl</code> what we want the downloaded file to be named</li>
  <li><code>-L https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz</code> This is the linked file we want to download</li></ul>
    
Next, we are going to index this file. While it is not necesscary for us to index this VCF file when using <code>SnpSift</code>, it does speed up the annotation. Importantly, <code>bcftools</code> <b><i>DOES REQUIRE</i></b> this index file and it also will require files that are indexed to also be compressed in <code>.gz</code> format. As a result we are going to leave, our file <code>GRCh38.p7.dbSNP.vcf.gz</code> compressed because:
    <ol><li>There is no need to uncompress it</li>
      <li>It will take up more space on the cluster if we do uncompress it</li>
      <li>Both <code>SnpSift</code> and <code>bcftools</code> can work from the compressed version</li></ol>

We are going to be using <code>tabix</code>, which is part of the <code>HTSlib</code> module. First, we will need to load the <code>HTSlib</code> module, which also requires us to load the <code>gcc</code> module as well:
    
<pre>
module load gcc/9.2.0
module load htslib/1.14
</pre>
    
In order to index our dbSNP file using <code>tabix</code>, we just need to run the following command:
    
<pre>
# YOU DO NOT NEED TO DO THIS
tabix GRCh38.p7.dbSNP.vcf.gz
</pre>

We can break down this command as:
    <ul><li><code>tabix</code> This is the command we use to call <code>tabix</code></li>
     <li><code>GRCh38.p7.dbSNP.vcf.gz</code> This is the VCF file we want to index.</li></ul>
    
After this finishes, we can see that we now have a <code>GRCh38.p7.dbSNP.vcf.gz.tbi</code> file in addition to our <code>GRCh38.p7.dbSNP.vcf.gz</code> file.
<hr />
</details>



We will be implementing these dbSNP annotations by using `SnpSift`, which is part of the `SnpEff` package.

Let's go ahead and add the `SnpSift` line that we are going to use to annotate our VCF file:

```
# Use dbSNP VCF to annotate our VCF
java -jar $SNPEFF/SnpSift.jar annotate \
  $DBSNP_DATABASE \
  -tabix \
  -noLog \
  $SNPEFF_ANNOTATED_VCF_FILE \
  > $DBSNP_ANNOTATED_VCF_FILE
```

Let's discuss each part of this command:

- `java -jar $SNPEFF/SnpSift.jar annotate` This calls the `annotate` function within the `SnpSift` package
- `$DBSNP_DATABASE` This is the dbSNP database that we want to use for annotations
  - We had set this variable earlier in the lesson to point to the database we had created for you `DBSNP_DATABASE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.dbSNP.vcf.gz`
- `-tabix` This is letting `SnpSift` know that the dbSNP database has a VCF index that it can use to speed up the annotation
- `-noLog` This does not report command usage to <code>SnpEff</code>'s server
- `$SNPEFF_ANNOTATED_VCF_FILE` This is the input VCF file we want to be annotated
- `> $DBSNP_ANNOTATED_VCF_FILE` This is the dbSNP annotated output VCF file

<details>
  <summary><b>Click here to see how to annotate our VCF file with the dbSNP annotation in <code>bcftools</code></b></summary>
The first thing we are going to need to do is load the modules that we will be using:<br>
<pre>
module load gcc/9.2.0
module load bcftools/1.14
</pre>
Assuming we have already indexed our dbSNP VCF file, the first thing that we are going to need to do is compress the VCF file that we wish to annotate with:
<pre>
gzip -c $SNPEFF_ANNOTATED_VCF_FILE &gt; ${SNPEFF_ANNOTATED_VCF_FILE}.gz
</pre>
Let's talk about this command:<br>
  <ul><li><code>gzip</code> This calls the function we are going to use for the compression</li>
    <li><code>-c</code> This is an option that keeps the original file. Otherwise, when the compression happens, the uncompressed file is removed</li>
    <li><code>$SNPEFF_ANNOTATED_VCF_FILE</code> This is the file that we want to compress</li>
    <li><code>&gt; ${SNPEFF_ANNOTATED_VCF_FILE}.gz</code> This is the compressed output file</li></ul><br>
Now we are ready to annotate the VCF file with the dbSNP annotations. The command to do so is:
<pre>
bcftools annotate \
  -a $DBSNP_DATABASE \
  -c ID,INFO \
  -O v \
  -o $DBSNP_ANNOTATED_VCF_FILE \
  ${SNPEFF_ANNOTATED_VCF_FILE}.gz
</pre>
Let's explain each part of this command:
  <ul><li><code>bcftools annotate</code> This calls the <code>annotate</code> function within <code>bcftools</code></li>
    <li><code>-a $DBSNP_DATABASE</code> This provides the dbSNP VCF to use for the annotations</li>
    <li><code>-c ID,INFO</code> These are the fileds we would like to add information to from our dbSNP file</li>
    <li><code>-O v</code> This is the output format. The options for this are:
      <table>
        <tr>
          <th>Abbreviation</th>
          <th>Output Format</th>
        </tr>
        <tr>
          <td>b</td>
          <td>Compressed BCF file</td>
        </tr>
        <tr>
          <td>u</td>
          <td>Uncompressed BCF file</td>
        </tr>
        <tr>
          <td>z</td>
          <td>Compressed VCF file</td>
        </tr>
        <tr>
          <td>v</td>
          <td>Uncompressed VCF file</td>
        </tr>
      </table></li>
    <li><code>-o $DBSNP_ANNOTATED_VCF_FILE</code> This is the dbSNP annotated output VCF file</li>
    <li><code>${SNPEFF_ANNOTATED_VCF_FILE}.gz</code> This is the input VCF file</li></ul>
<hr />
</details>


<details>
  <summary><b>Click here to see what our final <code>sbatch</code>code script for annotating our variants should look like</b></summary> 
  <pre>
#!/bin/bash
# Using SnpEff to annotate our variants<br>
# Assign sbatch directives
#SBATCH -p priority
#SBATCH -t 0-02:00:00
#SBATCH -c 1
#SBATCH --mem 8G
#SBATCH -o variant_annotation_syn3_normal_syn3_tumor_%j.out
#SBATCH -e variant_annotation_syn3_normal_syn3_tumor_%j.err<br>
# Load modules
module load gcc/9.2.0
module load bcftools/1.14
module load snpEff/5.2f<br>
# Assign variables
REPORTS_DIRECTORY=/home/$USER/variant_calling/reports/snpeff/
SAMPLE_NAME=mutect2_syn3_normal_syn3_tumor
REFERENCE_SEQUENCE_NAME=GRCh38.p7
CSV_STATS=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.csv"`
HTML_REPORT=`echo -e "${REPORTS_DIRECTORY}annotation_${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-effects-stats.html"`
REFERENCE_DATABASE=GRCh38.p7.RefSeq
DATADIR=/n/groups/hbctraining/variant_calling/reference/snpeff/data/
FILTERED_VCF_FILE=/n/scratch/users/${USER:0:1}/${USER}/variant_calling/vcf_files/${SAMPLE_NAME}_${REFERENCE_SEQUENCE_NAME}-pass-filt-LCR.vcf
PEDIGREE_HEADER_FILE=/home/$USER/variant_calling/scripts/syn3_normal_syn3_tumor_pedigree_header.txt
FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER=${FILTERED_VCF_FILE%.vcf}.pedigree_header.vcf
SNPEFF_ANNOTATED_VCF_FILE=${FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER%.vcf}.snpeff.vcf
DBSNP_DATABASE=/n/groups/hbctraining/variant_calling/reference/GRCh38.p7.dbSNP.vcf.gz
DBSNP_ANNOTATED_VCF_FILE=${SNPEFF_ANNOTATED_VCF_FILE%.vcf}.dbSNP.vcf<br>
# Create reports directory
mkdir -p $REPORTS_DIRECTORY<br>
# Append Header
bcftools annotate \
  --header-lines $PEDIGREE_HEADER_FILE \
  $FILTERED_VCF_FILE \
  > $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER<br>
# Run SnpEff
java -jar -Xmx4g $SNPEFF/snpEff.jar  eff \
  -dataDir $DATADIR \
  -cancer \
  -noLog \
  -csvStats $CSV_STATS \
  -s $HTML_REPORT \
  $REFERENCE_DATABASE \
  $FILTERED_VCF_FILE_WITH_PEDIGREE_HEADER \
  > $SNPEFF_ANNOTATED_VCF_FILE<br>
# Use dbSNP VCF to annotate our VCF
java -jar $SNPEFF/SnpSift.jar annotate \
  $DBSNP_DATABASE \
  -tabix \
  -noLog \
  $SNPEFF_ANNOTATED_VCF_FILE \
  > $DBSNP_ANNOTATED_VCF_FILE
</pre>
<hr />
</details>

Submit this script using:

```
sbatch variant_annotation_normal_tumor.sbatch
```

### Output

Let's take a look at our output from `SnpEff` now to see how our VCF file has been modified to include variant annotations. Let's open up our annotated VCF file in `less`:

```
less /n/scratch/users/${USER:0:1}/$USER/variant_calling/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.snpeff.dbSNP.vcf 
```

Scroll down to right before your variants are and you will notice that `SnpEff` has inserted two lines into your VCF file. First, it has inserted:

`##SnpEffVersion="4.3g (build 2016-11-28 08:31), by Pablo Cingolani"`

This is just a bit of extra metadata on your VCF file letting you know the version of `SnpEff` that was used to generate this VCF file. While it is best practice to have your code for annotating a VCF file archived in a script  like we have done, you could encounter a situation where you've either lost this script or you are analyzing someone else's work and need to know which version of `SnpEff` they used and you can find it here. 

```
##SnpEffCmd="SnpEff  -cancer -csvStats /home/$USER/variant_calling/reports/snpeff/annotation_mutect2_syn3_normal_syn3_tumor_GRCh38.p7-effects-stats.csv -s /home/$USER/variant_calling/reports/snpeff/annotation_mutect2_syn3_normal_syn3_tumor_GRCh38.p7-effects-stats.html GRCh38.p7.RefSeq /n/scratch/users/${USER:0:1}/$USER/variant_calling/vcf_files/mutect2_syn3_normal_syn3_tumor_GRCh38.p7-pass-filt-LCR.pedigree_header.vcf "
```

This tells you the command that was used by SnpEff to create this SnpEff annotated file. Thus, between the aforementioned version line and this line, you should be able to reproduce the annotation of any SnpEff annotation as long as the annotation database is part of the built-in databases.

Let's go ahead and compare the first entry of the filtered VCF file with the annotated VCF file:

#### Filtered VCF

```
1       1324300 .       G       A       .       PASS    AS_FilterStatus=SITE;AS_SB_TABLE=47,6|11,0;ClippingRankSum=0.390;DP=68;ECNT=1;FS=2.373;GERMQ=93;MBQ=27,27;MFRL=337,338;MMQ=60,60;MPOS=27;MQ=60.00;MQ0=0;MQRankSum=0.000;NALOD=1.54;NLOD=9.88;POPAF=6.00;ReadPosRankSum=-0.125;TLOD=24.65        GT:AD:AF:DP:F1R2:F2R1:SB        0/0:33,0:0.028:33:12,0:20,0:32,1,0,0    0/1:20,11:0.367:31:6,6:12,5:15,5,11,0
```

#### Annotated VCF

```
1       1324300 .       G       A       .       PASS    AS_FilterStatus=SITE;AS_SB_TABLE=47,6|11,0;ClippingRankSum=0.39;DP=68;ECNT=1;FS=2.373;GERMQ=93;MBQ=27,27;MFRL=337,338;MMQ=60,60;MPOS=27;MQ=60;MQ0=0;MQRankSum=0;NALOD=1.54;NLOD=9.88;POPAF=6;ReadPosRankSum=-0.125;TLOD=24.65;ANN=A|upstream_gene_variant|MODIFIER|CPTP|CPTP|transcript|NM_001029885.1|protein_coding||c.-2611G>A|||||463|,A|upstream_gene_variant|MODIFIER|CPTP|CPTP|transcript|XM_005244802.1|protein_coding||c.-3008G>A|||||456|,A|upstream_gene_variant|MODIFIER|CPTP|CPTP|transcript|XM_005244801.3|protein_coding||c.-2611G>A|||||843|,A|upstream_gene_variant|MODIFIER|CPTP|CPTP|transcript|XM_011542200.2|protein_coding||c.-2611G>A|||||1315|,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|XM_011541647.1|protein_coding|1/18|c.28+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|NM_001256456.1|protein_coding|1/18|c.-428+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|NM_001256460.1|protein_coding|1/17|c.-167+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|NM_001256462.1|protein_coding|1/14|c.28+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|NM_001256463.1|protein_coding|1/14|c.28+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|NM_017871.5|protein_coding|1/16|c.28+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|XM_017001558.1|protein_coding|1/18|c.-438+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|XM_017001557.1|protein_coding|1/17|c.-361+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|XM_011541648.1|protein_coding|1/18|c.-91+281C>T||||||,A|intron_variant|MODIFIER|CPSF3L|CPSF3L|transcript|XM_011541650.1|protein_coding|1/16|c.-254+281C>T||||||  GT:AD:AF:DP:F1R2:F2R1:SB        0/0:33,0:0.028:33:12,0:20,0:32,1,0,0    0/1:20,11:0.367:31:6,6:12,5:15,5,11,0
```

You can see that the only difference between these files is that we have appended onto the `INFO` field. This appendage contains lots of information such as the type of modification it is, the gene symbol of the gene modified and the accession number of the transcript modified. Because gene models overlap, you may have a single variant that alters multiple transcript models. In this case, each alteration will be separated by a `,`. We can also see that this SNP was not part of the dbSNP database that we had and thus it didn't have its `ID` field annotated with its dbSNP ID.

## Exercise

SnpEff generates an HTML report that gives a nice overview of the annotated variants. Use FileZilla to download the file  called `annotation_mutect2_syn3_normal_syn3_tumor_GRCh38.p7-effects-stats.html` from your reports directory. If you have trouble, the report can also be downloaded from [here](https://hbctraining.github.io/variant_analysis/html_reports/annotation_mutect2_syn3_normal_syn3_tumor_GRCh38.p7-effects-stats.html). Use this HTML report to answer the following questions:

**1)** How many variants are in the VCF file?

**2)** Which chromosome had the most variants?

**3)** How many variants were missense variants?

***

Now that we have successfully annotated our variants, let's talk about how we will use SnpSift to be able to prioritize our variants in the next lesson!


[Next lesson >>](10_variant_prioritization.md)

[Back to Schedule](../schedule/README.md)

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
