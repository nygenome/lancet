lancet
======

Lancet is a somatic variant caller (SNVs and indels) for short read data. Lancet uses a localized micro-assembly strategy to detect somatic mutation with high sensitivity and accuracy on a tumor/normal pair.
Lancet is based on the colored de Bruijn graph assembly paradigm where tumor and normal reads are jointly analyzed within the same graph. On-the-fly repeat composition analysis and self-tuning k-mer strategy are used together to increase specificity in regions characterized by low complexity sequences. Lancet requires the raw reads to be aligned with BWA (See [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) description for more info). Lancet is implemented in C++.

* Author: Giuseppe Narzisi, [New York Genome Center](https://www.nygenome.org)

Lancet is freely available for academic and non-commercial research purposes ([`LICENSE.txt`](https://github.com/nygenome/lancet/blob/master/LICENSE.txt)).  

### References

Narzisi G, Corvelo A, Arora K, Bergmann E, Shah M, Musunuri R, Emde AK, Robine N, Vacic V, Zody MC. *Genome-wide somatic variant calling using localized colored de Bruijn graphs.* 
<b>Communications Biology</b>, Nature Research publishing, volume 1, Article number: 20, 2018 (DOI:<a href=https://www.nature.com/articles/s42003-018-0023-9>10.1038/s42003-018-0023-9</a>). Also available at <b>CSHL bioRxiv</b> 196311; 2017 (DOI: <a href="https://doi.org/10.1101/196311">10.1101/196311</a>)

Rajeeva Musunuri, Kanika Arora, André Corvelo, Minita Shah, Jennifer Shelton, Michael C. Zody, Giuseppe Narzisi. *Somatic variant analysis of linked-reads sequencing data with Lancet.* <b>CSHL bioRxiv</b> 2020.07.04.158063; doi: https://doi.org/10.1101/2020.07.04.158063


### Downloading and building lancet

Building and running lancet from source requires a GNU-like environment with 

1. GCC (version >= 4.8.x)
2. GNU Make
3. GNU CMake (version >= 3.0)

Lancet can be built on most Linux installations. Most distributions already ship with all the c++ libraries that Lancet depends on: `lzma bz2 z dl pthread curl crypto deflate`. Compilation on MacOS require Xcode and Xcode command line tools installed. Lancet source code is available through github and can be obtained and compiled with the following command:

```sh
git clone git://github.com/nygenome/lancet.git
cd lancet
make
```
### Basic usage

A simple lancet command should look something like this:

```
lancet --tumor T.bam --normal N.bam --ref ref.fa --reg 22:1-51304566 --num-threads 8 > out.vcf
```

The command above detects somatic variants in a tumor/normal pair of bam files (*T.bam* and *N.bam*) for chromosome 22 using 8 threads and saves the variant calls in the output VCF file *out.vcf*. 

**NOTE**: a genomic region must be always specified via the --reg option with format *"chr:start-end"*. Single chromosome names are also supported (e.g., --reg 22).

### Genome-wide scan

Due to its pure local-assembly strategy, Lancet currently has longer runtimes compared to standard alignment-based variant callers. For whole-genome sequencing studies it is highly recommended to split the analysis by chromosome and then merge the results. Splitting the work by chromosome will also reduce the overall memory requirements to analyze the whole-genome data.

```
NUMBER_OF_AUTOSOMES=22
for chrom in `seq 1 $NUMBER_OF_AUTOSOMES` X Y; do
	qsub \
	-N lancet_chr${chrom} \
	-cwd \
	-pe smp 8 \
	-q dev.q \
	-j y \
	-b y \
	"lancet --tumor T.bam --normal N.bam --ref ref.fa --reg $chrom --num-threads 8 > ${chrom}.vcf"
done

// merge VCF files
```
The previous command shows an exemplary submission of multiple parallel lancet jobs, one for each human chromosome, to the Sun Grid Engine queuing system.

### Linked-Reads analysis

The recommended command line options for [10x Genomics](https://www.10xgenomics.com/) linked-reads analysis are:

```
lancet --linked-reads --primary-alignment-only --tumor T.bam --normal N.bam --ref ref.fa --reg chr1 --num-threads 8 > out.vcf
```
where:

* **linked-reads** activates the linked-reads mode to support LongRanger BAM format (e.g., BX and HP tags, etc.).
* **primary-alignment-only** forces the program to only use the primary alignment of each read for the analysis.

[LongRanger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) BAMs are directy supported, however, for improved accuarcy, we highly recommend to process the BAMs with the [MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) program from [Picard Tools](https://broadinstitute.github.io/picard/), which marks PCR duplicates more accurately than LongRanger.

### Output

Lancet generates in output the list of variants in VCF format (v4.1). All variants (SNVs and indels either shared, specific to the tumor, or specific to the normal) are exported in output. Following VCF conventions, high quality variants are flagged as **PASS** in the FILTER column. For non-PASS variants the FILTER info reports the list of filters that are not satisfied by each variant.

### Filters

The list of filters applied and the thresholds used for filtering are included in the VCF header section. For example:

```
##FILTER=<ID=LowFisherScore,Description="low Fisher's exact test score for tumor-normal allele counts (<5)">
```
The previous filter means that a variant flagged as **LowFisherScore** has not met the minimum Fisher's exact test score threshold for tumor-normal allele counts (default 5).

Below is the current list of filters:

1. **LowCovNormal**: low coverage in the normal
2. **HighCovNormal**: high coverage in the normal
3. **LowCovTumor**: low coverage in the tumor
4. **HighCovTumor**: high coverage in the tumor
5. **LowVafTumor**: low variant allele frequency in the tumor
6. **HighVafNormal**: high variant allele frequency in the normal
7. **LowAltCntTumor**: low alternative allele count in the tumor
8. **HighAltCntNormal**: high alternative allele count in the normal
9. **LowFisherScore**: low Fisher's exact test score for tumor-normal allele counts
10. **LowFisherSTR**: low Fisher's exact test score for tumor-normal STR allele counts
11. **StrandBias**: rejects variants where the vast majority of alternate alleles are seen in a single direction
12. **STR**: microsatellite mutation

### Visual inspection of the DeBruijn graph

The DeBruijn graph representation of a genomic region can be exported to file in [DOT](http://www.graphviz.org/doc/info/lang.html) format using the -A flag. 

**NOTE:** *The following procedure does not scale to large graphs. Please render a graph only to inspect a small genomic region of a few hundred base pairs. The -A flag must not be used during regular variant calling over large genomic regions.*

For example the following command:

```
lancet -A --tumor T.bam --normal N.bam --ref ref.fa --reg chr:start-end > out.vcf
```

will export the DeBruijn graph after every stage of the assembly (low covergae removal, tips removal, compression) to the following set of files:

1. chr:start-end.0.dot (initial graph)
2. chr:start-end.1l.cX.dot (after first low coverage nodes removal)
3. chr:start-end.2c.cX.dot (after compression)
4. chr:start-end.3l.cX.dot (after second low coverage nodes removal)
5. chr:start-end.4t.cX.dot (after tips removal)
6. chr:start-end.final.cX.dot (final graph)

Where X is the number of the correspending connected component (in most cases only one). 
These files can be rendered using the utilities available in the [Graphviz](http://www.graphviz.org/) visualization software package. Specifically we reccomand using the **sfdp** utlity which draws undirected graphs using the ``spring'' model and it uses a multi-scale approach to produce layouts of large graphs in a reasonably short time.

```
sfdp -Tpdf file.dot -O
```

For large graphs, Adobe Acrobat Reader may have troubles rendering the graph, in that case we recommend opening the PDF file using the "Preview" image viewer software available in MacOS.

An exemplary graph (before removal of low coverage nodes and tips) for a short region containing a somatic variant would look like this one:

![initial graph](https://github.com/nygenome/lancet/blob/master/doc/img/graph.insertion.png)

where the blue nodes are k-mers shared by both tumor and normal; the white nodes are k-mer with low support (e.g., sequencing errors); the red nodes are k-mers only present in the tumor node.

A clean bubble whitin a graph is displayed below:

![initial graph](https://github.com/nygenome/lancet/blob/master/doc/img/clean_bubble.png)

The final graph (after compression) containing one single variant is depicted below. Yellow and orange nodes are the source and sink nodes respectively 

<img src="https://github.com/nygenome/lancet/blob/master/doc/img/final_graph.png" width="400">

### Complete command-line options

```
  |                           |
  |      _` | __ \   __|  _ \ __|
  |     (   | |   | (     __/ |
 _____|\__,_|_|  _|\___|\___|\__|

Program: lancet (micro-assembly somatic variant caller)
Version: 1.1.0, October 18 2019
Contact: Giuseppe Narzisi <gnarzisi@nygenome.org>

Usage: lancet [options] --tumor <BAM file> --normal <BAM file> --ref <FASTA file> --reg <chr:start-end>
 [-h for full list of commands]

Required
   --tumor, -t              <BAM file>    : BAM file of mapped reads for tumor
   --normal, -n             <BAM file>    : BAM file of mapped reads for normal
   --ref, -r                <FASTA file>  : FASTA file of reference genome
   --reg, -p                <string>      : genomic region (in chr:start-end format)
   --bed, -B                <string>      : genomic regions from file (BED format)

Optional
   --min-k, k                <int>         : min kmersize [default: 11]
   --max-k, -K               <int>         : max kmersize [default: 101]
   --trim-lowqual, -q        <int>         : trim bases below qv at 5' and 3' [default: 10]
   --min-base-qual, -C       <int>         : minimum base quality required to consider a base for SNV calling [default: 17]
   --quality-range, -Q       <char>        : quality value range [default: !]
   --min-map-qual, -b        <int>         : minimum read mapping quality in Phred-scale [default: 15]
   --max-as-xs-diff, -Z      <int>         : maximum difference between AS and XS alignments scores [default: 5]
   --tip-len, -l             <int>         : max tip length [default: 11]
   --cov-thr, -c             <int>         : min coverage threshold used to select reference anchors from the De Bruijn graph [default: 5]
   --cov-ratio, -x           <float>       : minimum coverage ratio used to remove nodes from the De Bruijn graph [default: 0.01]
   --low-cov, -d             <int>         : low coverage threshold used to remove nodes from the De Bruijn graph [default: 1]
   --max-avg-cov, -u         <int>         : maximum average coverage allowed per region [default: 10000]
   --window-size, -w         <int>         : window size of the region to assemble (in base-pairs) [default: 600]
   --padding, -P             <int>         : left/right padding (in base-pairs) applied to the input genomic regions [default: 250]
   --dfs-limit, -F           <int>         : limit dfs/bfs graph traversal search space [default: 1000000]
   --max-indel-len, -T       <int>         : limit on size of detectable indel [default: 500]
   --max-mismatch, -M        <int>         : max number of mismatches for near-perfect repeats [default: 2]
   --num-threads, -X         <int>         : number of parallel threads [default: 1]
   --node-str-len, -L        <int>         : length of sequence to display at graph node (default: 100)

Filters
   --min-alt-count-tumor, -a  <int>        : minimum alternative count in the tumor [default: 3]
   --max-alt-count-normal, -m <int>        : maximum alternative count in the normal [default: 0]
   --min-vaf-tumor, -e        <float>      : minimum variant allele frequency (AlleleCov/TotCov) in the tumor [default: 0.04]
   --max-vaf-normal, -i       <float>      : maximum variant allele frequency (AlleleCov/TotCov) in the normal [default: 0]
   --min-coverage-tumor, -o   <int>        : minimum coverage in the tumor [default: 4]
   --max-coverage-tumor, -y   <int>        : maximum coverage in the tumor [default: 1000000]
   --min-coverage-normal, -z  <int>        : minimum coverage in the normal [default: 10]
   --max-coverage-normal, -j  <int>        : maximum coverage in the normal [default: 1000000]
   --min-phred-fisher, -s     <float>      : minimum fisher exact test score [default: 5]
   --min-phred-fisher-str, -E <float>      : minimum fisher exact test score for STR mutations [default: 25]
   --min-strand-bias, -f      <float>      : minimum strand bias threshold [default: 1]

Short Tandem Repeat parameters
   --max-unit-length, -U      <int>        : maximum unit length of the motif [default: 4]
   --min-report-unit, -N      <int>        : minimum number of units to report [default: 3]
   --min-report-len, -Y       <int>        : minimum length of tandem in base pairs [default: 7]
   --dist-from-str, -D        <int>        : distance (in bp) of variant from STR locus [default: 1]

Flags
   --linked-reads, -J            : linked-reads analysis mode
   --primary-alignment-only, -I  : only use primary alignments for variant calling
   --XA-tag-filter, -O           : skip reads with multiple hits listed in the XA tag (BWA only)
   --active-region-off, -W       : turn off active region module
   --kmer-recovery, -R           : turn on k-mer recovery (experimental)
   --print-graph, -A             : print graph (in .dot format) after every stage
   --verbose, -v                 : be verbose
   --more-verbose, -V            : be more verbose
```

### Funding

Informatics Technology for Cancer Research ([ITCR](https://itcr.cancer.gov/)) under the NCI R21 award [1R21CA220411-01A1](https://projectreporter.nih.gov/project_info_description.cfm?aid=9507249&icde=40260779&ddparam=&ddvalue=&ddsub=&cr=1&csb=default&cs=ASC&pball=).