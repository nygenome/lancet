lancet
======

Lancet is a somatic variant caller (SNVs and indels) for short read data. Lancet uses a localized micro-assembly strategy to detect somatic mutation with high sensitivity and accuracy on a tumor/normal pair.
Lancet is based on the colored de Bruijn graph assembly paradigm where tumor and normal reads are jointly analyzed within the same graph. On-the-fly repeat composition analysis and self-tuning k-mer strategy are used together to increase specificity in regions characterized by complex repeat structures. Lancet requires the raw reads to be aligned with BWA (See [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) description for more info). Lancet is implemented in C++.

* Version: 1.0.0
* Author: Giuseppe Narzisi, [New York Genome Center](https://www.nygenome.org)

Lancet is freely available for academic and non-commercial research purposes ([`LICENSE.txt`](https://github.com/nygenome/lancet/blob/master/LICENSE.txt)).  

### Downloading and building lancet

Building and running lancet from source requires a GNU-like environment with 

1. GCC
2. GNU Make
3. GNU CMake

It should be possible to build lancet on most Linux installations 
or on a Mac installation with [Xcode and Xcode command line tools] installed.
Lancet is available through github and can be obtained and compiled with following command:

```sh
git clone git://github.com/nygenome/lancet.git
cd lancet
make
```
### Basic usage

A simple lancet command should look something like this:

```
Lancet --tumor T.bam --normal N.bam --ref ref.fa --reg 22-1-51304566 --num-threads 8
```

The command above will detect somatic variants in the tumor/normal pair of bam files (T.bam and N.bam) for chromosome 22 using 8 threads.

### Output

Lancet generates in output the list of variants in VCF format (v4.1). All variants (SNVs and indels either shared, specific to the tumor, or specific to the normal) are exported in output. Following VCF conventions, high quality variants are flaged as **PASS** in the FILTER column. For non-PASS variants the FILTER info reports the list of filters that were not satified by the variant.

### Filters

The list of filters applied and the thresholds used for filtering are included in the header section. For example:

```
##FILTER=<ID=LowFisherScore,Description="low Fisher's exact test score for tumor-normal allele counts (<10)">
```
The previous filter means that a variant flagged as **LowFisherScore** has not met the minimum Fisher's exact test score threshold for tumor-normal allele counts (default 10).

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
10. **StrandBias**: rejects variants where the vast majority of alternate alleles are seen in a single direction

### Complete command-line options

```
Usage: Lancet [options] --tumor <BAM file> --normal <BAM file> --ref <FASTA file> --reg <chr:start-end>
 [-h for full list of commands]

Required
   --tumor, -t              <BAM file>    : BAM file of mapped reads for tumor
   --normal, -n             <BAM file>    : BAM file of mapped reads for normal
   --ref, -r                <FASTA file>  : FASTA file of reference genome
   --reg, -p                <string>      : genomic region (in chr:start-end format)

Optional
   --min-k, k                <int>         : min kmersize [default: 11]
   --max-k, -K               <int>         : max kmersize [default: 100]
   --trim-lowqual, -q        <int>         : trim bases below qv at 5' and 3' [default: 10]
   --quality-range, -Q       <char>        : quality value range [default: !]
   --min-map-qual, -b        <inr>         : minimum read mapping quality in Phred-scale [default: 0]
   --tip-len, -l             <int>         : max tip length [default: 11]
   --cov-thr, -c             <int>         : coverage threshold [default: 5]
   --cov-ratio, -x           <float>       : minimum coverage ratio [default: 0.01]
   --max-avg-cov, -u         <int>         : maximum average coverage allowed per region [default: 10000]
   --low-cov, -d             <int>         : low coverage threshold [default: 1]
   --window-size, -w         <int>         : window size of the region to assemble (in base-pairs) [default: 600]
   --dfs-limit, -F           <int>         : limit dfs search space [default: 1000000]
   --path-limit, -P          <int>         : limit on number of paths to report [default: 0]
   --max-indel-len, -T       <int>         : limit on size of detectable indel [default: 250]
   --max-mismatch, -M        <int>         : max number of mismatches for near-perfect repeats [default: 3]
   --num-threads, -X         <int>         : number of parallel threads [default: 1]
   --rg-file, -g             <string>      : read group file

Filters
   --min-alt-count-tumor, -a  <int>        : minimum alternative count in the tumor [default: 4]
   --max-alt-count-normal, -m <int>        : maximum alternative count in the normal [default: 0]
   --min-vaf-tumor, -e        <float>      : minimum variant allele frequency (AlleleCov/TotCov) in the tumor [default: 0.05]
   --max-vaf-normal, -i       <float>      : maximum variant allele frequency (AlleleCov/TotCov) in the normal [default: 0]
   --min-coverage-tumor, -o   <int>        : minimum coverage in the tumor [default: 4]
   --max-coverage-tumor, -y   <int>        : maximum coverage in the tumor [default: 1000000]
   --min-coverage-normal, -z  <int>        : minimum coverage in the normal [default: 10]
   --max-coverage-normal, -j  <int>        : maximum coverage in the normal [default: 1000000]
   --min-phred-fisher, -s     <float>      : minimum fisher exact test score [default: 10]

Flags
   -A            : print graph (in .dot format) after every stage
   -L <len>      : length of sequence to display at graph node (default: 100)
   -v            : be verbose
   -V            : be more verbose
```