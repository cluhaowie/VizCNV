Phased_BAF
================

<img src="BAF_phasing.svg" width="717" height="384" />

``` bash
git clone git@github.com:cluhaowie/VizCNV.git
cd VizCNV
chmod +x helper/Phased_BAF.R
```

## Usage and option summary

Usage:

``` bash
helper/Phased_BAF.R -I path/to/joint_genotyped.vcf -C auto
```

or:

``` bash
helper/Phased_BAF.R -I path/to/joint_genotyped.vcf -C chr1
```

| Options      | Description                                                                                                                                                                                                                                                                              |
|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| -I, –-input  | Specifies the path to a VCF file for SNP phasing and categorization. The VCF file should contain trio sample genotypes, and a joint genotyped VCF is recommended for optimal results. For joint genotyping, ref to \*\*Example Usage:\*\* \`-I path/to/joint_genotyped.vcf\`             |
| -O, –-output | Defines the path for the output file to store summarized information of each category. If no file path is provided, the information will be printed to standard output (stdout). The default behavior is to write to stdout. \*\*Example Usage:\*\* \`-O path/to/output.tsv\` (Optional) |
| -C, –-chr    | Specifies the chromosome(s) for analysis. You can specify a single chromosome with \`-C chr1\`, all chromosomes with \`-C all\`, or all autosomes with \`-C auto\`. \*\*Example Usage:\*\* \`-C chr1\`, \`-C all\`, \`-C auto\`                                                          |
| -R, –-ref    | Sets the reference genome to be used. The default reference genome is hg38. \*\*Example Usage:\*\* \`-R hg19\` (Optional; Default is \`-R hg38\`)                                                                                                                                        |
| -N, –-cores  | Determines the number of processor cores to be used for parallel processing. The default setting is 1 core. \*\*Example Usage:\*\* \`-N 4\` (Optional; Default is \`-N 1\`)                                                                                                              |

for example:

``` bash
$ helper/Phased_BAF.R
Usage: helper/Phased_BAF.R [options]
Options:
    -I CHARACTER, --input=CHARACTER
        Path to the joint genotyped VCF file

    -O CHARACTER, --output=CHARACTER
        Path to the output file

    -C CHARACTER, --chr=CHARACTER
        chr region for analysis, e.g. chr1, all, auto

    -R CHARACTER, --ref=CHARACTER
        Reference genome, e.g. hg38, hg19, [default -R hg38]

    -N CHARACTER, --cores=CHARACTER
        Number of cores for parallel analysis, [default -N 1]

    -h, --help
        Show this help message and exit
$ helper/Phased_BAF.R -I HG00405/fam1.filtered.snps.vcf.gz -C auto > fam1.phased_snp_summary.tsv
$ head -n 5 fam1.phased_snp_summary.tsv
B_InhFrom   count   total   freq    chrom
AOH_signal  214697  397205  0.540519379161894   chr1
ME  3767    397205  0.00948376782769603 chr1
Notphased   36127   397205  0.0909530343273624  chr1
fam1-2_HG00404  73751   397205  0.18567490338742    chr1
fam1-3_HG00403  68863   397205  0.173368915295628   chr1
```

## Visualization of Biallelic SNV Phasing and Categorization Metrics for Trio Samples: HG00405 (Child), HG00404 (Mother), HG00403 (Father) on Autosomes

<figure>
<img src="Phased_BAF_files/figure-gfm/unnamed-chunk-1-1.png"
alt="HG00405_metrics" />
<figcaption aria-hidden="true">HG00405_metrics</figcaption>
</figure>

## Visualization of Biallelic SNV Phasing and Categorization Metrics for Trio Samples with genomewide parental UPD: BG1383 (Child), BG1384 (Mother), BG1385 (Father) on Autosomes

<figure>
<img src="Phased_BAF_files/figure-gfm/unnamed-chunk-2-1.png"
alt="UPD_case" />
<figcaption aria-hidden="true">UPD_case</figcaption>
</figure>

<img src="UPD_case.png" width="724" />
