# 🎉 ImmunoTyper2 🧬

**ImmunoTyper2** is a powerful tool for Immunoglobulin Variable Gene genotyping and CNV analysis from whole genome sequencing (WGS) short reads using ILP Optimization. Check out our [paper here](https://www.cell.com/cell-systems/fulltext/S2405-4712(22)00352-0?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471222003520%3Fshowall%3Dtrue) for more details.

📢 **New Feature:** 
* **Prefix Consistency Analysis (Gurobi Solver Only)**: Evaluate the stability and confidence of allele calls using prefix consistency metrics, outputting detailed accuracy insights. (See the new output file description below for more details).
- New optional ILP solver! The free-to-use [Google OR Tools](https://developers.google.com/optimization) with the CP SAT solver can be used with the `--solver or-tools` param ! 🎉 Note it is best suited to lower coverage or less complex loci (e.g. does not work on 30x IGHV).
- Now supporting V gene calling for **all IG and TR loci!** (IGH, IGL, IGK, TRA, TRB, TRG, TRD) 
- **Novel Variant Calling**: Identify novel variants for all gene calls using FreeBayes and WhatsHap. Variants are listed in `<prefix>-<gene_type>-novel-variants.txt` and phased VCFs are available in `<prefix>-<gene_type>-novel_variant_vcfs/<gene_id>_variants.vcf`.

## 🚀 Installation

### Singularity / Docker

For the easiest installation, **we recommend using Singularity** by pulling the Docker image available on DockerHub at `cdslsahinalp/immunotyper-sr`.

To run the image with Singularity (commonly used on HPCs), use the following command:

```sh
singularity pull docker://cdslsahinalp/immunotyper-sr
singularity run -B <BAM_DIRECTORY>:<BAM_DIRECTORY> -B <OUTPUT_PATH>:/output immunotyper-sr_latest.sif <OPTIONAL ARGUMENTS> <BAM_DIRECTORY>/<BAM_FILE> 
```

You can also run the image with Docker, however **this method has not been tested**:

```sh
docker pull cdslsahinalp/immunotyper-sr
docker run -v <BAM_DIRECTORY>:<BAM_DIRECTORY> -v <OUTPUT_PATH>:/output immunotyper-sr <OPTIONAL ARGUMENTS> <BAM_DIRECTORY>/<BAM_FILE> 
```

### Conda + Pip

If you already have BWA installed and prefer not to create a new environment, you can download the latest release binary (see right toolbar) and install it with pip:

```
pip install <binary.whl>
```

For the best experience, we recommend setting up a clean environment first:

```
conda create -n immunotyper-SR -c bioconda python=3.8 bwa bowtie2 freebayes whatshap
conda install -y -n immunotyper-SR -c gurobi gurobi
conda install -y -n immunotyper-SR -c conda-forge samtools
conda activate immunotyper-SR
pip install <binary.whl>
```

### Environment and Dependencies

Installing ImmunoTyper-SR with pip will automatically install these dependencies:

- [biopython](https://biopython.org/)
- [dill](https://pypi.org/project/dill/)
- [gurobipy](https://www.gurobi.com/documentation/9.5/quickstart_mac/cs_grbpy_the_gurobi_python.html)
- [logbook](https://logbook.readthedocs.io/en/stable/)
- [pysam](https://pysam.readthedocs.io/en/latest/api.html)
- [ortools](https://developers.google.com/optimization)
- [wurlitzer](https://pypi.org/project/wurlitzer/)

In addition to the above, you will need 

1.  [BWA mem](http://bio-bwa.sourceforge.net/bwa.shtml) mapper. We recommend using a new conda environment for the installation, which you can also use to install BWA:

```
conda create -n immunotyper-SR -c bioconda python=3.8 bwa samtools
conda activate immunotyper-SR
pip install <binary.whl>
```

2.  [Gurobi](https://www.gurobi.com/) solver configured with a valid license. Note OR-Tools is also available as an alternative ILP solver, and is configured as the default solver for TCR loci. You can get a free academic license for Gurobi [here](https://www.gurobi.com/academia/academic-program-and-licenses/). 

To check that gurobi is correctly configured, run `gurobi_cl` from a shell.

### Installing from source

If the binary fails to install, you can build the tool from source:


```
conda create -n immunotyper-SR -c bioconda python=3.8 bwa bowtie2 freebayes whatshap
conda install -y -n immunotyper-SR -c gurobi gurobi
conda install -y -n immunotyper-SR -c conda-forge samtools
conda activate immunotyper-SR
git clone git@github.com:algo-cancer/ImmunoTyper-SR.git ./ImmunoTyper-SR
cd ImmunoTyper-SR
python -m pip install --upgrade  build
python -m build
pip install dist/<.tar.gz or .whl build>
```

## 🛠️ Running ImmunoTyper-SR:

### 🔬 Prefix Consistency Confidence Metric

* **What is Prefix Consistency?**
    * Prefix consistency is a score (typically 0-5) indicating how consistently an allele appears with the same copy number across increasingly relaxed Integer Linear Programming (ILP) solutions compared to the primary optimal solution. Higher values generally indicate greater confidence in the allele call.

* **Source of Accuracy Metrics & Optimal Threshold Calculation**:
    * The accuracy metrics ("Probability of TP", "Optimal Threshold Gene Accuracy") that inform this analysis and are presented in the associated output file (see below) are derived from a pre-computed data resource. These reference metrics were established by analyzing 40 samples from the 1000 Genomes Project, using Human Pangenome Reference Consortium (HPRC) assemblies as ground truth.
    * The `optimal_threshold` for each gene, used in this report, was determined by maximizing the F-beta score (with $\beta=0.5$, thereby prioritizing precision over recall) on this reference dataset. This methodology aims to provide high-confidence allele identifications. Further details on this validation and the prefix consistency metric are covered in our upcoming paper (currently under review).


After installing with pip, use the command immunotyper-SR. The only required input is a BAM file. Outputs are generated in the current working directory, where <prefix> is the input BAM filename without the extension:

- <prefix>-<gene_type>_functional_allele_calls.txt: List of functional alleles called.
- <prefix>-<gene_type>_allele_calls.txt: Includes pseudogenes.
- <prefix>-<gene_type>-immunotyper-debug.log: Log file.
- <prefix>-<gene_type>-novel-variants.tsv: Novel variants called using FreeBayes and WhatsHap. Format: `gene position ref alt`.
- <prefix>-<gene_type>-novel_variant_vcfs/<gene_id>_variants.vcf: Phased VCFs variants. NOTE these include all variants called relative to wildtype, so include any variants from non-wildtype called alleles as well as any present novel alleles.
- <prefix>-<gene_type>-read_assignment: Contains FASTA and BAM files of reads assigned to each called allele.
- <prefix>-<GENE_TYPE>-solutions-with-consistency-confidence.tsv: Provides a detailed analysis of allele call confidence based on prefix consistency for alleles found in the optimal solution. File Columns: `Allele`, `Functional status`, `Consistency value`, `Probability of TP`, `Passes Optimal Threshold`, `Optimal Threshold Gene Accuracy`. **Note: This output is currently only generated when using the Gurobi solver.** For a detailed explanation of Prefix Consistency and how the accuracy metrics are derived, please see the "🔬 Prefix Consistency Analysis Details" section.


IMPORTANT: If your BAM was mapped to GRCh37 use the `--hg37` flag. 

```
$ immunotyper-SR --help
usage: immunotyper-SR [-h] [--gene_type {ighv,iglv,trav,igkv,trbv,trdv,trgv}] [--output_dir OUTPUT_DIR] [--ref REF] [--hg37] [--solver {gurobi,or-tools}] [--bwa BWA] [--max_copy MAX_COPY]
                      [--landmarks_per_group LANDMARKS_PER_GROUP] [--landmark_groups LANDMARK_GROUPS] [--stdev_coeff STDEV_COEFF] [--seq_error_rate SEQ_ERROR_RATE] [--solver_time_limit SOLVER_TIME_LIMIT]
                      [--debug_log_path DEBUG_LOG_PATH] [--write_cache_path WRITE_CACHE_PATH] [--threads THREADS] [--no_coverage_estimation] [--save_extracted_reads] [--solution_precision SOLUTION_PRECISION]
                      [--no_vcf] [--no_read_assignment] [--multi_band_solutions]
                      bam_path

ImmunoTyper-SR: Ig Genotyping using Short Read WGS

positional arguments:
  bam_path              Input BAM file

optional arguments:
  -h, --help            show this help message and exit
  --gene_type {ighv,iglv,trav,igkv,trbv,trdv,trgv}
                        Specify which genes to target
  --output_dir OUTPUT_DIR
                        Path to output directory. Outputs txt file of allele calls with prefix matching input BAM file name.
  --ref REF             Path to the reference FASTA to decode CRAM files. Option is not used if bam_path is not a CRAM.
  --hg37                Flag if BAM mapped to GRCh37 not GRCh38
  --solver {gurobi,or-tools}
                        Choose ilp solver
  --bwa BWA             path to bwa executible if not in $PATH
  --max_copy MAX_COPY   Maximum number of allele copies to call
  --landmarks_per_group LANDMARKS_PER_GROUP
                        Number of landmarks per group to use (default = 6)
  --landmark_groups LANDMARK_GROUPS
                        Number of landmark groups to use (default = 6)
  --stdev_coeff STDEV_COEFF
                        Standard deviation scaling coefficient (default = 1.5)
  --seq_error_rate SEQ_ERROR_RATE
                        Expected sequence error rate (default = 0.02)
  --solver_time_limit SOLVER_TIME_LIMIT
                        Time limit for ILP solver in hours
  --debug_log_path DEBUG_LOG_PATH
                        Path to write log
  --write_cache_path WRITE_CACHE_PATH
                        Specific location and name of allele db sam mapping cache
  --threads THREADS     Max number of threads to use
  --no_coverage_estimation
                        Disables empirical coverage
  --save_extracted_reads
                        Save the extracted reads FASTA file instead of deleting it after use
  --solution_precision SOLUTION_PRECISION
                        Optimality gap parameter for the solver, interger value: 1e-{solution_precision}
  --no_vcf              Do not write VCF files for novel variants
  --no_read_assignment  Do not write read assignment files
  --multi_band_solutions
                        Calculate and write multi band solutions
```