# WGBS Pipeline

The WGBS analysis pipeline processes raw FASTQ data through a series of steps including adapter trimming, quality control, genome mapping and DMR. Using Singularity for reproducibility and supporting batch analysis of multiple samples.

## Workflow Diagram

![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/WGBS_pipeline.png)

## Requirements

1. **Recommended System Configuration**:

     * 8-core CPU
     * 80 GB RAM

2. **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

     * **Step 1: Install System Dependencies**

       ```bash
       # Update package lists and install dependencies
       sudo apt-get update
       sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
            libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
       curl wget git
       ```
     
     * **Step 2: Install Go Language**
     
       ```bash
       # Download and install Go
       wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
       sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
       rm go1.21.3.linux-amd64.tar.gz
     
       # Configure Go environment variables and apply them
       echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
       echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
       source ~/.bashrc
       ```
     
     * **Step 3: Download, Build, and Install Singularity**
     
       ```bash
       # Note: The script navigates to /mnt/share/software. 
       # You can change this to your preferred directory for source code.
       cd /mnt/share/software
     
       # Download the Singularity CE source code
       wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz
     
       # Extract the archive and clean up
       tar -xvzf singularity-ce-4.0.1.tar.gz
       rm singularity-ce-4.0.1.tar.gz
       cd singularity-ce-4.0.1
     
       # Configure the build
       ./mconfig
     
       # Build Singularity (this can be time-consuming)
       cd builddir
       make
     
       # Install Singularity to the system
       sudo make install
       ```
     
     * **Step 4: Verify the Installation**
     
       ```bash
       # Check the installed version
       singularity --version
       
       # Display help information
       singularity -h
       ```
3. **snakemake**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. 

```bash
pip install snakemake
```

4. **Pipeline Files**:

     * `WGBS.smk`
     * `wgbs.sif` (The Singularity container)

5. **Reference Data**: A directory containing all necessary reference files (e.g., bitmapperBS indices and GTF annotation, etc.).

**Note on Sequencing Type:**
This pipeline supports paired-end (PE) data. The example below shows the format for paired-end data.

### 1. Prepare the Reference Genome

The pipeline requires several pre-built reference files. Below are the steps to generate them for the Solanum lycopersicum Micro-Tom genome using NCBI annotations.

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p data/GCA_040143515.1/
cd data/GCA_040143515.1/
```

#### Common Reference Files

We require the following base files:

**Download Reference Files:**

![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/ncbi_genome.png)

##### Build bitmapperBS Indices:
```bash
singularity exec wgbs.sif bitmapperBS --index data/GCA_040143515.1/GCA_040143515.1_ASM4014351v1_genomic.fna
```

### 2. Prepare the test fastq data

The data run by this pipeline is from SRR34756619,SRR34756620,SRR34756625 and SRR34756626 in the SRA database.The specific processing method is as follows

#### Create a dedicated directory for the test sra data:

```bash
mkdir -p data/samples
cd data/samples
```
#### Download the test sra data
```bash
prefetch SRR34756619
prefetch SRR34756620
prefetch SRR34756625
prefetch SRR34756626
```
#### Convert sra data to fastq data
```bash
fastq-dump --split-files SRR34756619\SRR34756619.sra
fastq-dump --split-files SRR34756620\SRR34756620.sra
fastq-dump --split-files SRR34756625\SRR34756625.sra
fastq-dump --split-files SRR34756626\SRR34756626.sra
```

#### Randomly sample fastq data:

```bash
singularity exec WGBS.sif seqtk sample -s100 SRR34756619_1.fastq 27780 > SRR34756619_R1.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756619_2.fastq 27780 > SRR34756619_R2.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756620_1.fastq 27780 > SRR34756620_R1.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756620_2.fastq 27780 > SRR34756620_R2.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756625_1.fastq 27780 > SRR34756625_R1.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756625_2.fastq 27780 > SRR34756625_R2.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756626_1.fastq 27780 > SRR34756626_R1.fastq
singularity exec WGBS.sif seqtk sample -s100 SRR34756626_2.fastq 27780 > SRR34756626_R2.fastq
pigz -p 8 SRR34756619_R1.fastq
pigz -p 8 SRR34756619_R2.fastq
pigz -p 8 SRR34756620_R1.fastq
pigz -p 8 SRR34756620_R2.fastq
pigz -p 8 SRR34756625_R1.fastq
pigz -p 8 SRR34756625_R2.fastq
pigz -p 8 SRR34756626_R1.fastq
pigz -p 8 SRR34756626_R2.fastq
```
### 3. Check the WGBS snakemake workflow

The specific snakemake workflow is as follows

#### Dry-run:

Here /mnt/liuq/WGBS/ represents the root directory.

```bash
snakemake -np \
          -s WGBS.smk \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/WGBS:/root/"
```
#### Draw a detailed DAG picture

Here /mnt/liuq/WGBS represents the root directory.

```bash
snakemake -s WGBS.smk \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/WGBS:/root/" \
          --dag  | \
dot -Tsvg > dag.svg
```

![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/dag.png)
### 4. Check the current working directory

The initial working structure for snakemake in /mnt/liuq/WGBS/:
## Input Structure and Interpretation

Before the pipeline, the input directory contain several files and directories. Below is a detailed explanation of what each file you should supply and how it can be used.
```bash
├── wgbs.sif
├── WGBS.smk
├── config.yaml
├── dag.svg
├── data
│   ├── GCA_040143515.1
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.fai
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.bwt
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.occ
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.sa
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.pac
│   │   └── GCA_040143515.1_ASM4014351v1_genomic.fna.index.methy
│   └── samples
│       ├── SRR34756619_R1.fastq.gz
│       ├── SRR34756619_R2.fastq.gz
│       ├── SRR34756620_R1.fastq.gz
│       ├── SRR34756620_R2.fastq.gz
│       ├── SRR34756625_R1.fastq.gz
│       ├── SRR34756625_R2.fastq.gz
│       ├── SRR34756626_R1.fastq.gz
│       └── SRR34756626_R2.fastq.gz
├── samplesheet.txt
└── scripts
    └── methylSig_WGBS.r
```
## Running

After adding the corresponding parameters in config.yaml,you can execute the pipeline using a single command, the only important parameter is thread (**-j**).
### Example Commands

Here /mnt/liuq/WGBS/ represents the root directory.

```bash
snakemake -s WGBS.smk \
		  -j 30 \
          --use-singularity \
          --singularity-args "--bind /mnt/liuq/WGBS/:/root/"
```

Then delete the intermediate files and folders

```bash
rm -rf qc mapped data/multiqc_report/*/multiqc_data
rm text/mBias*.txt text/*.log text/*.svg data/*/*.zip
```
## Output Structure and Interpretation

After the pipeline completes, the output directory will contain several files and directories. Below is a detailed explanation of what each file is and how it can be used.
```bash
├── wgbs.sif
├── WGBS.smk
├── config.yaml
├── dag.svg
├── data
│   ├── GCA_040143515.1
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.fai
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.bwt
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.occ
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.index.sa
│   │   ├── GCA_040143515.1_ASM4014351v1_genomic.fna.index.bs.pac
│   │   └── GCA_040143515.1_ASM4014351v1_genomic.fna.index.methy
│   ├── multiqc_report
│   │   ├── SRR34756619
│   │   │   └── multiqc_report.html
│   │   ├── SRR34756620
│   │   │   └── multiqc_report.html
│   │   ├── SRR34756625
│   │   │   └── multiqc_report.html
│   │   └── SRR34756626
│   │       └── multiqc_report.html
│   ├── samples
│   │   ├── SRR34756619_R1.fastq.gz
│   │   ├── SRR34756619_R2.fastq.gz
│   │   ├── SRR34756620_R1.fastq.gz
│   │   ├── SRR34756620_R2.fastq.gz
│   │   ├── SRR34756625_R1.fastq.gz
│   │   ├── SRR34756625_R2.fastq.gz
│   │   ├── SRR34756626_R1.fastq.gz
│   │   └── SRR34756626_R2.fastq.gz
│   ├── SRR34756619
│   │   ├── SRR34756619_R1_val_1_fastqc.html
│   │   └── SRR34756619_R2_val_2_fastqc.html
│   ├── SRR34756620
│   │   ├── SRR34756620_R1_val_1_fastqc.html
│   │   └── SRR34756620_R2_val_2_fastqc.html
│   ├── SRR34756625
│   │   ├── SRR34756625_R1_val_1_fastqc.html
│   │   └── SRR34756625_R2_val_2_fastqc.html
│   └── SRR34756626
│       ├── SRR34756626_R1_val_1_fastqc.html
│       └── SRR34756626_R2_val_2_fastqc.html

├── samplesheet.txt
├── scripts
│   └── methylSig_WGBS.R
└── text
    └── diff_simple_gr.txt
```

- **`multiqc_report.html`**
  Open multiqc_report.html in a web browser to explore all sections interactively.
  
   - **Application**: This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .
  
       - **General Statistics**: A combined table summarizing important metrics for each sample:
       
       ![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/general_statistic.png)
       
       - **FastQC**: Quality-control metrics on raw and trimmed reads, including  
        'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
           'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
           'Sequence Length Distribution', 'Sequence Duplication Levels',  
           'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.
         
       - **Sequence Quality Histograms**: The mean quality value across each base position in the read. 
       
       
       ![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/fastqc_per_base_sequence_quality_plot.png)
       
      - **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.  
      
       ![](https://github.com/Haolab-BIG/WGBS-Processing-Pipeline/raw/main/picture/fastqc_per_sequence_quality_scores_plot.png)
      
  
- **`fastqc.html(zip)`**
  
  - **Content**: Open fastqc.html in a web browser to explore all sections interactively.Similar to the above multiqc results.
  
  - **Application**:This is the first file you should check to assess the overall quality of your sequencing data. It helps identify problematic samples (e.g., high duplication) .Similar to the above multiqc results.


- **`diff_simple_gr.txt`**
  - **Content**: MethylSig differential m5C analysis results.
  - **Application**:  It used to examine differences m5C sites  between the case and the control samples.

```bash
seqnames        start   end     width   strand  meth_case       meth_control    meth_diff       direction       stat    pvalue  fdr
CM079145.1      33865278        33865278        1       *       0       0       0       case    -0.023764693866117      0.981040302304091       1
CM079145.1      33865293        33865293        1       *       0       0       0       case    -0.023764693866117      0.981040302304091       1
CM079145.1      33865338        33865338        1       *       0       0       0       case    -0.023764693866117      0.981040302304091       1
CM079145.1      33900747        33900747        1       *       0       0       0       case    0.0634846256951526      0.949380601465285       1
CM079145.1      33907232        33907232        1       *       0       0       0       case    0.00183197533327639     0.998538295983419       1
```

## Video Tutorials

Watch this video tutorial to see a complete walk through of running the pipeline:

https://github.com/Haolab-BIG/WGBS/raw/main/snakemake.mp4
