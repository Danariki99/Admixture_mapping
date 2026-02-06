# ADMIXTURE MAPPING PROJECT

## abstract
n recent years, large-scale genome-wide association studies (GWAS) and sequencing projects have identified over 70,000 associations between genetic variants and human traits or diseases. However, the majority of these studies rely on individuals who are genetically similar to the European superpopulation reference panels, reflecting the composition of biobanks enriched for individuals with predominantly European-like genetic backgrounds. This bias limits insights into the full spectrum of human genetic diversity and its impact on disease.  Admixture mapping offers a complementary approach to GWAS by leveraging differences in haplotype frequencies across populations represented in reference panels to identify risk loci for complex traits. Here, we perform large-scale admixture mapping across the UK Biobank and All of Us Research Program datasets, examining associations between local haplotype ancestry and 172 phenotypes. Our approach not only replicates previously reported associations but also identifies novel loci with important implications for human disease. Additionally, we apply a fine-mapping strategy to the 31 identified genomic regions by analyzing subsets of individuals stratified by their predominant ancestry components. This allows the same genomic region to yield distinct association signals across ancestry-defined subsamples. Overall, this approach localized 55 association signals and highlighted the importance of modeling ancestry structure to refine genetic associations and quantify their contribution to trait variability.

## Authors 
Riccardo Smeriglio<sup>2</sup> Sonia Moreno-Grau<sup>1</sup>, Daniel Mas-Montserrat<sup>1,3</sup> Christophe Thomassin<sup>1,7</sup>, Guhan Ventakaraman<sup>1</sup>, Caterina Fuses<sup>4,6</sup>, Manuel A. Rivas<sup>1</sup>, Alessandro Savino<sup>2</sup>, Jordi Abante<sup>4,6</sup>, Stefano Di Carlo<sup>2</sup> , Alexander G. Ioannidis<sup>1,3,7</sup>.
Affiliations
1. Department of Biomedical Data Science, Stanford University School of Medicine, Stanford, CA, USA.  
2. Control and Computer Engineering Department, Politecnico di Torino, Torino, Italy
3. Galatea Bio, Inc, Miami Lakes, FL, USA.
4. Department of Biomedical Sciences, School of Medicine, Universitat de Barcelona, Barcelona, Spain 
5. Creatio, School of Medicine, Universitat de Barcelona, Barcelona, Spain
6. Institute of Neurosciences, Universitat de Barcelona, Barcelona, Spain
7. Department of Biomolecular Engineering, University of California, Santa Cruz.

## Code testing
This repository includes the code developed for the manuscript:

**"Large-scale admixture mapping unveils new genetic insights into human disease"**

Due to access restrictions, reproducing the results presented in the manuscript requires access to the UK Biobank (UKBB) and All of Us datasets, which are not publicly available.  
However, we provide a **testing pipeline** that can be run on a small synthetic VCF file to validate the code structure and functionality.



---
### 1) clone the repository
clone the repo here:
```bash
git clone https://github.com/Danariki99/Admixture_mapping

```

### 2) Install the Requirements

All the codes have been executed with python:3.8.20 

Install the necessary Python packages using:

```bash
cd Admixture_mapping
pip install -r requirements.txt
pip install --no-deps scikit-allel==1.3.1

```

### 3) install plink2
Install here plink2
```bash
cd ../
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20250701.zip
unzip plink2_linux_x86_64_20250701.zip
chmod +x plink2

```

### 4) install gnomix
clone gnomix
```bash
git clone https://github.com/AI-sandbox/gnomix

```

### 5) install bcftools

```bash
    sudo apt install -y bcftools

```

### 6) install Rscript

```bash
    sudo apt-get update && apt-get install -y --no-install-recommends r-base

    Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
    Rscript -e 'BiocManager::install(version = "3.21", ask = FALSE)'
    Rscript -e 'BiocManager::install("biomaRt")'
    Rscript -e 'install.packages(c("data.table", "optparse"), repos="https://cloud.r-project.org")'
    Rscript -e 'install.packages("dbplyr", repos = "https://cloud.r-project.org")'

```


### 7) Execute the pipeline:
The pipeline includes all the steps performed after the execution of Gnomix for Local Ancestry Inference (LAI).
Since the original input files (such as VCFs and reference panels) used in the study cannot be shared, we provide a minimal example .vcf.gz file to illustrate the full pipeline structure. you can find the data folder here: https://drive.cloud.polito.it/index.php/s/mkaNL3pidDZXa7f

To run the pipeline, use the following command:

```bash
cd Admixture_mapping
./code_test.sh /path/to/data/folder path/to/your/desired/output/folder
```
Where:

- </path/to/data/folder> is the path to your data folder.

- <path/to/your/desired/output/folder> is the path to the folder you want to put results in

All results and plots will be automatically saved in the path/to/your/desired/output/folder folder. 

## Container
In case of problems reproducing the results, here we provide a guide on how to run the experiments on a Singularity container, both interactively and by executing a runscript.

## Experimental setup

Follow these steps to setup for reproducing the experiments provided in the paper
### 1) Install `Singularity` from https://docs.sylabs.io/guides/3.0/user-guide/installation.html:
	* Install `Singularity` release 3.10.2, with `Go` version 1.18.4
	* Suggestion: follow instructions provided in _Download and install singularity from a release_ section after installing `Go`
	* Install dependencies from: https://docs.sylabs.io/guides/main/admin-guide/installation.html

### 2) Clone the repository in your home folder

```bash
git clone https://github.com/Danariki99/Admixture_mapping

```

### 3) install plink2
Install here plink2
```bash
wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20250701.zip
unzip plink2_linux_x86_64_20250701.zip
chmod +x plink2

```

### 4) install gnomix
clone gnomix
```bash
git clone https://github.com/AI-sandbox/gnomix

```

### 5) Move to the source subfolder, and build the Singularity container with 
```bash
cd Admixture_mapping
sudo singularity build singularity.sif singularity.def
```
or using fake root privileges
```bash
cd Admixture_mapping
singularity build --fakeroot singularity.sif singularity.def
```

## Reproducing the analysis interactively within the Singularity container

To run testing, manually launch the Singularity container.

First of all, launch the Singularity container
```bash
singularity shell singularity.sif
```
This will run a shell within the container, and the following prompt should appear:
```bash
Singularity>
```

Be carefull, in the singularity container you will just see the folders of the direct path to the Admixture_mapping repository. If you want to put the data folder and the results folder in another path in your PC, you can bind the path in singularity running the command:

```bash
singularity shell --bind /path/to/your/folder:/linked/path/in/singularity singularity.sif
```

Now execute the whole code runnig this command:

```bash
./code_test.sh /linked/path/in/singularity/data/folder /linked/path/in/singularity/output/folder
```

be careful, the data folder and the output folders need to be inside the folder that you link into singularity, otherwise the container will not be able to see them


## Reproducing the analysis running the Singularity container

To reproduce the analysis from this paper, you can also run the Singularity container singularity container.sif in this way:

Move to the `source` folder and run the `singularity.sif` file
```bash
cd Admixture_mapping
singularity run -- bind --bind /path/to/your/folder:/linked/path/in/singularity singularity.sif /linked/path/in/singularity/data/folder /linked/path/in/singularity/output/folder
```

## Disclaimer

Since both the UK Biobank (UKBB) and All of Us datasets cannot be publicly shared, the test pipeline has been adapted to run on a small synthetic dataset, which you can find at (https://drive.cloud.polito.it/index.php/s/mkaNL3pidDZXa7f). Although the images generated by the pipeline are similar to those presented in the paper, the results should not be considered significant. Additionally, some significance thresholds have been adjusted to ensure that the pipeline can run successfully on the example data.