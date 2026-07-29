# ADMIXTURE MAPPING PROJECT

## abstract
Genome-wide association studies have successfully identified thousands of genetic associations, yet their predominant reliance on European-descent populations limits insights into the full spectrum of human genetic diversity and its impact on disease. Admixture mapping offers a powerful, complementary approach by leveraging differences in haplotype frequencies across ancestral backgrounds to identify risk loci for complex traits. Here, we perform a large-scale, multi-ancestry admixture mapping study across 415,792 unrelated individuals in the UK Biobank, examining associations between local haplotype ancestry and 108 phenotypes. Our approach identifies 13 genome-wide significant ancestry-phenotype associations, recovering previously reported signals while uncovering four novel ancestry-associated findings, including new risk loci for atrial fibrillation, dermatitis, and angina pectoris. To overcome the limited resolution of traditional admixture mapping, we implemented a conditional fine-mapping framework, which enabled us to localize four putatively causal variants. In silico variant effect prediction and eQTL integration revealed regulatory and missense effects predominantly localized to lung, and immune tissues, aligning with captured phenotypes such as asthma, dermatitis, and hypothyroidism. Notably, our findings demonstrate striking genetic heterogeneity, revealing how the same clinical phenotype can arise through distinct genetic pathways depending on the ancestral background. Overall, this work highlights the critical importance of modeling local ancestry structure to refine genetic associations, uncover novel disease mechanisms, and improve the equitable translation of genomic medicine.

## Authors 
Riccardo Smeriglio<sup>1,a</sup>, Sonia Moreno-Grau<sup>2,3,a,b</sup>, Daniel Mas Montserrat<sup>2</sup>, Guhan Venkataraman<sup>2</sup>, David Bonet<sup>2,4-8</sup>, Caterina Fuses<sup>5-8</sup>, Manuel A. Rivas<sup>2</sup>, Alessandro Savino<sup>1</sup>, Stefano Di Carlo<sup>1</sup>, Jordi Abante<sup>5-8,b</sup>, Alexander G. Ioannidis<sup>2,4,9,b</sup>


1. Control and Computer Engineering Department, Politecnico di Torino, Torino, Italy 
2. Department of Biomedical Data Science, Stanford University School of Medicine, Stanford, CA, USA
3. Faculty of Health Sciences, Universidad Europea de Valencia, Valencia, Spain
4. Genomics Institute, University of California, Santa Cruz, Santa Cruz, CA, USA
5. Department of Biomedical Sciences, School of Medicine, Universitat de Barcelona, Barcelona, Spain
6. Institute of Neurosciences, Universitat de Barcelona, Barcelona, Spain 
7. Institut d'Investigacions Biomèdiques August Pi i Sunyer (IDIBAPS), Barcelona, Spain
8. Centro de Investigación Biomédica en Red Sobre Enfermedades Neurodegenerativas (CIBERNED), Instituto de Salud Carlos III, Madrid, Spain
9. Institute for Computational and Mathematical Engineering, Stanford University, Stanford, CA, USA  
a. These authors contributed equally  
b. Corresponding authors


## Code testing
This repository includes the code developed for the manuscript:

**"Multi-ancestry admixture mapping reveals ancestry-associated disease loci in the UK Biobank"**

Due to access restrictions, reproducing the results presented in the manuscript requires access to the UK Biobank (UKBB) which is not publicly available.  
However, we provide a **testing pipeline** that can be run on a small synthetic VCF file to validate the code structure and functionality.



---
### 1) clone the repository
clone the repo here:
```bash
git clone https://github.com/Danariki99/Admixture_mapping

```

### 2) Install the Requirements

All the codes have been executed with python: 3.8.20 

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

### 4) install RFMix (local ancestry inference — default)
RFMix v2 is the default LAI step of the pipeline. Clone and compile it from
source **as a sibling of `Admixture_mapping/`** (needs `gcc/g++`, `make`,
`autoreconf`):
```bash
cd ../
git clone https://github.com/slowkoni/rfmix.git
cd rfmix
autoreconf --force --install
./configure
make
rm -rf .git          # keep only the files/binary, no nested git repo
./rfmix               # sanity check: prints "RFMIX v2.03 ..."
cd ../Admixture_mapping
```
This produces the `rfmix/rfmix` binary that the LAI wrapper
(`LAI/rfmix_test.py`) calls via the relative path `../rfmix/rfmix`.

### 5) install bcftools and tabix
`tabix` is required to index the per-chromosome VCFs for RFMix.
```bash
    sudo apt install -y bcftools tabix

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
The pipeline includes all the steps performed.
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

### 4) install RFMix (default LAI)
Clone and compile RFMix as a sibling of the repository (the pipeline calls it via
`../rfmix/rfmix`):
```bash
git clone https://github.com/slowkoni/rfmix.git
cd rfmix && autoreconf --force --install && ./configure && make && rm -rf .git && cd ..

```

### 5) Move to the `Admixture_mapping` subfolder, and build the Singularity container with 
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

To reproduce the analysis from this paper, you can also run the `singularity.sif`
container directly (its runscript calls `code_test.sh` with the two arguments):

Move to the `Admixture_mapping` folder and run the `singularity.sif` file
```bash
cd Admixture_mapping
singularity run --bind /path/to/your/folder:/linked/path/in/singularity singularity.sif /linked/path/in/singularity/data/folder /linked/path/in/singularity/output/folder
```

## Disclaimer

Since the UK Biobank (UKBB) cannot be publicly shared, the test pipeline has been adapted to run on a small synthetic dataset, which you can find at (https://drive.cloud.polito.it/index.php/s/mkaNL3pidDZXa7f). Although the images generated by the pipeline are similar to those presented in the paper, the results should not be considered significant. Additionally, some significance thresholds have been adjusted to ensure that the pipeline can run successfully on the example data.