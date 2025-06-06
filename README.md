# ADMIXTURE MAPPING PROJECT

## abstract
In recent years, large-scale genome-wide association studies (GWAS) and sequencing projects have identified over 70,000 associations between genetic variants and human traits or diseases. However, the majority of these studies rely on participants from European ancestry dominated biobanks, limiting insights into the full spectrum of human genetic diversity and its impact on disease. Admixture mapping offers a complementary approach to GWAS by leveraging differences in disease prevalence across ancestral populations to identify risk loci for complex traits. Here, we present an extension of the snputils software suite to perform admixture mapping analysis at a large scale. We used this approach to analyze large-scale GWAS cohorts, namely the UK Biobank and the All of Us Research Program datasets, and studied associations between local ancestry and a total of 172 phenotypes. Our method not only replicates previously reported associations but also identifies novel loci with important implications for human disease. Additionally, we incorporate a fine-mapping strategy within the 31 identified regions, allowing us to genetically pinpoint the signal in 26 cases. Our findings highlight key genetic markers by considering local ancestry and quantifying their contribution to trait variability.

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

### 1) Install the Requirements

Install the necessary Python packages using:

```bash
pip install -r requirements.txt

```

### 2) Execute the pipeline:
The pipeline includes all the steps performed after the execution of Gnomix or RFMix. Since the original input files (such as VCFs and ancestry panels) used in the study cannot be shared, we provide a minimal example .msp file located in the data/ folder to illustrate the pipeline structure.

To run the test pipeline, use:

```bash
./code_test.sh data/test.msp
```

All results and plots will be automatically saved in the results/ folder.