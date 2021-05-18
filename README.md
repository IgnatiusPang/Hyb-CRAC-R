# Hyb-CRAC-R

Analysis of Cross-linking, ligation and sequencing of hybrids (CLASH) experimental data  to identify RNA-RNA interactions using the pyCRAC (https://git.ecdf.ed.ac.uk/sgrannem/pycrac) and hyb (https://github.com/gkudla/hyb) pipelines and custom R scripts

This pipeline was designed to run on the Katana PBS cluster at UNSW. It is most likely that modifications to the scripts are required (e.g. removal of or edits to the 'module' commands) to run on another compute platform. 

This repository is https://github.com/IgnatiusPang/Hyb-CRAC-R

# Installation of Software Prerequisites

## bbtools
To install bbmerge with conda run:
```{bash}
conda install -c agbiome bbtools
```

## FastQC
To install FastQC with conda run one of the following:
```{bash}
conda install -c bioconda fastqc
conda install -c bioconda/label/broken fastqc
conda install -c bioconda/label/cf201901 fastqc
```

## MultiQC
To install MultiQC with conda run one of the following:
```{bash}
conda install -c bioconda multiqc
conda install -c bioconda/label/cf201901 multiqc
```

      module load python/2.7.15 

## blat

To install blat with conda run one of the following:
```{bash}
conda install -c bioconda blat
conda install -c bioconda/label/cf201901 blat
```

## Blast+
To install this package with conda run one of the following:
```{bash}
conda install -c bioconda blast
conda install -c bioconda/label/cf201901 blast
```

## Bowtie2
To install Bowtie2 with conda run one of the following:
```{bash}
conda install -c bioconda bowtie2
conda install -c bioconda/label/broken bowtie2
conda install -c bioconda/label/cf201901 bowtie2
```

## Novoalign
For installation on how to install novoaligner please refer to the following:
http://www.novocraft.com/support/download/

One such verion for Linux 3.1 may be downloaded using wget and subsequently unzipped: 
```{bash}
wget http://www.novocraft.com/support/download/download.php?filename=V4.03.01/novocraftV4.03.01.Linux3.10.0.tar.gz
tar -xvf novocraftV4.03.01.Linux3.10.0.tar.gz
```

## fastx
To install this package with conda run one of the following:
```{bash}
conda install -c bioconda fastx_toolkit
conda install -c bioconda/label/cf201901 fastx_toolkit
```

## pblat
To install this package with conda run:
```{bash}
conda install -c bioconda pblat
```

## flexbar
To install this package with conda run one of the following:
```{bash}
conda install -c bioconda flexbar
conda install -c bioconda/label/cf201901 flexbar
```

## perl
To install this package with conda run one of the following:
```{bash}
conda install -c conda-forge perl
conda install -c conda-forge/label/gcc7 perl
conda install -c conda-forge/label/broken perl
conda install -c conda-forge/label/cf201901 perl
conda install -c conda-forge/label/cf202003 perl
```

## viennarna
To install this package with conda run one of the following:
```{bash}
conda install -c bioconda viennarna
conda install -c bioconda/label/cf201901 viennarna
```

## unafold
```{bash}
conda install -c bioconda oligoarrayaux
```

## hyb

Follow installation instructions from this GitHub repository: https://github.com/gkudla/hyb

For example run the script:
```{bash}
./INSTALL
```

## pyCRAC installation


* pyCRAC is installed locally using a virual environment. It's installed in the dir

```{bash}
my_python_dir/
```

* The installation required that we remove python from .bashrc setup so that python/2.7.10 can be loaded - the virualenv command has changed in python3 so after restarting the katana session

```{bash}
mkdir my_python_dir
cd my_python_dir
module load python/2.7.15
virtualenv --system-site-packages my_python_env
. my_python_env/bin/activate           #The environment needs to be activated
wget https://bitbucket.org/sgrann/pycrac/get/cc0dcde5ea3f.zip
unzip cc0dcde5ea3f.zip
cd sgrann-pycrac-cc0dcde5ea3f/
python setup.py install
```

* Have had to install some additional packages: 

```{bash}
wget http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip
unzip pp-1.6.4.zip
cd pp-1.6.4
python setup.py install
```

* pip seems to be a better way to install the following

```{bash}
pip install greenlet
pip install curtsies
pip install requests
pip install pygments
pip install ruffus --upgrade
``

* py ReadCoutners seems to work now
* to exit the session

```{r}
deactivate
```


# Usage Notes:

## Run Novoindex to create an index for novoaligner:
```{r}
./novoindex Hyb-CRAC-R/Data/Genomic/JKD6009/Sa_JKD6009.nix \
  Hyb-CRAC-R/Data/Genomic/JKD6009/Sa_JKD6009.fasta
```

## Run example analysis

* To run an example analysis, we have to do the following: 

1. Open the Snakemake file "Source/Demultiplex/Script_Per_Dataset/dm_miseq_github.smk", change the necessary parameters inside this file.
2. Inside this file change the CRAC_ENVIRONMENT_ACTIVATE variable to point to your installation of the pyCRAC installation:   (For example: CRAC_ENVIRONMENT_ACTIVATE=/srv/scratch/z1234567/my_python_dir/my_python_env/bin/activate) 
3. The file "Source/Demultiplex/Script_Per_Dataset/run_smk_dm_miseq_github.sh" contains the script that runs the snakemake command. Remove any 'module' command if needed. Configure the snakemake command parameters (e.g. remove the --cluster-config and --cluster parameter if not using PBS queue, remove the --use-envmodules parameter if not using any PBS 'module'.)
4. The "Source/Demultiplex/Script_Per_Dataset/common_rules.smk" file contains other snakemake rules. Remove any 'module' command if needed. Correct the path to executables and referenced scripts as required. 
5. Run the script "bash Source/Demultiplex/Script_Per_Dataset/run_smk_dm_miseq_github.sh"
6. Reference sequences such as Barcodes and Adaptors are included in the folder (Hyb-CRAC-R/Data/Ref_Seq)
7. For more information on how to run the analysis please see the following set of slides: "Hyb-CRAC-R usage information.pdf"
8. Good luck!
