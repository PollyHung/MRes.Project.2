## Difficulty in Installation of MAJIQ     

### Why download MAJIQ?    
We cannot identify any signficant alternative splicing event using SUPPA2, upon querying [assessment papers](https://academic.oup.com/bib/article/21/6/2052/5648232) for differential splicing analysis we found that under the same dataset, althoguh SUPPA2 is as effective as MAJIQ for identifying qPCR confirmed AS events (i.e., these two packages are most **accurate** in identifying biological existing AS events. SUPPA2 identifies far less AS events compared to most other available packages on the market, and their FDR is the highest for some unknown reason? Like out of all the analysis performed they produces the highest FDR for AS events. Thus statistically challenging for subsequent analysis. SUPPA2 is recommended by PacBio and like SUPPA2, MAJIQ is also an event based differential AS event analysis package. Therefore we thought to substitute SUPPA2 with MAJIQ. In addition, we will also perform exon-based differential AS event analysis using DEXSeq.        

### How to install?     
MAJIQ instlalation requires the following lib/software to be installed in your system.
1. Python version 3.10; Lower or higher versions may work, but 3.10 is the currently tested and supported version        
2. C++11 compiler with openMP. GCC includes that by default, but clang can be updated to include those (Clang/OMP). MAJIQ/VOILA has been tested to work with GNU GCC>=7.2, RedHat GCC>=4.7.2.      
3. HTSlib library. This is a C library for reading/writing high-throughput sequencing data developed by Samtools organisation. MAJIQ installation assumes the library and its header files are present in the Unix default locations (/usr/local/lib, /usr/local/include). If that is not the case the appropiate locations can be specified setting the following environment variables.

*DO NOT TRY ON LOCAL COMPUTER, USE HPC INSTEAD*     
Step 1: activate an environment, python version 3.11.6      
```
python3 -m venv env
source env/bin/activate
```
Step 2: download and install latest htslib     
```
curl -OL https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar -xf htslib-1.13.tar.bz2  # extract archive
cd htslib-1.13  # change into htslib source directory
./configure --prefix=$HOME/install/htslib-1.13
make
make install
```
Step 3: download the majiq packages from the bitbucket website    
*do not directly use their recommended code pip install git+https://bitbucket.org/biociphers/majiq_academic.git, it won't work*
From https://bitbucket.org/biociphers/majiq_academic/downloads/ download the biociphers-majiq_academic-8ff62a93f786.zip
download, unzip, and rename the folder to majiq_academic
```
cd ~/MRes.project.2/codes/0_packages/
wget https://bitbucket.org/biociphers/majiq_academic/get/8ff62a93f786.zip
unzip biociphers-majiq_academic-8ff62a93f786.zip
mv biociphers-majiq_academic-8ff62a93f786 majiq_academic
```
Step 4: Install the majiq packages     
remember prior to installation you have to export the htslib pathway to the package you installed in step 1. 
```
cd majiq_academic/
export HTSLIB_LIBRARY_DIR=$HOME/install/htslib-1.13/lib
export HTSLIB_INCLUDE_DIR=$HOME/install/htslib-1.13/include
pip install .
```

### Academic License 
Starting with MAJIQ version 2.5, academic and commercial versions require providing a license file for use. Please note that the software will not function without providing the license file.            
In order to support usage of the license with minimal disruption to various workflow use cases, we provide many methods to provide the license when running majiq or voila:      
1. Provide the switch --license to majiq or voila with an explicit path to the license file    
2. Set the environment variable MAJIQ_LICENSE_FILE with the explicit path to the license file    
3. Check the current working directory for any file which begins with "majiq_license"    
4. Check the users's home directory for any file which begins with "majiq_license"    
License is downloaded to the `majiq_academic` folder in the name of `majiq_license_academic_official.lic`      


