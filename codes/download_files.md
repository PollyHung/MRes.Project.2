## File Download Details    
This is a script showing how and were packages and reference files are acquired, most actions took place in command line, but if otherwise it will be stated. 

### Download Packages 
#### SMRTLINK v13 
Website: https://www.pacb.com/support/software-downloads/   
```
wget https://downloads.pacbcloud.com/public/software/installers/smrtlink-release-sequel2_13.1.0.221970.zip
unzip smrtlink-release-sequel2_13.1.0.221970.zip
chmod +x smrtlink-release-sequel2_13.1.0.221970_linux_x86-64_libc-2.17_anydistro.run
./smrtlink-release-sequel2_13.1.0.221970_linux_x86-64_libc-2.17_anydistro.run --rootdir smrtlink --smrttools-only
```

#### BGI functions
Website: https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline/tree/master for bulk sequencing
Website: https://github.com/shizhuoxing/scISA-Tools/tree/master for single cell sequencing
```
git clone https://github.com/shizhuoxing/BGI-Full-Length-RNA-Analysis-Pipeline.git
mv BGI-Full-Length-RNA-Analysis-Pipeline bgi_commands
ls bgi_commands/bin
## here is where the commands are located
## classify_by_primer.fullpa.pl  classify_by_primer.pl  flnc2sam.pl  PolymeraseReads.stat.pl  SubReads.stat.pl
```

#### SQUANTI3
Website: https://github.com/ConesaLab/SQANTI3
```
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.2.1.tar.gz
tar -xf v5.2.1.tar.gz
cd SQANTI3-5.2.1
conda env create -f SQANTI3.conda_env.yml
conda activate SQANTI3.env
## once activated the virtual environemnt, you should see prompt changing to:
(SQANTI3.env)$
```


