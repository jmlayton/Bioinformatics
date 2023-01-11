# Janja
A Bioinformatic Pipeline for Processing 16S Sequence Data

This pipeline requires two different software installations for successful completion.

The first is the [Figaro software](https://github.com/Zymo-Research/figaro) developed by Zymo Biomics.

Finally, the [Dada2 package](https://benjjneb.github.io/dada2/dada-installation.html) developed by Dr. Benjamin Callahan must be installed as an R library.

**To install Janja run the following code**

```
git clone https://github.com/jmlayton/Bioinformatics.git
```

The repository contains one master directory called 'Janja.' This is where the Figaro, and Dada2 sub executables are located. These directories are also where you will find the bioinformatic outputs of each software. All sequence files, taxa.csv, and ASV_table.csv, will remain adjacent to the defined multiplexed sequencES directory. 

The directory of demultiplexed sequences must contain unmergerd forward and reverse sequences. 

The pipeline requires 6 predetermined inputs, found in the janja.par file. Set each variable as desired. The maxEE and maxN variables are used in the Dada2 software during filtering. More information can be found [here](https://benjjneb.github.io/dada2/tutorial.html) under "Filter and Trim."

To initalize and run the Janja Pipeline execute the following code:

```
chmod +x ./janja.sh
./janja.sh
```

Updates Needed:
- Silva classifier version is hardcoded. Script should instead check for latest clasifier and download it.
