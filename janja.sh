#!/bin/bash
SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
Janja_master_dir=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

echo "--------------------Beginning Janja Script ----------------------------" | tee $Janja_master_dir/janja_error.out $Janja_master_dir/janja_output.out

exec 2>> $Janja_master_dir/janja_error.out 1>> $Janja_master_dir/janja_output.out

chmod +x $Janja_master_dir/Janja/Qiime2/demux.sh

cd $Janja_master_dir/Janja/Qiime2
./demux.sh

chmod +x $Janja_master_dir/Janja/Figaro/figaro.sh

cd $Janja_master_dir/Janja/Figaro
./figaro.sh

#Dada2 configuration file build
rm -r $Janja_master_dir/Janja/Dada2/Dada2_output; mkdir -p $Janja_master_dir/Janja/Dada2/Dada2_output
source  $Janja_master_dir/janja.par

#Write configuration from janaj.par
echo "f_primer_len <- $Forward_Primer_Length" | tee $Janja_master_dir/Janja/Dada2/Dada2_Config.R 
echo "r_primer_len <- $Reverse_Primer_Length" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R 
echo "maxEE <- $maxEE" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R
echo "maxN <- $maxN" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R

#Write configuration from Forward_Reverse_Trim.par
source $Janja_master_dir/Janja/Figaro/figaro_out/Forward_Reverse_Trim.par
echo "f_truncLen <- $Forward_Trim_Length" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R
echo "r_truncLen <- $Reverse_Trim_Length" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R
echo "Janja_master_dir <- '$Janja_master_dir'" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R
echo "Multiplexed_Seqs_Directory <-  '$Multiplexed_Seqs_Directory'" | tee -a $Janja_master_dir/Janja/Dada2/Dada2_Config.R

cd $Janja_master_dir/Janja/Dada2
Rscript dada2.R

echo "--------------------End of Janja Script ----------------------------" | tee -a $Janja_master_dir/janja_error.out