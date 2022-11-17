#!/bin/bash

SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
Figaro_dir=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

cd $Figaro_dir
cd ../..
janja_dir=$(pwd)

echo "--------------------Beginning Figaro Script ----------------------------" | tee -a $janja_dir/janja_error.out 

exec 2>> $janja_dir/janja_error.out 1>> $janja_dir/janja_output.out

# cd
# source ~/miniconda3/etc/profile.d/conda.sh

# wget http://john-quensen.com/wp-content/uploads/2020/03/figaro.yml
# conda env create -n figaro -f figaro.yml

# wget https://github.com/Zymo-Research/figaro/archive/master.zip
# unzip master.zip
# rm master.zip
# mv figaro-master figaro
# cd figaro
# chmod 755 figaro.py


# Activate the FIGARO environment
source ~/miniconda3/etc/profile.d/conda.sh # If necessary for your conda installation.
conda activate figaro

cd $Figaro_dir
cd ../..
janja_dir=$(pwd)

source  $janja_dir/janja.par
cd $Multiplexed_Seqs_Directory
cd ../
Parent=$(pwd)
Demultiplexed_Seqs_Dir=$(pwd)/Demultiplexed_Seqs

rm -r $Figaro_dir/figaro_out; mkdir $Figaro_dir/figaro_out

cd $Figaro_dir/figaro_out

 # Run FIGARO
 # cd to installation folder
 cd ~/figaro/figaro
 python figaro.py -i $Demultiplexed_Seqs_Dir -o $Figaro_dir/figaro_out \
    -f $Forward_Primer_Length -r $Reverse_Primer_Length -a $Amplicon_Length -F illumina

F_R=$(sed 'x;$!d' <$janja_dir/janja_output.out); Forward_Trim=$(cut -b 19-21 <<< $F_R); Reverse_Trim=$(cut -b 24-26 <<< $F_R)

echo "Forward_Trim_Length=$Forward_Trim" | tee $Figaro_dir/figaro_out/Forward_Reverse_Trim.par
echo "Reverse_Trim_Length=$Reverse_Trim" | tee -a $Figaro_dir/figaro_out/Forward_Reverse_Trim.par

conda deactivate

echo "--------------------End of Figaro Script ----------------------------" | tee -a $janja_dir/janja_error.out