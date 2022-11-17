#!/bin/bash

SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
Qiime2_dir=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

cd $Qiime2_dir
cd ../..
janja_dir=$(pwd)

echo "--------------------Beginning Qiime 2 Demultiplexing Script ----------------------------" | tee -a $janja_dir/janja_error.out 

exec 2>> $janja_dir/janja_error.out 1>> $janja_dir/janja_output.out

eval "$(conda shell.bash hook)"
conda activate qiime2-2022.2

source $janja_dir/janja.par

cd $Multiplexed_Seqs_Directory

Sequence_Batch=$Multiplexed_Seqs_Directory #folder of forward.fastq.gz, reverse.fastq.gz, barcodes.fastq.gz, metadata.txt   
mkdir -p $Multiplexed_Seqs_Directory/temp1
mv *.fastq.gz $Multiplexed_Seqs_Directory/temp1/

input=$Multiplexed_Seqs_Directory/temp1

mapfile=$Sequence_Batch/metadata.txt 

cd ../

Parent=$(pwd)
dataset=$(pwd)/Demultiplexed_Seqs

rm -r $(pwd)/Demultiplexed_Seqs

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path $input \
  --output-path $input'.qza' 

qiime demux emp-paired \
  --i-seqs $input'.qza' \
  --m-barcodes-file $mapfile \
  --m-barcodes-column BarcodeSequence \
  --o-per-sample-sequences $Parent/"demux.qza" \
  --p-no-golay-error-correction \
  --output-dir $dataset

qiime tools export \
  --input-path $Parent/"demux.qza" \
  --output-path $dataset/ 

cd $Multiplexed_Seqs_Directory

mv $input/*.fastq.gz $Multiplexed_Seqs_Directory/
rm -r $Multiplexed_Seqs_Directory/temp1; rm -r $Multiplexed_Seqs_Directory/temp1.qza 
rm $Parent/demux.qza


echo "--------------------End of Qiime 2 Demultiplexing Script ----------------------------" | tee -a $janja_dir/janja_error.out
