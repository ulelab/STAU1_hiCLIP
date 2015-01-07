mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/hiCLIP
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/ribosome_profiling
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/mRNASeq
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/hnRNPc
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/STAU1_iCLIP

mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/hiCLIP
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/ribosome_profiling
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/mRNASeq
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/hnRNPc
mkdir -p /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/STAU1_iCLIP

cp /teraraid/yoichiro/LUs27_HiCLIP/196nts/rowReads/LUs27_196.fq /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/hiCLIP/
cp /teraraid/yoichiro/LUs27_HiCLIP/196nts/barcode/LUs27_196.barcode /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/hiCLIP/

# hnRNPc data: 20110307_LUjh25 id 5 iCLIP_HnRNPC_Hela_Rt9Clip_hu_NNNTGGCNN_20110307_LUjh25_5.fq.gz
cp /teraraid/yoichiro/transcriptome_mapping/hnRNPc/rawReads/*.fq /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/hnRNPc/
cp /teraraid/yoichiro/transcriptome_mapping/hnRNPc/barcode/*.barcode /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/hnRNPc/

cp /teraraid/yoichiro/transcriptome_mapping/STAU1_iCLIP/rawReads/*.fq /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/STAU1_iCLIP/
cp /teraraid/yoichiro/transcriptome_mapping/STAU1_iCLIP/barcode/*.barcode /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/STAU1_iCLIP/

cp /teraraid/yoichiro/transcriptome_mapping/mRNASeq/rowReads/*.fq /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/mRNASeq/
cp /teraraid/yoichiro/transcriptome_mapping/mRNASeq/barcode/*.barcode /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/mRNASeq/

cp /teraraid/yoichiro/transcriptome_mapping/RibosomeProfiling/rowReads/*.fq /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/fastq/ribosome_profiling/
cp /teraraid/yoichiro/transcriptome_mapping/RibosomeProfiling/barcode/*.barcode /netscr/yoichiro/STAU1_hiCLIP/data/unprocessed/barcode/ribosome_profiling/

qsub -j y -b y -pe smp 8 ~/Python-2.7.1/python /netscr/yoichiro/STAU1_hiCLIP/exec/run_analysis.py -d /netscr/yoichiro/STAU1_hiCLIP -p ~/Python-2.7.1 -r ~/R-2.15.1




