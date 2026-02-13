#!/bin/bash
python -m venv .venv
source .venv/bin/activate # deactivate

pip install -U pip
pip install -e . # .toml 존재하는 위치에서
# pip uninstall -y trias

trias bed-prepare \
  --bed hg38.bed.gz \
  --fasta /home/sycw12475/01_projects/241204_IlBinKim_ASD/00_public_data/Homo_sapiens_assembly38.fasta \
  --outdir ./01_hg38 \
  --threads 70 \
  --chunk 60000 \
  --motif-col 10 \
  --has-header

trias bed-prepare \
  --bed T2T.bed.gz \
  --fasta /home/junkim/01_projects/00_reference_genome/230823_CHM13_v2.0/chm13v2.0.fa \
  --outdir ./02_T2T \
  --threads 70 \
  --chunk 60000 \
  --motif-col 10 \
  --mapping t2t_chr_name.txt \
  --has-header
