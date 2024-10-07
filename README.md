<!-- vim-markdown-toc GFM -->

* [Set up environment](#set-up-environment)
* [Run](#run)
* [Basecalling](#basecalling)

<!-- vim-markdown-toc -->

Set up environment
==================

```
mamba create --yes -n dnascent-fork-detection
mamba activate dnascent-fork-detection
mamba install -n dnascent-fork-detection --yes --file requirements.txt
```

Run
===

```
snakemake --rerun-trigger mtime -p -n -j 10 \
    --default-resources "partition='nodes'" "mem='5G'" "time='4:00:00'" "cpus='4'" \
    --cluster 'sbatch -A project0015 --partition={resources.partition} --time {resources.time} --cpus-per-task={resources.cpus} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel \
    -C sample_sheet=$PWD/sample_sheet.tsv \
       species=$PWD/species.tsv \
    -d output -- dnascent
```

Basecalling
===========

Memo: To find which configuration settings you need, look in file
`final_summary_<run id>.txt` inside the run directory. This file should have
line saying which sequencer and kit was used, like:

```
protocol=sequencing/sequencing_MIN106_DNA:FLO-MIN106:SQK-LSK109
```

```
nohup snakemake -s guppy.smk -p -j 1 \
    -C guppy_dir=/export/projects/III-data/wcmp_bioinformatics/db291g/applications/ont-guppy/ont-guppy-gpu_6.4.2_linux64 \
       guppy_config=dna_r9.4.1_450bps_fast.cfg \
       run_dirs='Gab_dnacent_60m_15_2_23/no_sample/20230215_1415_MN23371_FAU06237_e2778d89
        Gab_DNAcent_60m_b_7_3_23/no_sample/20230307_1205_MN23371_FAU92692_b96dfc11
        Gab_dnacent_1mM_28_11_22/no_sample/20220906_0616_MC-114086_FAU08788_7bbffcec
        Gab_DNAcent_25_5_23/Gabriel_DNAcent_LeishTC_25_5_23/no_sample/20230525_1145_X1_FAW63215_0877a6bf'
```
