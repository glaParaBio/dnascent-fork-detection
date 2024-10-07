import sys
import shutil
import pandas
import numpy
import glob
from pathlib import Path
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()

EVENTS = ['origin', 'leftFork', 'rightFork', 'termination']

LmjF_CHROMS = ['LmjF.%02d' % x for x in range(1, 37)]

def srcdirn(x):
    """We need this because of
    https://bitbucket.org/snakemake/snakemake/issues/1132/different-strings-pointing-to-the-same
    """
    return os.path.normpath(srcdir(x))

def find_fastq(directory):
    fx = list(Path(directory).rglob('*.fastq'))
    gz = list(Path(directory).rglob('*.fastq.gz'))
    tmp = [re.sub('\.gz$', '', str(x)) for x in gz]
    if len(fx + tmp) != len(set(fx + tmp)):
        raise Exception('There are duplicate fastq files (compressed and uncompressed) in directory %s' % directory)
    return fx + gz

def get_fastq_for_sample(ss, sample_id):
    dat = ss[(ss.sample_id == sample_id) & (pandas.isna(ss.fastq_dir) == False)]
    fq = []
    for i, row in dat.iterrows():
        d = os.path.join(row.run_dir, row.fastq_dir)
        fx = find_fastq(d)
        for x in fx:
            fq.append(os.path.normpath(x))
    assert len(fq) > 0
    return sorted(set(fq))

def get_fastq_for_experiment_id(wc):
    dat = ss[ss.experiment_id == wc.experiment_id]
    fq = []
    for i, row in dat.iterrows():
        d = os.path.join(row.run_dir, row.fastq_dir)
        fx = find_fastq(d)
        for x in fx:
            fq.append(os.path.normpath(x))
    assert len(fq) > 0
    assert len(fq) == len(set(fq))
    return sorted(fq)

ss = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#')
ss = ss[ss.use == 'yes']
run_id = []
for idx,row in ss.iterrows():
    if pandas.isna(row.run_dir):
        run_id.append(numpy.nan)
    else:
        rid = os.path.basename(row.run_dir) + '_' + re.sub('/', '_', row.fastq_dir) 
        # assert rid not in run_id
        run_id.append(rid)
ss['run_id'] = run_id

species_ss = pandas.read_csv(config['species'], sep= '\t', comment= '#')
assert all(x in list(species_ss.species) for x in ss.species)

sample_ss = ss[['sample_id', 'species', 'use']].drop_duplicates()

wildcard_constraints:
    experiment_id = '|'.join([re.escape(x) for x in ss[pandas.isna(ss.experiment_id) == False].experiment_id.unique()]),
    pair_id = '|'.join([re.escape(x) for x in ss[pandas.isna(ss.pair_id) == False].pair_id.unique()]),
    sample_id = '|'.join([re.escape(x) for x in ss.sample_id.unique()]),
    species = '|'.join([re.escape(x) for x in ss.species.unique()]),
    chrom = '|'.join([re.escape(x) for x in LmjF_CHROMS]),
    event =  '|'.join([re.escape(x) for x in EVENTS]),

include: 'workflows/dnascent.smk'


rule dnascent:
    input:
        expand('{species}/dnascent/{sample_id}.fork.bed', zip, sample_id= ss.sample_id, species= ss.species),
        expand('{species}/dnascent/{sample_id}.forkSense.bedgraph.gz', zip, sample_id= ss[ss.sample_type.isin(['BrdU', 'ctrl'])].sample_id, species= ss[ss.sample_type.isin(['BrdU', 'ctrl'])].species),
        expand('{species}/dnascent/{sample_id}.brdu.bedgraph.gz', zip, sample_id= ss[ss.sample_type.isin(['BrdU', 'ctrl'])].sample_id, species= ss[ss.sample_type.isin(['BrdU', 'ctrl'])].species),
        os.path.join(workflow.basedir, 'results/event_summary.tsv'),


rule download_fasta:
    output:
        fa= 'ref/{species}_Genome.fasta',
        fai= 'ref/{species}_Genome.fasta.fai',
    params:
        url= lambda wc: species_ss[species_ss.species == wc.species].fasta.iloc[0],
        gunzip = lambda wc: '| gunzip' if species_ss[species_ss.species == wc.species].fasta.iloc[0].endswith('.gz') else '',
    shell:
        r"""
        curl -L -s {params.url} {params.gunzip} > {output.fa}
        samtools faidx {output.fa}
        """


rule download_gff:
    output:
        gff= 'ref/{species}.gff',
    params:
        url= lambda wc: species_ss[species_ss.species == wc.species].gff.iloc[0],
        gunzip = lambda wc: '| gunzip' if species_ss[species_ss.species == wc.species].gff.iloc[0].endswith('.gz') else '',
    shell:
        r"""
        curl -L -s {params.url} {params.gunzip} > {output.gff}
        """

rule nanoplot:
    input:
        fq= get_fastq_for_experiment_id, 
    output:
        tsv= 'nanoplot/{experiment_id}.NanoStats.tsv',
        txt= 'nanoplot/{experiment_id}.NanoStats.txt',
    run:
        with open(output.tsv + '.files', 'w') as fq:
            for x in input.fq:
                fq.write(x + ' ')
        
        shell(r"""
            fq=`cat {output.tsv}.files`
            NanoPlot --fastq_rich $fq --prefix {wildcards.experiment_id}. --outdir nanoplot
            sed 's/: \+/:\t/; s/\s\+$//' {output.txt} \
            | awk -v FS='\t' '{{if(NF == 1){{print $0 "\t."}}else{{print $0}}}}' > {output.tsv}
            """)

        os.remove(output.tsv + '.files')


rule minimap2:
    input:
        fa= 'ref/{species}_Genome.fasta',
        fq= lambda wc: get_fastq_for_sample(ss, wc.sample_id),
    output:
        bam= '{species}/minimap2/{sample_id}.bam',
        bai= '{species}/minimap2/{sample_id}.bam.bai',
    run:
        with open(output.bam + '.files', 'w') as fq:
            for x in input.fq:
                fq.write(x + ' ')

        shell(r"""
            fq=`cat {output.bam}.files`
            minimap2 --MD -R '@RG\tID:{wildcards.sample_id}\tSM:{wildcards.sample_id}' -a -x map-ont -t 12 {input.fa} $fq \
            | samtools sort > {output.bam}
            samtools index {output.bam}
            """)

        os.remove(output.bam + '.files')

rule bam_stat:
    input:
        bam= '{species}/minimap2/{sample_id}.bam',
    output:
        stat= '{species}/minimap2/{sample_id}.stat',
    shell:
        r"""
        samtools stat -@ 4 {input.bam} > {output.stat}
        """
