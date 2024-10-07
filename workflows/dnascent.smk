rule dnascent_install:
    output:
        'DNAscent/bin/DNAscent',
    shell:
        r"""
        rm -rf DNAscent
        git clone --recursive https://github.com/MBoemo/DNAscent
        cd DNAscent
        git checkout 7e4be09
        make
        """


rule prepare_sequencing_summary:
    output:
        seq_sum= temp('{run_id}.sequencing_summary.tmp'),
    run:
        sample_dat = ss[(ss.run_id == wildcards.run_id) & (pandas.isna(ss.fast5_dir) == False) & (pandas.isna(ss.fastq_dir) == False)][['run_dir', 'fastq_dir', 'fast5_dir']].drop_duplicates()
        assert len(sample_dat) == 1
        run_dir = sample_dat.run_dir.iloc[0]
        fastq_dir = sample_dat.fastq_dir.iloc[0]

        # sum_file = glob.glob(os.path.join(run_dir, 'sequencing_summary_*.txt'))
        sum_file = list(set(glob.glob(os.path.join(run_dir, fastq_dir, 'sequencing_summary*.txt'))))
        assert len(sum_file) == 1
        sum_file = sum_file[0]

        # Get actual fast5 files - there should be no need for doing this!
        fast5 = glob.glob(os.path.join(run_dir, sample_dat.fast5_dir.iloc[0], '*.fast5'))
        fast5 = [os.path.basename(x) for x in fast5]

        # usecols = ['filename_fast5', 'read_id']
        usecols = ['filename', 'read_id']
        sum_out = pandas.read_csv(sum_file, sep= '\t', usecols= usecols)[usecols]
        sum_out = sum_out[sum_out.filename.isin(fast5)]
        sum_out.to_csv(output.seq_sum, sep= '\t', index= False)

def get_fast5_dir(wc):
    dat = ss[(ss.run_id == wc.run_id) & (pandas.isna(ss.fast5_dir) == False)][['run_dir', 'fast5_dir']].drop_duplicates()
    assert len(dat) == 1
    fast5 = os.path.abspath(os.path.join(dat.run_dir.iloc[0], dat.fast5_dir.iloc[0]))
    return fast5

rule dnascent_index:
    input:
        dnascent='DNAscent/bin/DNAscent',
        seq_sum= '{run_id}.sequencing_summary.tmp',
        fast5_dir= get_fast5_dir,
    output:
        idx= temp('{species}/dnascent/{run_id}/index.dnascent'),
    shell:
        r"""
        LD_LIBRARY_PATH=${{CONDA_PREFIX}}/lib \
        {input.dnascent} index --output {output.idx} --files {input.fast5_dir} --sequencing-summary {input.seq_sum}
        """


rule dnascent_concat_index:
    input:
        idx= lambda wc: temp(expand('{{species}}/dnascent/{run_id}/index.dnascent', run_id=ss[ss.sample_id == wc.sample_id].run_id.unique())),
    output:
        idx= temp('{species}/dnascent/{sample_id}/index.dnascent'),
    shell:
        r"""
        cat {input.idx} | grep '^#' | uniq > {output.idx}
        cat {input.idx}  | grep -v '^#' >> {output.idx}
        """


rule dnascent_detect:
    input:
        dnascent='DNAscent/bin/DNAscent',
        idx= '{species}/dnascent/{sample_id}/index.dnascent',
        bam= '{species}/minimap2/{sample_id}.bam',
        ref= 'ref/{species}_Genome.fasta',
    output:
        detect= '{species}/dnascent/{sample_id}.detect',
    params:
        mapq=20,
        rlen=1000,
        gpu=lambda wc: ss[ss.sample_id == wc.sample_id].gpu.iloc[0],
    resources:
        partition='gpu',
        mem='20G',
        time='96:00:00',
        cpus='20',
    shell:
        r"""
        LD_LIBRARY_PATH=${{CONDA_PREFIX}}/lib \
        {input.dnascent} detect --bam {input.bam} --reference {input.ref} \
            --index {input.idx} --output {output.detect} --threads {resources.cpus} --GPU {params.gpu} \
            --quality {params.mapq} --length {params.rlen}
        """

rule dnascent_forksense:
    input:
        dnascent='DNAscent/bin/DNAscent',
        detect= '{species}/dnascent/{sample_id}.detect',
    output:
        forkSense= '{species}/dnascent/{sample_id}/dnascent.forkSense',
        bed= ['{species}/dnascent/{sample_id}/leftForks_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/origins_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/rightForks_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/terminations_DNAscent_forkSense.bed'],
    shell:
        r"""
        module load compilers/gcc/13.1.0 || true
        pwd=$PWD
        cd `dirname {output.forkSense}`
        LD_LIBRARY_PATH=${{CONDA_PREFIX}}/lib \
        $pwd/{input.dnascent} forkSense --detect $pwd/{input.detect} --output `basename {output.forkSense}` --markOrigins --markTerminations --markForks -t 20
        """

rule cat_fork_beds:
    input:
        bed= ['{species}/dnascent/{sample_id}/leftForks_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/origins_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/rightForks_DNAscent_forkSense.bed',
              '{species}/dnascent/{sample_id}/terminations_DNAscent_forkSense.bed'],
    output:
        bed= temp('{species}/dnascent/{sample_id}.fork.tmp'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

fin <- strsplit('{input.bed}', ' ')[[1]]

bed <- list()
for(x in fin) {{
    if(sum(!grepl('#', readLines(x))) == 0) {{
        next
    }}
    feature <- sub('_DNAscent_forkSense.bed', '', basename(x))
    feature <- sub('s$', '', feature)
    xin <- fread(sprintf('grep -v "^#" %s', x), header= FALSE, col.names= c('#chrom', 'start', 'end', 'rname', 'chrom2', 'aln_start', 'aln_end', 'strand'))
    stopifnot(identical(xin$`#chrom`, xin$chrom2))
    xin[, chrom2 := NULL]
    xin[, feature := feature]
    bed[[length(bed) + 1]] <- xin
}}
bed <- rbindlist(bed)

colorder <- c('#chrom', 'start', 'end', 'rname', 'feature', 'strand', 'aln_start', 'aln_end')
if(nrow(bed) == 0) {{
    bed <- setNames(data.table(matrix(nrow= 0, ncol= length(colorder))), colorder)
}} else {{
    bed[, strand := ifelse(strand == 'fwd', '+', ifelse(strand == 'rev', '-', NA))]
    stopifnot(!is.na(bed$strand))
    bed <- bed[order(`#chrom`, start, end)]

    setcolorder(bed, colorder)
}}

write.table(x= bed, file= '{output.bed}', sep= '\t', row.names= FALSE, quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


checkpoint split_detect:
    input:
        detect= '{species}/dnascent/{sample_id}.detect',
    output:
        outdir= temp(directory('{species}/dnascent/split_detect_{sample_id}')),
    shell:
        r"""
        mkdir -p {output.outdir}
        {workflow.basedir}/scripts/split_dnascent.py {input.detect} --prefix {output.outdir}/ -n 10000 -s '.detect'
        """


def aggregate_detect(wc):
    checkpoint_output = checkpoints.split_detect.get(**wc).output.outdir
    return expand('{species}/dnascent/split_detect_{sample_id}/{n}.bedgraph',
        species=wc.species,
        sample_id=wc.sample_id,
        n=glob_wildcards(os.path.join(checkpoint_output, "{n}.detect")).n)


rule cat_detect:
    input:
        bdg=aggregate_detect,
    output:
        bdg= '{species}/dnascent/{sample_id}.bedgraph.gz',
    resources:
        mem='20G'
    shell:
        r"""
        LC_ALL=C sort --parallel 8 -S 5G -T . -m -k1,1 -k2,2n {input.bdg} \
        | awk 'NR == 1 || $0 !~ "^#"' \
        | bgzip -@ 4 > {output.bdg}
        tabix -f -p bed {output.bdg}
        """


rule brdu_stats_for_forks:
    input:
        bed= '{species}/dnascent/{sample_id}.fork.tmp',
        bdg= '{species}/dnascent/{sample_id}.bedgraph.gz',
    output:
        bed= '{species}/dnascent/{sample_id}.fork.bed',
    shell:
        r"""
        zcat {input.bdg} \
        | {workflow.basedir}/scripts/stat_forks.py -c 0.5 -f {input.bed} -b - > {output.bed}
        """


rule brdu_bedgraph:
    input:
        bdg= '{species}/dnascent/{sample_id}.bedgraph.gz',
    output:
        bdg= '{species}/dnascent/{sample_id}.brdu.bedgraph.gz',
    params:
        brdu_cutoff= 0.8
    shell:
        r"""
        zcat {input.bdg} \
        | awk -v cutoff={params.brdu_cutoff} -v OFS='\t' '$0 !~ "^#" {{if($4 > cutoff) {{$4=1}} else {{$4=0}} print $0}}' \
        | bedtools groupby -g 1,2,3 -c 4,4 -o count,sum \
        | awk -v OFS='\t' 'BEGIN{{print "#chrom", "start", "end", "pct_brdu", "depth", "n_brdu"}} 
            {{print $1, $2, $3, sprintf("%.1f", 100 * $5/$4), $4, $5}}' \
        | bgzip -@ 4 > {output.bdg}
        tabix -p bed {output.bdg}
        """

checkpoint split_forkSense:
    input:
        forkSense= '{species}/dnascent/{sample_id}/dnascent.forkSense',
    output:
        outdir= temp(directory('{species}/dnascent/split_forkSense_{sample_id}')),
    shell:
        r"""
        mkdir -p {output.outdir}
        {workflow.basedir}/scripts/split_dnascent.py {input.forkSense} --prefix {output.outdir}/ -n 10000 -s '.forkSense'
        """


def aggregate_forkSense(wc):
    checkpoint_output = checkpoints.split_forkSense.get(**wc).output.outdir
    return expand('{species}/dnascent/split_forkSense_{sample_id}/{n}.forkSense.bedgraph',
        species=wc.species,
        sample_id=wc.sample_id,
        n=glob_wildcards(os.path.join(checkpoint_output, "{n}.forkSense")).n)


rule cat_forkSense:
    input:
        bdg=aggregate_forkSense,
    output:
        bdg= '{species}/dnascent/{sample_id}.forkSense.bedgraph.gz',
    resources:
        mem='20G'
    shell:
        r"""
        LC_ALL=C sort --parallel 8 -S 5G -T . -m -k1,1 -k2,2n {input.bdg} \
        | awk 'NR == 1 || $0 !~ "^#"' \
        | bgzip -@ 4 > {output.bdg}
        tabix -f -p bed {output.bdg}
        """


rule count_reads:
    input:
        bam= '{species}/minimap2/{sample_id}.bam',
    output:
        txt= temp('{species}/dnascent/{sample_id}.detect.tmp'),
    params:
        rlen=rules.dnascent_detect.params.rlen,
        mapq=rules.dnascent_detect.params.mapq,
    shell:
        r"""
        samtools view -u -q {params.mapq} {input.bam} \
        | bamToBed \
        | awk '($3 - $2) >= {params.rlen}' \
        | awk -v OFS='\t' '{{n += 1; len += ($3-$2)}} END {{print n, len}}' > {output.txt}
        """

rule summarise_events:
    input:
        forks= set(expand('{species}/dnascent/{sample_id}.fork.bed', zip, species= ss.species, sample_id= ss.sample_id)),
        reads= set(expand('{species}/dnascent/{sample_id}.detect.tmp', zip, species= ss.species, sample_id= ss.sample_id)),
        ss= os.path.join(workflow.basedir, 'sample_sheet.tsv')
    output:
        tsv= os.path.join(workflow.basedir, 'results/event_summary.tsv'),
    shell:
        r"""
cat <<'EOF' > {rule}.$$.tmp.R
library(data.table)

ss <- fread('{input.ss}')
ff <- strsplit('{input.forks}', ' ')[[1]]
rr <- strsplit('{input.reads}', ' ')[[1]]

forks <- list()
for(x in ff) {{
    txt <- fread(x)
    sample_id <- sub('\\.fork\\.bed$', '', basename(x))
    txt[, sample_id := sample_id]
    txt <- txt[, list(n_events= .N), list(sample_id, feature)]
    forks[[length(forks) + 1]] <- txt
}}
forks <- rbindlist(forks)

reads <- list()
for(x in rr) {{
    txt <- fread(x)
    sample_id <- sub('\\.detect\\.tmp$', '', basename(x))
    txt[, sample_id := sample_id]
    reads[[length(reads) + 1]] <- txt
}}
reads <- rbindlist(reads)
setnames(reads, c('V1', 'V2'), c('n_reads', 'aln_len'))

smry <- merge(forks, reads, by= 'sample_id', all= TRUE)
smry[, events_per_1Mb_aln := 1000000 * n_events / aln_len]

smry <- dcast(data= smry, sample_id + n_reads + aln_len ~ feature, value.var= c('events_per_1Mb_aln', 'n_events'), fill= 0)
smry[, n_events_NA := NULL]
smry[, events_per_1Mb_aln_NA := NULL]
setnames(smry, names(smry), sub('n_events_', '', names(smry)))

idx <- grep('events_per_1Mb_aln', names(smry))

for(i in idx) {{
    old <- names(smry)[i]
    new <- sub('events_per_1Mb_aln_', '', old) 
    new <- sprintf('%s_per_1Mb_aln', new)
    setnames(smry, old, new)
}}

smry <- merge(unique(ss[, list(sample_id, species, sample_type)]), smry, by= 'sample_id')
smry <- smry[order(species, sample_type, sample_id)]

# Convert floating points to string for printing
for(i in grep('_per_1Mb_aln', names(smry))) {{
    smry[[i]] <- sprintf('%.3f', smry[[i]])
}}
setnames(smry, 'aln_len', 'aln_len_Mb')
smry[, aln_len_Mb := sprintf('%.1f', aln_len_Mb/1000000)]

write.table(smry, '{output.tsv}', sep= '\t', row.names= FALSE, quote= FALSE)

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """


