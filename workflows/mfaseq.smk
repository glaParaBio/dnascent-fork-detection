rule bwa_index:
    input:
        fa= 'ref/{species}_Genome.fasta',
    output:
        bwt= 'ref/{species}_Genome.fasta.bwt',
    shell:
        r"""
        bwa index {input.fa}
        """

def get_mfa_fastq(wc):
    fq = [mfa_ss[mfa_ss.sample_id == wc.mfa_id].R1.iloc[0]]
    fq2 = mfa_ss[mfa_ss.sample_id == wc.mfa_id].R2.iloc[0]
    if not pandas.isna(fq2):
        fq.append(fq2)
    if config['remote']:
        for i in range(len(fq)):
            if fq[i] in remote_fastq:
                fq[i] = FTP.remote(fq[i], immediate_close= False, keep_local= False)
    return fq

rule bwa_align:
    input:
        fq= get_mfa_fastq,
        bwt= 'ref/{species}_Genome.fasta.bwt',
        fa= 'ref/{species}_Genome.fasta',
    output:
        bam= 'bwa/{mfa_id}.bam',
        md= 'bwa/{mfa_id}.md',
    shell:
        r"""
        bwa mem -R '@RG\tID:{wildcards.mfa_id}\tSM:{wildcards.mfa_id}\tLB:{wildcards.mfa_id}\tPL:ILLUMINA\tPU:NA' \
            -t 16 {input.fa} {input.fq} \
        | samtools fixmate -m -@ 4 - - \
        | samtools sort -@ 4 \
        | samtools markdup -l 300 -f {output.md} -@ 4 - {output.bam}
        samtools index -@ 4 {output.bam}
        """

rule bigwig:
    input:
        bam= 'bwa/{mfa_id}.bam',
    output:
        bw= 'bigwig/{mfa_id}.bw',
    shell:
        r"""
        bamCoverage -b {input.bam} -o {output.bw} \
            --minMappingQuality 30 \
            --binSize 25 \
            --normalizeUsing BPM \
            --numberOfProcessors 16
        """
