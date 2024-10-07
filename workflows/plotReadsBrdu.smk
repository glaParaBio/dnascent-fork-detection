rule plotReadsBrdu:
    input:
        expand('TriTrypDB-45_LmajorFriedlin/dnascent/plotReadsBrdu/any/{sample_id}/{chrom}.pdf',
            sample_id=ss[(ss.use == 'yes') & (ss.species == 'TriTrypDB-45_LmajorFriedlin')].sample_id, 
               chrom=LmjF_CHROMS),

        expand('TriTrypDB-45_LmajorFriedlin/dnascent/plotReadsBrdu/only_with_fork/{sample_id}/{chrom}.pdf',
            sample_id=ss[(ss.use == 'yes') & (ss.species == 'TriTrypDB-45_LmajorFriedlin')].sample_id, 
               chrom=LmjF_CHROMS),


rule extend_mfa_peaks:
    input:
        mfa_peaks= config['mfaseq_peaks'],
        fai= 'ref/TriTrypDB-45_LmajorFriedlin_Genome.fasta.fai',
    output:
        bed= temp('mfaseq_slop.bed'),
    shell:
        r"""
        slopBed -b 50000 -i {input.mfa_peaks} -g {input.fai} > {output.bed}
        """


rule brdu_intersecting_mfa_peaks:
    input:
        mfa_peaks= 'mfaseq_slop.bed',
        bdg= '{species}/dnascent/{sample_id}.bedgraph.gz',
    output:
        bdg= temp('{species}/dnascent/plotReadsBrdu/{sample_id}/{chrom}.bedgraph'),
    shell:
        r"""
        tabix -h {input.bdg} {wildcards.chrom} \
        | awk -v OFS='\t' '$4 > 0.2 {{print $1, $8, $9, $4, $5, $6, $7, $2, $3}}' \
        | intersectBed -a - -b {input.mfa_peaks} -header -wa > {output.bdg}
        """


rule plot_brdu_mfa_peaks:
    input:
        bed= '{species}/dnascent/{sample_id}.fork.bed',
        bdg= '{species}/dnascent/plotReadsBrdu/{sample_id}/{chrom}.bedgraph',
        ref_points= config['mfaseq_peaks'],
    output:
        plot= '{species}/dnascent/plotReadsBrdu/any/{sample_id}/{chrom}.pdf',
    shell:
        r"""
        {workflow.basedir}/scripts/plotReadsBrdu.r -i {input.bdg} -o {output.plot} \
            -r {input.ref_points} \
            -c {wildcards.chrom} \
            -n {wildcards.sample_id} \
            --downsample 100 \
            --annotate {input.bed} \
            --min-aln-len 2000
        """


rule plot_brdu_mfa_peaks_with_fork:
    input:
        bed= '{species}/dnascent/{sample_id}.fork.bed',
        bdg= '{species}/dnascent/plotReadsBrdu/{sample_id}/{chrom}.bedgraph',
        ref_points= config['mfaseq_peaks'],
        forks= '{species}/dnascent/{sample_id}.fork.bed',
    output:
        plot= '{species}/dnascent/plotReadsBrdu/only_with_fork/{sample_id}/{chrom}.pdf',
    shell:
        r"""
        awk '$1 == "{wildcards.chrom}" || $1 == "#chrom"' {input.forks} \
        | cut -f 4 \
        | sort \
        | uniq \
        | grep -w -F -f - {input.bdg} \
        | {workflow.basedir}/scripts/plotReadsBrdu.r -i - -o {output.plot} \
            -r {input.ref_points} \
            -c {wildcards.chrom} \
            -n {wildcards.sample_id} \
            --downsample 100 \
            --annotate {input.bed} \
            --min-aln-len 2000
        """


