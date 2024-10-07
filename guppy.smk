FASTQ_DIR = 'fastq'

run_dirs = config['run_dirs'].split(' ')

fastq = [os.path.join(x, f'{FASTQ_DIR}/sequencing_summary.txt') for x in run_dirs]

rule all:
    input:
        fastq,   


rule guppy_basecaller:
    input:
        run_dir='{run_dir}',
    output:
        '{run_dir}/%s/sequencing_summary.txt' % FASTQ_DIR,
    params:
        guppy_dir=os.path.abspath(config['guppy_dir']),
        guppy_config=config['guppy_config'],
    shell:
        r"""
        n_gpu=`nvidia-smi -L | wc -l`
        device='cuda:'
        for ((i=0; i<$n_gpu; i++))
        do
            device=${{device}}"$i,"
        done
        device=`echo $device | sed 's/,$//'`

        cd {input.run_dir} 

        {params.guppy_dir}/bin/guppy_basecaller \
            -c {params.guppy_dir}/data/{params.guppy_config} \
            -i . \
            -s ./fastq \
            --num_callers 50 \
            --recursive \
            --compress_fastq \
            --device "$device" \
            --disable_qscore_filtering
        """
