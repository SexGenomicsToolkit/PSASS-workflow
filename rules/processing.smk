
rule bwa_index:
    '''
    '''
    input:
        config['assembly']
    output:
        f'{config["assembly"]}.amb',
        f'{config["assembly"]}.ann',
        f'{config["assembly"]}.bwt',
        f'{config["assembly"]}.pac',
        f'{config["assembly"]}.sa'
    benchmark:
        'benchmarks/bwa_index.tsv'
    log:
        'logs/bwa_index.tsv'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('bwa_index')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('bwa_index', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('bwa_index', attempt)
    shell:
        'bwa index {input} 2> {log}'


def get_reads(wildcards):
    '''
    '''
    reads = config['reads'][wildcards.pool][wildcards.lane]
    if 'R1' in reads:
        return [reads['R1'], reads['R2']]
    else:
        return reads


rule bwa_mem:
    '''
    '''
    input:
        index_files = rules.bwa_index.output,
        assembly_file = config['assembly'],
        reads_files = get_reads
    output:
        temp('output/{pool}_{lane}.cram')
    benchmark:
        'benchmarks/bwa_mem_{pool}_{lane}.tsv'
    log:
        'logs/bwa_mem_{pool}_{lane}.txt'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('bwa_mem')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('bwa_mem', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('bwa_mem', attempt)
    shell:
        'bwa mem -t {threads} {input.assembly_file} {input.reads_files} 2> {log} | '
        'samtools view -C -h -T {input.assembly_file} -o {output} 2>> {log}'


rule samtools_sort:
    '''
    '''
    input:
        rules.bwa_mem.output
    output:
        temp('output/{pool}_{lane}.sorted.cram')
    benchmark:
        'benchmarks/samtools_sort_{pool}_{lane}.tsv'
    log:
        'logs/samtools_sort_{pool}_{lane}.txt'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('samtools_sort')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('samtools_sort', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('samtools_sort', attempt)
    shell:
        'samtools sort -@ {threads} -o {output} {input} 2> {log}'


def merge_pool_input(wildcards):
    '''
    '''
    pool = wildcards.pool
    lanes = [k for k in config['reads'][pool]]
    return expand('output/{pool}_{lane}.sorted.cram', pool=pool, lane=lanes)


rule samtools_merge:
    '''
    '''
    input:
        merge_pool_input
    output:
        'output/{pool}.cram'
    benchmark:
        'benchmarks/samtools_merge_{pool}.tsv'
    log:
        'logs/samtools_merge_{pool}.txt'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('samtools_merge')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('samtools_merge', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('samtools_merge', attempt)
    shell:
        'samtools merge -r -@ {threads} {output} {input} 2> {log}'


rule samtools_rmdup:
    '''
    '''
    input:
        rules.samtools_merge.output
    output:
        'output/{pool}.no_duplicates.cram'
    benchmark:
        'benchmarks/samtools_rmdup_{pool}.tsv'
    log:
        'logs/samtools_rmdup_{pool}.txt'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('samtools_rmdup')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('samtools_rmdup', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('samtools_rmdup', attempt)
    shell:
        'samtools rmdup {input} {output} 2> {log}'


rule psass_pileup:
    '''
    '''
    input:
        pool1 = 'output/{pool1}.no_duplicates.cram',
        pool2 = 'output/{pool2}.no_duplicates.cram',
        reference = config['assembly']
    output:
        'output/{pool1}_{pool2}.sync'
    benchmark:
        'benchmarks/psass_pileup_{pool1}_{pool2}.tsv'
    log:
        'logs/psass_pileup_{pool1}_{pool2}.txt'
    conda:
        '../envs/psass.yaml'
    threads: get_threads('psass_pileup')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('psass_pileup', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('psass_pileup', attempt)
    params:
        min_quality = config['psass_pileup']['min_quality']
    shell:
        'psass pileup -r {input.reference} -q {params.min_quality} '
        '-o {output} {input.pool1} {input.pool2} 2> {log}'
