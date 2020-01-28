
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
        'benchmarks/bwa/index.tsv'
    log:
        'logs/bwa_index.tsv'
    conda:
        '../envs/psass.yaml'
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
        'benchmarks/align_{pool}_{lane}.tsv'
    log:
        'logs/align_{pool}_{lane}.txt'
    conda:
        '../envs/psass.yaml'
    threads:
        config['resources']['alignment']['threads']
    resources:
        memory = lambda wildcards, attempt: config['resources']['alignment']['memory'] * attempt
    params:
        runtime = config['resources']['alignment']['runtime']
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
        'benchmarks/sort_{pool}_{lane}.tsv'
    log:
        'logs/sort_{pool}_{lane}.txt'
    conda:
        '../envs/psass.yaml'
    threads:
        config['resources']['sort']['threads']
    resources:
        memory = lambda wildcards, attempt: config['resources']['sort']['memory'] * attempt
    params:
        runtime = config['resources']['sort']['runtime']
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
        'benchmarks/merge_{pool}.tsv'
    log:
        'logs/merge_{pool}.txt'
    conda:
        '../envs/psass.yaml'
    threads:
        config['resources']['merge']['threads']
    resources:
        memory = lambda wildcards, attempt: config['resources']['merge']['memory'] * attempt
    params:
        runtime = config['resources']['merge']['runtime']
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
        'benchmarks/rmdup_{pool}.tsv'
    log:
        'logs/rmdup_{pool}.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['rmdup']['memory'] * attempt
    params:
        runtime = config['resources']['rmdup']['runtime']
    shell:
        'samtools rmdup {input} {output} 2> {log}'


rule samtools_mpileup:
    '''
    '''
    input:
        pool1 = 'output/{pool1}.no_duplicates.cram',
        pool2 = 'output/{pool2}.no_duplicates.cram',
        reference = config['assembly']
    output:
        'output/{pool1}_{pool2}_nucleotides.tsv'
    benchmark:
        'benchmarks/mpileup_{pool1}_{pool2}.tsv'
    log:
        'logs/mpileup_{pool1}_{pool2}.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['mpileup']['memory'] * attempt
    params:
        runtime = config['resources']['mpileup']['runtime'],
        min_quality = config['mpileup']['min_quality']
    shell:
        'samtools mpileup -f {input.reference} -Q {params.min_quality} -aa '
        '{input.pool1} {input.pool2} 2> {log} | '
        'psass convert - > {output}'
