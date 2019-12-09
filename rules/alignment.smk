
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
        'logs/bwa/index.tsv'
    conda:
        '../envs/psass.yaml'
    shell:
        'bwa index {input} 2> {log}'


def get_reads(wildcards):
    '''
    '''
    reads = config['reads'][wildcards.sample][wildcards.lane]
    return [reads['R1'], reads['R2']]


rule bwa_mem:
    '''
    '''
    input:
        index_files = rules.bwa_index.output,
        assembly_file = config['assembly'],
        reads_files = get_reads
    output:
        temp('output/{sample}_{lane}.cram')
    benchmark:
        'benchmarks/align_{sample}_{lane}.tsv'
    log:
        'logs/align_{sample}_{lane}.txt'
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
        temp('output/{sample}_{lane}.sorted.cram')
    benchmark:
        'benchmarks/sort_{sample}_{lane}.tsv'
    log:
        'logs/sort_{sample}_{lane}.txt'
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


def merge_sample_input(wildcards):
    '''
    '''
    sample = wildcards.sample
    lanes = [k for k in config['reads'][sample]]
    return expand('output/{sample}_{lane}.sorted.cram', sample=sample, lane=lanes)


rule samtools_merge:
    '''
    '''
    input:
        merge_sample_input
    output:
        'output/{sample}.cram'
    benchmark:
        'benchmarks/merge_{sample}.tsv'
    log:
        'logs/merge_{sample}.txt'
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
        'output/{sample}.no_duplicates.cram'
    benchmark:
        'benchmarks/rmdup_{sample}.tsv'
    log:
        'logs/rmdup_{sample}.txt'
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
        sample1 = 'output/{sample1}.no_duplicates.cram',
        sample2 = 'output/{sample2}.no_duplicates.cram',
        reference = config['assembly']
    output:
        'output/{sample1}_{sample2}.pileup'
    benchmark:
        'benchmarks/mpileup_{sample1}_{sample2}.tsv'
    log:
        'logs/mpileup_{sample1}_{sample2}.txt'
    conda:
        '../envs/psass.yaml'
    threads:
        config['resources']['mpileup']['threads']
    resources:
        memory = lambda wildcards, attempt: config['resources']['rmdup']['memory'] * attempt
    params:
        runtime = config['resources']['rmdup']['runtime']
    shell:
        'bcftools mpileup -d 500 -f {input.reference} -Q 0 --threads {threads} '
        '-o {output} {input.sample1} {input.sample2} 2> {log}'
