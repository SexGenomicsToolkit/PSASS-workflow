

rule psass:
    '''
    '''
    input:
        rules.samtools_mpileup.output
    output:
        'output/psass_{sample1}_{sample2}/{preset}_depth.tsv',
        'output/psass_{sample1}_{sample2}/{preset}_fst_position.tsv',
        'output/psass_{sample1}_{sample2}/{preset}_fst_window.tsv',
        'output/psass_{sample1}_{sample2}/{preset}_snps_position.tsv',
        'output/psass_{sample1}_{sample2}/{preset}_snps_window.tsv'
    benchmark:
        'benchmarks/psass_{sample1}_{sample2}/{preset}.tsv'
    log:
        'logs/psass_{sample1}_{sample2}/{preset}.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['psass']['memory'] * attempt
    params:
        runtime = config['resources']['psass']['runtime'],
        psass = ' '.join(f'{k} {v}' for k, v in config['psass']['{preset}']),
        prefix = 'output/psass/{sample1}_{sample2}/{preset}'
    shell:
        'psass analyze --input-file {input} --output-prefix {params.prefix} '
        '{params.psass} 2> {log}'


rule circos_plot:
    '''
    '''
    input:
        rules.psass.output
    output:
        'output/psass_{sample1}_{sample2}/{preset}.png'
    benchmark:
        'benchmarks/psass_{sample1}_{sample2}/{preset}_plot.tsv'
    log:
        'logs/psass_{sample1}_{sample2}/{preset}_plot.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['psass-vis']['memory'] * attempt
    params:
        runtime = config['resources']['psass-vis']['runtime'],
        prefix = 'output/psass/{sample1}_{sample2}/{preset}'
    script:

