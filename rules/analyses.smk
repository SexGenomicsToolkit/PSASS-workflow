

rule psass:
    '''
    '''
    input:
        rules.psass_pileup.output
    output:
        window = 'output/psass_{pool1}_{pool2}/{preset}_window.tsv',
        fst = 'output/psass_{pool1}_{pool2}/{preset}_fst.tsv',
        snps = 'output/psass_{pool1}_{pool2}/{preset}_snps.tsv',
    benchmark:
        'benchmarks/psass_{pool1}_{pool2}/{preset}.tsv'
    log:
        'logs/psass_{pool1}_{pool2}/{preset}.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['psass']['memory'] * attempt
    params:
        runtime = config['resources']['psass']['runtime'],
        psass = lambda wildcards: ' '.join(f'--{k} {v}' for k, v in config['psass'][wildcards.preset].items())
    shell:
        'psass analyze {input} {output.window} '
        '--snps-file {output.snps} --fst-file {output.fst} '
        '{params.psass} 2> {log}'


rule circos_plot:
    '''
    '''
    input:
        rules.psass.output
    output:
        'output/psass_{pool1}_{pool2}/{preset}.png'
    benchmark:
        'benchmarks/psass_{pool1}_{pool2}/{preset}_plot.tsv'
    log:
        'logs/psass_{pool1}_{pool2}/{preset}_plot.txt'
    conda:
        '../envs/psass.yaml'
    resources:
        memory = lambda wildcards, attempt: config['resources']['psass-vis']['memory'] * attempt
    params:
        runtime = config['resources']['psass-vis']['runtime'],
        prefix = 'output/psass/{pool1}_{pool2}/{preset}'
    shell:
        'echo nothing'

