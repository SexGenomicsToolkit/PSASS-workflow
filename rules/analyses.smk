
def generate_psass_param_string(wildcards):
    '''
    '''
    params = ' '.join(f'--{k}' if isinstance(v, bool) else f'--{k} {v}' for
                      k, v in config['psass'][wildcards.preset].items())
    return params


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
    threads: get_threads('psass')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('psass', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('psass', attempt)
    params:
        psass_settings = generate_psass_param_string
    shell:
        'psass analyze {input} {output.window} '
        '--snp-file {output.snps} --fst-file {output.fst} '
        '{params.psass_settings} 2> {log}'


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
    threads: get_threads('circos_plot')
    resources:
        mem_mb = lambda wildcards, attempt: get_mem('circos_plot', attempt),
        runtime_s = lambda wildcards, attempt: get_runtime('circos_plot', attempt)
    params:
        prefix = 'output/psass/{pool1}_{pool2}/{preset}'
    shell:
        'echo nothing'
