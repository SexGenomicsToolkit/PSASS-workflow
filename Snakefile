import itertools

wildcard_constraints:
    pool = '[^.]+'

configfile: 'config.yaml'

include: 'rules/resources.smk'
include: 'rules/processing.smk'
include: 'rules/analyses.smk'


def all_pairs(wildcards):
    '''
    '''
    pools = [s for s in config['reads']]
    input_files = []
    for pool1, pool2 in itertools.combinations(pools, 2):
        input_files.append(rules.psass_pileup.output[0].format(pool1=pool1, pool2=pool2))
        input_files += expand(rules.psass.output.window, pool1=pool1, pool2=pool2, preset=config['psass'])
    return input_files


rule all:
    input:
        all_pairs
