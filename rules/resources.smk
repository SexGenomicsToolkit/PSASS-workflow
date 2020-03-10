'''
Functions and variables shared between several snakefiles.
'''

import collections
import logging
import os
import re
import yaml


def flatten(l):
    '''
    '''
    for el in l:
        if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


def allocate_resources():
    '''
    '''
    parsed_resources = collections.defaultdict(dict)
    cfg_resources = yaml.safe_load(open(config['resources_info']))
    threads_default = config['presets']['threads'][cfg_resources['default']['threads']]
    mem_mb_default = config['presets']['mem_mb'][cfg_resources['default']['mem_mb']]
    runtime_default = config['presets']['runtime'][cfg_resources['default']['runtime']]
    for rule, cfg in cfg_resources.items():
        if rule == 'default':
            continue
        try:
            parsed_resources[rule]['threads'] = config['presets']['threads'][cfg['threads']]
        except KeyError:
            logging.warning(f'Could not find threads resources preset <{cfg["threads"]}> for rule <{k}>. Using default ({threads_default})')
            parsed_resources[rule]['threads'] = threads_default
        try:
            parsed_resources[rule]['mem_mb'] = config['presets']['mem_mb'][cfg['mem_mb']]
        except KeyError:
            logging.warning(f'Could not find mem_mb resources preset <{cfg["mem_mb"]}> for rule <{k}>. Using default ({mem_mb_default})')
            parsed_resources[rule]['mem_mb'] = mem_mb_default
        try:
            parsed_resources[rule]['runtime'] = config['presets']['runtime'][cfg['runtime']]
        except KeyError:
            logging.warning(f'Could not find runtime resources preset <{cfg["runtime"]}> for rule <{k}>. Using default ({runtime_default})')
            parsed_resources[rule]['runtime'] = runtime_default
    if config['resources']:
        for rule, specs in config['resources'].items():
            try:
                for spec, value in specs.items():
                    parsed_resources[rule][spec] = value
            except KeyError:
                logging.warning(f'Invalid resource <{spec}> or rule name <{rule}> in resources section of config file')
    config['resources'] = parsed_resources
    config['resources']['default'] = cfg_resources['default']


def get_threads(rule):
    '''
    '''
    try:
        threads = config['resources'][rule]['threads']
    except KeyError:
        threads = config['resources']['default']['threads']
    return threads


def get_mem(rule, attempt):
    '''
    Memory increased 1.5x per attempt
    '''
    try:
        mem_mb = config['resources'][rule]['mem_mb']
    except KeyError:
        mem_mb = config['presets']['mem_mb'][config['resources']['default']['mem_mb']]
    if isinstance(mem_mb, (int, float)):
        mem_mb = int(mem_mb * (1 + (attempt - 1) / 2))
    elif mem_mb.isdigit():
        mem_mb = int(int(mem_mb) * (1 + (attempt - 1) / 2))
    elif mem_mb[-1] in ('G', 'M', 'K'):  # Careful, this cannot be used in resources (not an int)
        tmp = float(mem_mb[:-1]) * (1 + (attempt - 1) / 2)
        mem_mb = f'{tmp}{mem_mb[-1]}'
    return mem_mb


def get_runtime(rule, attempt):
    '''
    '''
    try:
        runtime = config['resources'][rule]['runtime']
    except KeyError:
        runtime = config['presets']['runtime'][config['resources']['default']['runtime']]
    if isinstance(runtime, int):
        time = runtime
    else:
        try:
            d, h, m, s = (int(f) for f in re.split(':|-', runtime))
            time = ((((d * 24) + h) * 60) + m) * 60 + s
        except ValueError:
            logging.warning(f'Invalid runtime format for rule <{rule}>: <{runtime}>')
            time = 3600
    time = int(time * (1 + (attempt - 1) / 2))
    return time


# Parse all resources specification sources and populate the config dictionary with resources for each rule
allocate_resources()
