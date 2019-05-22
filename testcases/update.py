#!/usr/bin/env python

from __future__ import print_function

import os
import tempfile
import glob
import argparse
import yaml

# Try to load pyfabm.
# If available, we'll use it to automatically add comments to fabm.yaml.
try:
    import pyfabm.complete_yaml
except ImportError:
    print('Unable to load pyfabm. Automatic addition of comments to fabm.yaml will be disabled. To fix this, build and install pyfabm by running cmake+make for source directory <FABMDIR>/src/drivers/python.')
    pyfabm = None

# ------------------------------------------
# Hook into PyYAML to make it preserve the order of dictionary elements (i.e., modules in fabm.yaml).
try:
    import collections
    def dict_representer(dumper, data):
        return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.items())
    def dict_constructor(loader, node):
        return collections.OrderedDict(loader.construct_pairs(node))
    yaml.add_representer(collections.OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
except ImportError:
    pass
# ------------------------------------------

def enumerateModels(fabm_yaml):
    for name,data in fabm_yaml.get('instances',{}).items():
        if isinstance(data,dict):
            yield name,data.get('model',None),data.get('parameters',{}),data.get('coupling',{})

def update(path):
    print('Updating %s...' % path)
    with open(path,'rU') as f:
        fabm_yaml = yaml.load(f)

    # Add ersem/benthic_erosion module if required
    needs_erosion = False
    has_erosion = False
    for name,model,parameters,coupling in enumerateModels(fabm_yaml):
        if model=='ersem/benthic_column_particulate_matter':
            needs_erosion = needs_erosion or parameters.get('resuspension', False)
            parameters.pop('er', None)
            parameters.pop('vel_crit', None)
        elif model=='ersem/benthic_erosion':
            has_erosion = True
    if needs_erosion and not has_erosion:
        print('  Adding benthic erosion module.')
        fabm_yaml['instances']['erosion'] = {'model':'ersem/benthic_erosion'}

    # Update ersem/benthic_bacteria to deal with variable number of food sources
    for name,model,parameters,coupling in enumerateModels(fabm_yaml):
        if model!='ersem/benthic_bacteria' or 'nfood' in parameters: continue
        print('  Updating benthic bacteria %s.' % name)
        nfood = 0
        if parameters.get('suQ1',0)!=0:
            nfood += 1
            parameters['su%i' % nfood] = parameters.pop('suQ1',None)
            coupling['food%i' % nfood] = coupling['Q1']
        if parameters.get('suQ6s',0)!=0 or parameters.get('suQ6f',0)!=0:
            nfood += 1
            if 'suQ6s' in parameters: parameters['su%i'    % nfood] = parameters.pop('suQ6s',None)
            if 'suQ6f' in parameters: parameters['suf%i'   % nfood] = parameters.pop('suQ6f',None)
            if 'pue6'  in parameters: parameters['pue%i'   % nfood] = parameters.pop('pue6',None)
            if 'puinc' in parameters: parameters['puinc%i' % nfood] = parameters.pop('puinc',None)
            coupling['food%i' % nfood] = coupling['Q6']
        if parameters.get('suQ7',0)!=0:
            nfood += 1
            parameters['su%i' % nfood] = parameters.pop('suQ7')
            if 'pue7' in parameters: parameters['pue%i' % nfood] = parameters.pop('pue7',None)
            coupling['food%i' % nfood] = coupling['Q7']
        parameters['nfood'] = nfood
        parameters.pop('suQ1',None)
        parameters.pop('suQ6s',None)
        parameters.pop('suQ6f',None)
        parameters.pop('suQ7',None)
        coupling.pop('Q7',None)

    if pyfabm is not None:
        # Write YAML to temporary file.
        with tempfile.NamedTemporaryFile('w',delete=False) as f:
            yaml.dump(fabm_yaml, f, default_flow_style=False, indent=2)

        print('  Cleaning YAML...')
        pyfabm.complete_yaml.processFile(f.name,path)

        os.remove(f.name) 
    else:
        with open(path,'w') as f:
            yaml.dump(fabm_yaml, f, default_flow_style=False, indent=2)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='This script updates a fabm.yaml ERSEM configuration to make it compatible with the latest version of the ERSEM code.')
    parser.add_argument('path',nargs='+',help='Path to your fabm.yaml file that will be updated.',default=[])
    args = parser.parse_args()
    for path in args.path:
        for p in glob.glob(path): update(p)