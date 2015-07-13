#!/usr/bin/env python

# This script takes a legacy ERSEM configuration and converts it into a
# YAML-based configuration file (typically "fabm.yaml") used by FABM-ERSEM.
#
# The legacy ERSEM configuration is described by namelist files in two directories:
# BioParams and include. The paths to these directories must be provided as 1st and
# 2nd commend line arguments. The mapping between old and new parameters is defined
# in fabm-ersem-template.yaml.
#
# Output of the script is a new "fabm.yaml" file (name is controllable with the 3rd command line argument),
# with parameter values taken from the legacy namelist files.
#
# It is recommended you run with the -c/--clean switch, which will create a standardized, self-documenting
# YAML file. That does require FABM's python driver to have been built AND installed.

import sys,os,optparse,tempfile,fnmatch,re
import xml.etree.ElementTree
import yaml

template = os.path.join(os.path.dirname(__file__),'fabm-ersem-template.yaml')

parser = optparse.OptionParser(usage='usage: %prog [options] BioParams_dir include_dir [output_file]')
parser.add_option('--docdyn', action='store_true',help='enable DOC dynamics (equivalent to legacy ERSEM with DOCDYN defined)',default=False)
parser.add_option('--iop', action='store_true',help='enable IOP based light model (equivalent to legacy ERSEM with IOPMODEL defined)',default=False)
parser.add_option('-c', '--clean',    action='store_true', help='clean-up YAML file (sort parameters, add descriptions)',default=False)
parser.add_option('-v', '--verbose',  action='store_true', help='detailed output',default=False)
parser.add_option('-t', '--template', type='string',help='path to template file (default: %s)' % template,default=template)
parser.add_option('--subtract_background', action='store_true',help='subtract the background value from the initial state')
parser.add_option('--ibenxin', type='int',help='force benthic configuration (0: none, 1: benthic returns, 2: Oldenburg benthos, default = autodetect)',default=None)
options, args = parser.parse_args()

# Check arguments
if len(args)<2:
    print 'Two arguments must be provided:'
    print ' - path to BioParams directory'
    print ' - path to include directory'
    print 'A third argument (file path to save the FABM-ERSEM configuration to) is optional; it defaults to fabm.yaml'
    sys.exit(2)

outfile = 'fabm.yaml'

# Get arguments
bioparamspath,includepath = args[:2]
if len(args)>2: outfile = args[2]

if not os.path.isdir(bioparamspath):
    print 'First argument (BioParams) must be the path to a directory.'
    sys.exit(2)
if not os.path.isdir(includepath):
    print 'Second argument (include) must be the path to a directory.'
    sys.exit(2)

ersemdict = {}
ersemsources = {}

def processNamelistsInDirectory(dirpath,dirlabel):
    for name in os.listdir(dirpath):
        if not name.lower().endswith('.nml'): continue
        path = os.path.join(dirpath,name)
        with open(path,'rU') as f:
            for l in f:
                l = l.split('!',1)[0].strip()
                if l and '=' in l:
                    key,value = map(str.strip,map(str.lower,l.split('=',1)))
                    if key in ersemdict: print 'WARNING: duplicate entries for %s found in original namelists.' % key
                    ersemdict[key] = value
                    ersemsources[key] = '%s/%s' % (dirlabel,name)

processNamelistsInDirectory(bioparamspath,'BioParams')
processNamelistsInDirectory(includepath,'include')

ersemdict_used = set()

# ------------------------------------------
# Hook into PyYAML to make it preserve the order of dictionary elements.
try:
    import collections
    def dict_representer(dumper, data):                                                            
        return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())                                                                                         
    def dict_constructor(loader, node):                                                            
        return collections.OrderedDict(loader.construct_pairs(node))                               
    yaml.add_representer(collections.OrderedDict, dict_representer)                                
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
except ImportError:
    pass
# ------------------------------------------

# Read YAML template
print 'Loading template %s...' % options.template
f = open(options.template,'rU')
root = yaml.load(f)
f.close()

def setValue(path,value):
    comp = path.split('/')
    d = root
    for c in comp[:-1]:
        if c not in d: d[c] = {}
        d = d[c]
    d[comp[-1]] = value
    if options.verbose: print 'setting %s to %s' % (path,value)

def deleteValue(path,d=root,prefix='',exp=None):
    if exp is None: exp = re.compile(fnmatch.translate(path))
    n = 0
    for name in d.keys():
        value = d[name]
        if exp.match(prefix+name):
            del d[name]
            if options.verbose: print 'deleting %s' % (prefix+name)
            n += 1
        elif isinstance(value,dict):
            n += deleteValue(path,d=value,prefix=prefix+name+'/',exp=exp)
    assert d is not root or n>0,'No matching keys matching "%s" found.' % path
    return n

def replaceValue(old,new,d=root,path=''):
    for name in d.keys():
        value = d[name]
        if isinstance(value,dict):
            replaceValue(old,new,d=value,path=path+'/'+name)
        elif value==old:
            if options.verbose: print 'Setting %s to %s' % (path+'/'+name,new)
            d[name] = new

# -------------------------------------------------
# Configure ERSEM
# -------------------------------------------------

if options.ibenxin is None: options.ibenxin = int(ersemdict['ibenxin'])
if options.ibenxin==0:
    print 'Configuration without benthos is currently not supported.'
    sys.exit(1)
elif options.ibenxin==1:
    # Benthic return: remove benthic modules from YAML
    deleteValue('instances/Q17')     # buried matter
    deleteValue('instances/K?')      # benthic nutrients
    deleteValue('instances/G?')      # benthic dissolved gases
    deleteValue('instances/H*')      # benthic bacteria
    deleteValue('instances/Y*')      # benthic fauna
    deleteValue('instances/ben_col') # benthic column setup
    deleteValue('instances/Q6s_*')   # Silicate remineralization
    deleteValue('instances/ben_nit')
    deleteValue('instances/Q7/coupling/burial_target')
    deleteValue('instances/Q?/initialization/pen_depth_?')
    deleteValue('instances/Q7/parameters/burial')
    setValue('instances/Q6/model','ersem/benthic_base')
    setValue('instances/Q7/model','ersem/benthic_base')
    replaceValue('Q6/surface','Q6')
    replaceValue('Q7/surface','Q7')
elif options.ibenxin==2:
    # Full benthos: disable benthic returns
    deleteValue('instances/Q?/parameters/remin')
    deleteValue('instances/Q?/parameters/pN3')
    deleteValue('instances/Q?/coupling/O3c')
    deleteValue('instances/Q?/coupling/N??')
else:
    print 'ERROR: unknown value %i used for ibenxin.' % options.ibenxin
    sys.exit(2)

if options.docdyn:
    # replicate legacy ERSEM with DOCDYN defined
    setValue('instances/B1/model','ersem/bacteria_docdyn')
    deleteValue('instances/B1/parameters/redfield')
    deleteValue('instances/B1/parameters/rR2R1')
    deleteValue('instances/B1/parameters/puRP?')
    deleteValue('instances/B1/parameters/R1R2')
    setValue('instances/R1/parameters/c0',0.0034)
    setValue('instances/R2/parameters/c0',0.0033)
    setValue('instances/P1/parameters/docdyn',True)
    setValue('instances/P2/parameters/docdyn',True)
    setValue('instances/P3/parameters/docdyn',True)
    setValue('instances/P4/parameters/docdyn',True)
    deleteValue('instances/P?/parameters/R1R2')
else:
    # replicate legacy ERSEM without DOCDYN defined
    deleteValue('instances/R3')
    deleteValue('instances/B1/coupling/R3c')
    deleteValue('instances/B1/parameters/sRP?R1')
    deleteValue('instances/B1/parameters/rR?')
    deleteValue('instances/B1/parameters/frR3')

if options.iop:
    # replicate legacy ERSEM with IOPMODEL defined
    setValue('instances/light/model','ersem/light_iop')
    deleteValue('instances/light/parameters/EPS0r')
    deleteValue('instances/light/parameters/EPSESS')
    deleteValue('instances/P?/parameters/EPS')
    deleteValue('instances/R?/parameters/EPS')
else:
    # replicate legacy ERSEM without IOPMODEL defined
    deleteValue('instances/zenithAngle')
    deleteValue('instances/light/parameters/a0w')
    deleteValue('instances/light/parameters/b0w')
    deleteValue('instances/P?/parameters/iopABS')
    deleteValue('instances/R?/parameters/iopABS')
    deleteValue('instances/P?/parameters/iopBBS')
    deleteValue('instances/R?/parameters/iopBBS')

print 'NOTE: legacy setup has pCO2 = %s. In FABM this is set through forcing.' % float(ersemdict['pco2a3'])
if 'essxr' in ersemdict.keys(): print 'NOTE: legacy setup has silt = %s (for default light model). In FABM this is set through forcing.' % float(ersemdict['essxr'])
if 'abessxr' in ersemdict.keys(): print 'NOTE: legacy setup has silt absorption = %s (for IOP based light model). In FABM this is set through forcing.' % float(ersemdict['abessxr'])

# Filter out empty dictionaries - PyYAML dumps these as {}, which FABM does not understand.
def removeEmptyDictionaries(d=root):
    for name in d.keys():
        value = d[name]
        if isinstance(value,dict):
            removeEmptyDictionaries(d=value)
            if not value: del d[name]
removeEmptyDictionaries()

def processDictionary(d,depth=0,path=''):
    prefix = ' '*(depth*2)
    for name in sorted(d.keys()):
        value = d[name]

        if isinstance(value,dict):
            if options.verbose: print '%s%s' % (prefix,name)
            processDictionary(value,depth+1,path=path+'/'+name)
            continue
        elif not isinstance(value,basestring) or not value.startswith('$'):
            continue

        ersemvar_lower = value[1:].lower()
        if ersemvar_lower not in ersemdict:
            print '%s%s/%s: ERROR: legacy ERSEM variable %s not found.' % (prefix,path,name,value)
            sys.exit(1)
        ersemdict_used.add(ersemvar_lower)
        value2 = ersemdict[ersemvar_lower]
        if '.' in value2 or 'e' in value2.lower():
            value2 = float(value2)
        else:
            value2 = int(value2)
        if options.verbose: 
            print '%s%s = %s (%s/%s)' % (prefix,name,value2,ersemsources[ersemvar_lower],value[1:])
        d[name] = value2

processDictionary(root)

print 'Saving YAML...'

if options.clean:
    # Write YAML to temporary file.
    f = tempfile.NamedTemporaryFile('w',delete=False)
    yaml.dump(root,f, default_flow_style=False, indent=2)
    f.close()

    # Clean-up file
    try:
       import pyfabm.complete_yaml
    except ImportError:
       print 'Unable to load pyfabm. Please build and install python-fabm by running cmake+make for source directory $FABMDIR/src/drivers/python.'
       sys.exit(1)
    print 'Cleaning YAML...'
    pyfabm.complete_yaml.processFile(f.name,outfile,subtract_background=options.subtract_background)

    os.remove(f.name) 
else:
    with open(outfile,'w') as f:
        yaml.dump(root,f, default_flow_style=False, indent=2)

print 'Done. Written %s.' % outfile
