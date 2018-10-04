import os
import subprocess
import json
import re
import tempfile
from functools import reduce
import hashlib
from flask import current_app
import yaml
from shlex import quote
from shutil import copyfile, copytree
from .database import int_pattern
from .files import list_input
from .utils import fullresponse
from lib.prediction_class import *

def resolve_filepath(filepath):
    data = current_app.config['storage']['loader']()
    if int_pattern.match(filepath) and str(filepath) in data['input']:
        filepath = data['input'][str(filepath)]['fullname']
    if not os.path.isfile(filepath):
        filepath = os.path.join(
            current_app.config['files']['data-dir'],
            'input',
            filepath
        )
        if not os.path.isfile(filepath):
            return None
    return filepath

def hashfile(filepath):
    try:
        reader = open(filepath, mode='rb')
    except FileNotFoundError:
        return None
    chunk = reader.read(4096)
    hasher = hashlib.md5()
    while chunk:
        hasher.update(chunk)
        chunk = reader.read(4096)
    reader.close()
    return hasher.digest()

def precheck(configObj, data):
    """Check if a process with these same parameters has already been run successfully"""
    #This is a temprary stand-in until https://github.com/griffithlab/pVAC-Seq/pull/292 is merged
    input_hash = hashfile(configObj['input'])
    for i in range(data['processid']+1):
        key = 'process-%d'%i
        if key in data:
            reader = open(os.path.join(
                data[key]['output'],
                'config.json'
            ))
            current_config = json.load(reader)
            reader.close()
            #input
            if len(set(configObj)^set(current_config)):
                # if there are keys in one set and not the other
                continue
            #otherwise, check every key except input and output
            failed = False
            for param in configObj:
                if param in {'input', 'output'}:
                    #Skip for now, because these will need a longer check
                    continue
                if configObj[param] != current_config[param]:
                    failed = True
                    break
            if not failed:
                #do longer checks on input
                current_hash = hashfile(current_config['input'])
                if current_hash != input_hash:
                    continue
                return i
    return None

def staging(parameters):
    """Stage input for a new pVAC-Seq run.  Generate a unique output directory and \
    save uploaded files to temporary locations (and give pVAC-Seq the filepaths). \
    Then forward the command to start()"""
    input_file = parameters['input']
    data = current_app.config['storage']['loader']()
    samplename = re.sub(r'[^\w\s.]', '_', parameters['samplename'])
    list_input()  # update the manifest stored in current_app
    # input_manifest = current_app.config['storage']['manifest']
    current_path = os.path.join(
        current_app.config['files']['data-dir'],
        '.processes',
        samplename
    )

    if os.path.exists(current_path):
        i = 1
        while os.path.exists(current_path + "_" + str(i)):
            i += 1
        current_path += "_" + str(i)

    temp_path = tempfile.TemporaryDirectory()

    input_path = resolve_filepath(input_file)
    if not input_path:
        return (
            {
                'status': 400,
                'message': 'Unable to locate the given file: %s' % input_file,
                'fields': 'input'
            }, 400
        )
    phased_proximal_variants_vcf = parameters.pop('phased_proximal_variants_vcf', "")
    if len(phased_proximal_variants_vcf):
        phased_proximal_variants_vcf_path = resolve_filepath(phased_proximal_variants_vcf)
        if not phased_proximal_variants_vcf_path:
            return (
                {
                    'status': 400,
                    'message': 'Unable to locate the given file: %s' % phased_proximal_variants_vcf,
                    'fields': 'phased_proximal_variants_vcf'
                }, 400
            )
    else:
        phased_proximal_variants_vcf_path = ""

    if 'epitope_lengths' in parameters:
        epitope_lengths = [int(item) for item in parameters['epitope_lengths'].split(',')]
    else:
        epitope_lengths = ""

    configObj = {
        'input': input_path,  # input
        'phased_proximal_variants_vcf': phased_proximal_variants_vcf_path,
        'samplename': samplename,  # samplename
        'alleles': parameters['alleles'].split(','),
        'output': current_path,
        'epitope_lengths': epitope_lengths,
        'prediction_algorithms': parameters['prediction_algorithms'].split(','),
        'peptide_sequence_length': parameters.pop('peptide_sequence_length', 21),
        'net_chop_method': parameters.pop('net_chop_method', ""),
        'netmhc_stab': bool(parameters.pop('netmhc_stab', False)),
        'pass_only': bool(parameters.pop('pass_only', False)),
        'allele_specific_cutoffs': bool(parameters.pop('allele_specific_cutoffs', False)),
        'top_score_metric': parameters.pop('top_score_metric', 'median'),
        'binding_threshold': parameters.pop('binding_threshold', 500),
        'minimum_fold_change': parameters.pop('minimum_fold_change', 0),
        'normal_cov': parameters.pop('normal_cov', 5),  # normal_cov
        'tdna_cov': parameters.pop('tdna_cov', 10),  # tdna_cov
        'trna_cov': parameters.pop('trna_cov', 10),  # trna_cov
        'normal_vaf': parameters.pop('normal_vaf', 0.02),  # normal_vaf
        'tdna_vaf': parameters.pop('tdna_vaf', 0.25),  # tdna_vaf
        'trna_vaf': parameters.pop('trna_vaf', 0.25),  # trna_vaf
        'expn_val': parameters.pop('expn_val', 1),  # expn_val
        'net_chop_threshold': parameters.pop('net_chop_threshold', 0.5),
        'fasta_size': parameters.pop('fasta_size', 200),
        'maximum_transcript_support_level': parameters.pop('maximum_transcript_support_level', 1),
        'iedb_retries': parameters.pop('iedb_retries', 5),
        'iedb_install_dir': parameters.pop('iedb_install_dir', ""),
        'keep_tmp_files': bool(parameters.pop('keep_tmp_files', False)),
        'downstream_sequence_length': parameters.pop('downstream_sequence_length', 'full')
    }
    force = bool(parameters.pop('force', False))
    checkOK = precheck(configObj, data) if not force else None
    if checkOK is None:
        copytree(temp_path.name, current_path)
        writer = open(os.path.join(
            os.path.abspath(current_path),
            'config.json'
        ), 'w')
        json.dump(configObj, writer, indent='\t')
        writer.close()
        temp_path.cleanup()
        new_id = start(**configObj)

        return ({
            'status': 201,
            'message': "Process started.",
            'processid': new_id
        }, 201)
    
    return (
        {
            'status': 400,
            'message': "The given parameters match process %d" % checkOK,
            'fields': "N/A"
        }, 400)


def start(input, phased_proximal_variants_vcf, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, net_chop_method, netmhc_stab, pass_only, top_score_metric,
          binding_threshold, allele_specific_cutoffs, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf, maximum_transcript_support_level,
          expn_val, net_chop_threshold, fasta_size, iedb_retries, iedb_install_dir,
          downstream_sequence_length, keep_tmp_files):
    """Build the command for pVAC-Seq, then spawn a new process to run it"""
    if type(epitope_lengths) == list:
        epitope_lengths = ','.join(str(item) for item in epitope_lengths)
    if type(alleles) == list:
        alleles = ','.join(str(item) for item in alleles)
    if type(prediction_algorithms) == list:
        prediction_algorithms = ','.join(str(item) for item in prediction_algorithms)
    command = [
        'pvacseq',
        'run',
        input,
        samplename,
        alleles
    ]
    for algo in prediction_algorithms.split(','):
        command.append(algo)
    command += [
        output,
        '-e', epitope_lengths,
        '-l', str(peptide_sequence_length),
        '-m', top_score_metric,
        '-b', str(binding_threshold),
        '-c', str(minimum_fold_change),
        '--normal-cov', str(normal_cov),
        '--tdna-cov', str(tdna_cov),
        '--trna-cov', str(trna_cov),
        '--normal-vaf', str(normal_vaf),
        '--tdna-vaf', str(tdna_vaf),
        '--trna-vaf', str(trna_vaf),
        '--expn-val', str(expn_val),
        '--maximum-transcript-support-level', str(maximum_transcript_support_level),
        '-s', str(fasta_size),
        '-r', str(iedb_retries),
        '-d', str(downstream_sequence_length)
    ]
    if len(net_chop_method):
        command += [
            '--net-chop-method', net_chop_method,
            '--net-chop-threshold', str(net_chop_threshold)
        ]
    if netmhc_stab:
        command.append('--netmhc-stab')
    if allele_specific_cutoffs:
        command.append('--allele-specific-binding-thresholds')
    if keep_tmp_files:
        command.append('-k')
    if pass_only:
        command.append('--pass-only')
    if len(iedb_install_dir):
        command += [
            '--iedb-install-directory',
            iedb_install_dir
        ]
    if len(phased_proximal_variants_vcf):
        command +=[
            '--phased-proximal-variants-vcf',
            phased_proximal_variants_vcf,
        ]

    # stdout and stderr from the child process will be directed to this file
    logfile = os.path.join(output, 'pVAC-Seq.log')
    with current_app.config['storage']['synchronizer']:
        data = current_app.config['storage']['loader']()
        data['processid']+=1
        os.makedirs(os.path.dirname(logfile), exist_ok = True)
        current_app.config['storage']['children'][data['processid']] = subprocess.Popen(
            command,
            stdout=open(logfile, 'w'),  # capture stdout in the logfile
            stderr=subprocess.STDOUT,
            # isolate the child in a new process group
            # this way it will remainin running no matter what happens to this process
            preexec_fn=os.setpgrp
        )
        # Store some data about the child process
        data.addKey(
            'process-%d'%(data['processid']),
            {
                'command': " ".join([quote(token) for token in command]),
                'logfile':logfile,
                'pid':current_app.config['storage']['children'][data['processid']].pid,
                'status': 0,
                'files':{},
                'output':os.path.abspath(output)
            },
            current_app.config['files']['processes']
        )
        if 'reboot' not in data:
            data.addKey(
                'reboot',
                current_app.config['reboot'],
                current_app.config['files']['processes']
            )
        data.save()
        return data['processid']

def test():
    """Return the submission page (a stand-in until there is a proper ui for submission)"""
    reader = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_start.html'))
    data = reader.read()
    reader.close()
    return data


def check_allele(allele):
    """Checks if the requested allele is supported by pVAC-Seq or not"""
    data = current_app.config['storage']['loader']()
    if 'allele_file' not in current_app.config:
        current_app.config['allele_file'] = list(PredictionClass.all_valid_allele_names())
    for line in current_app.config['allele_file']:
        if line.strip() == allele:
            return True
    return False

# returns map of alleles to prediction algorithms it can be used with
def valid_alleles(page, count, prediction_algorithms=None, name_filter=None):
    data = PredictionClass.allele_info(prediction_algorithms, name_filter)
    return fullresponse(data, page, count)

# takes in comma delimited string of prediction algorithms,
# returns map of algorithms to valid alleles for that algorithm
def valid_alleles_per_algorithm(prediction_algorithms):
    valid_allele_list = {}
    for algorithm in prediction_algorithms.split(","):
        prediction_class = globals()[algorithm]
        alleles = prediction_class().valid_allele_names()
        # alleles sometimes returns as dict_keys instead of an array, so must specify as list 
        valid_allele_list[algorithm] = list(alleles)
    return valid_allele_list

# naming prediction_algorithms to keep consistent with the pVac-Seq documentation
def prediction_algorithms():
    return PredictionClass.prediction_methods()
