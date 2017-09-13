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
from lib.prediction_class import *

def resolve_filepath(filepath):
    if int_pattern.match(filepath) and int(filepath) <= len(current_app.config['storage']['manifest']):
        filepath = current_app.config['storage']['manifest'][int(filepath)]
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
    additional_hashes = {}
    if configObj['additional_input_file_list'] != "":
        additional_hashes = {
            key:hashfile(path)
            for (key, path) in yaml.load(
                open(configObj['additional_input_file_list'])
            ).items()
        }
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
                #if there are keys in one set and not the other
                continue
            #otherwise, check every key except input, output, and additional files
            failed = False
            for param in configObj:
                if param in {'input', 'output', 'additional_input_file_list'}:
                    #Skip for now, because these will need a longer check
                    continue
                if configObj[param] != current_config[param]:
                    failed = True
                    break
            if not failed:
                #do longer checks on input and additional files
                current_hash = hashfile(current_config['input'])
                if current_hash != input_hash:
                    continue
                current_input_hashes = {}
                if current_config['additional_input_file_list'] != "":
                    current_input_hashes = {
                        key:hashfile(path)
                        for (key, path) in yaml.load(
                            open(current_config['additional_input_file_list'])
                        ).items()
                    }
                if len(set(additional_hashes)^set(current_input_hashes)):
                    continue
                if reduce(
                    lambda x,y: x and additional_hashes[y]==current_input_hashes[y],
                    current_input_hashes,
                    True
                ):
                    return i
    return None

def staging(input, samplename, alleles, epitope_lengths, prediction_algorithms,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          normal_snvs_coverage_file, normal_indels_coverage_file,
          tdna_snvs_coverage_file, tdna_indels_coverage_file,
          trna_snvs_coverage_file, trna_indels_coverage_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold, fasta_size, iedb_retries, iedb_install_dir,
          downstream_sequence_length, keep_tmp_files, force):
    """Stage input for a new pVAC-Seq run.  Generate a unique output directory and \
    save uploaded files to temporary locations (and give pVAC-Seq the filepaths). \
    Then forward the command to start()"""
    input_file = input
    data = current_app.config['storage']['loader']()
    samplename  = re.sub(r'[^\w\s.]', '_', samplename)
    list_input() #update the manifest stored in current_app
    # input_manifest = current_app.config['storage']['manifest']
    current_path = os.path.join(
        current_app.config['files']['data-dir'],
        'results',
        samplename
    )
    if os.path.exists(current_path):
        i = 1
        while os.path.exists(current_path+"_"+str(i)):
            i += 1
        current_path += "_"+str(i)

    temp_path = tempfile.TemporaryDirectory()

    input_path = resolve_filepath(input_file)
    if not input_path:
        return (
            {
                'code':400,
                'message':'Unable to locate the given file: %s'%input_file,
                'fields':'input'
            },400
        )

    additional_input_file_list = open(os.path.join(temp_path.name, "additional_input_file_list.yml"), 'w')

    if gene_expn_file:
        gene_expn_file_path = resolve_filepath(gene_expn_file)
        if not gene_expn_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%gene_expn_file,
                    'fields':'gene_expn_file'
                },400
            )
        if os.path.getsize(gene_expn_file_path):
            yaml.dump({"gene_expn_file": gene_expn_file_path}, additional_input_file_list, default_flow_style=False)

    if transcript_expn_file:
        transcript_expn_file_path = resolve_filepath(transcript_expn_file)
        if not transcript_expn_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%transcript_expn_file,
                    'fields':'transcript_expn_file'
                },400
            )
        if os.path.getsize(transcript_expn_file_path):
            yaml.dump({"transcript_expn_file" :transcript_expn_file_path}, additional_input_file_list, default_flow_style=False)

    if normal_snvs_coverage_file:
        normal_snvs_coverage_file_path = resolve_filepath(normal_snvs_coverage_file)
        if not normal_snvs_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%normal_snvs_coverage_file,
                    'fields':'normal_snvs_coverage_file'
                },400
            )
        if os.path.getsize(normal_snvs_coverage_file_path):
            yaml.dump({"normal_snvs_coverage_file" :normal_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    if normal_indels_coverage_file:
        normal_indels_coverage_file_path = resolve_filepath(normal_indels_coverage_file)
        if not normal_indels_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%normal_indels_coverage_file,
                    'fields':'normal_indels_coverage_file'
                },400
            )

        if os.path.getsize(normal_indels_coverage_file_path):
            yaml.dump({"normal_indels_coverage_file" :normal_indels_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    if tdna_snvs_coverage_file:
        tdna_snvs_coverage_file_path = resolve_filepath(tdna_snvs_coverage_file)
        if not tdna_snvs_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%tdna_snvs_coverage_file,
                    'fields':'tdna_snvs_coverage_file'
                },400
            )
        if os.path.getsize(tdna_snvs_coverage_file_path):
            yaml.dump({"tdna_snvs_coverage_file" :tdna_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    if tdna_indels_coverage_file:
        tdna_indels_coverage_file_path = resolve_filepath(tdna_indels_coverage_file)
        if not tdna_indels_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%tdna_indels_coverage_file,
                    'fields':'tdna_indels_coverage_file'
                },400
            )

        if os.path.getsize(tdna_indels_coverage_file_path):
            yaml.dump({"tdna_indels_coverage_file" :tdna_indels_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    if trna_snvs_coverage_file:
        trna_snvs_coverage_file_path = resolve_filepath(trna_snvs_coverage_file)
        if not trna_snvs_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%trna_snvs_coverage_file,
                    'fields':'trna_snvs_coverage_file'
                },400
            )

        if os.path.getsize(trna_snvs_coverage_file_path):
            yaml.dump({"trna_snvs_coverage_file" :trna_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    if trna_indels_coverage_file:
        trna_indels_coverage_file_path = resolve_filepath(trna_indels_coverage_file)
        if not trna_indels_coverage_file_path:
            return (
                {
                    'code':400,
                    'message':'Unable to locate the given file: %s'%trna_indels_coverage_file,
                    'fields':'trna_indels_coverage_file'
                },400
            )

        if os.path.getsize(trna_indels_coverage_file_path):
            yaml.dump({"trna_indels_coverage_file" :trna_indels_coverage_file_path}, additional_input_file_list, default_flow_style=False)

    additional_input_file_list.flush()

    configObj = {
        'input': input_path, #input
        'samplename':samplename, #samplename
        'alleles':alleles.split(','),
        'output':current_path,
        'epitope_lengths':[int(item) for item in epitope_lengths.split(',')],
        'prediction_algorithms':prediction_algorithms.split(','),
        'peptide_sequence_length':peptide_sequence_length,
        'additional_input_file_list':(
            additional_input_file_list.name if additional_input_file_list.tell() else '' #check if any data was written to file
        ),
        'net_chop_method':net_chop_method,
        'netmhc_stab':bool(netmhc_stab),
        'top_result_per_mutation':bool(top_result_per_mutation),
        'top_score_metric':top_score_metric,
        'binding_threshold':binding_threshold,
        'minimum_fold_change':minimum_fold_change,
        'normal_cov':normal_cov, #normal_cov
        'tdna_cov':tdna_cov, #tdna_cov
        'trna_cov':trna_cov, #trna_cov
        'normal_vaf':normal_vaf, #normal_vaf
        'tdna_vaf':tdna_vaf, #tdna_vaf
        'trna_vaf':trna_vaf, #trna_vaf
        'expn_val':expn_val, #expn_val
        'net_chop_threshold':net_chop_threshold,
        'fasta_size':fasta_size,
        'iedb_retries':iedb_retries,
        'iedb_install_dir':iedb_install_dir,
        'keep_tmp_files':bool(keep_tmp_files),
        'downstream_sequence_length':downstream_sequence_length
    }
    checkOK = precheck(configObj, data) if not force else None
    if checkOK is None:
        copytree(temp_path.name, current_path)
        print(additional_input_file_list.tell())
        if configObj['additional_input_file_list']:
            configObj['additional_input_file_list'] = os.path.join(
                current_path,
                os.path.basename(additional_input_file_list.name)
            )
        writer = open(os.path.join(
            os.path.abspath(current_path),
            'config.json'
        ),'w')
        json.dump(configObj, writer, indent='\t')
        writer.close()
        temp_path.cleanup()
        new_id = start(**configObj)

        return ({
            'code':201,
            'message': "Process started.",
            'processid': new_id
        },
        201)
    
    return (
        {
            'code':400,
            'message':"The given parameters match process %d"%checkOK,
            'fields':"N/A"
        },
        400
    )


def start(input, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, additional_input_file_list,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
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
        '-s', str(fasta_size),
        '-r', str(iedb_retries),
        '-d', str(downstream_sequence_length)
    ]
    if len(additional_input_file_list):
        command += ['-i', additional_input_file_list]
    if len(net_chop_method):
        command += [
            '--net-chop-method', net_chop_method,
            '--net-chop-threshold', str(net_chop_threshold)
        ]
    if netmhc_stab:
        command.append('--netmhc-stab')
    if top_result_per_mutation:
        command.append('--top-result-per-mutation')
    if keep_tmp_files:
        command.append('-k')
    if len(iedb_install_dir):
        command += [
            '--iedb-install-directory',
            iedb_install_dir
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

# takes in comma delimited string of prediction algorithms,
# returns map of algorithms to valid alleles for that algorithm
def valid_alleles(prediction_algorithms):
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
