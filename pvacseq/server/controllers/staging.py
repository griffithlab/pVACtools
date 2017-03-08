import os
import subprocess
import json
import re
import tempfile
from flask import current_app
from yaml import dump
from shlex import quote
from shutil import copyfile
from .database import int_pattern
from .files import list_input

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

def staging(input, samplename, alleles, epitope_lengths, prediction_algorithms,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          normal_snvs_coverage_file, normal_indels_coverage_file,
          tdna_snvs_coverage_file, tdna_indels_coverage_file,
          trna_snvs_coverage_file, trna_indels_coverage_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold, fasta_size, iedb_retries, downstream_sequence_length, keep_tmp_files):
    """Stage input for a new pVAC-Seq run.  Generate a unique output directory and \
    save uploaded files to temporary locations (and give pVAC-Seq the filepaths). \
    Then forward the command to start()"""
    ['input','gene_expn_file', 'transcript_expn_file',
    'normal_snvs_coverage_file', 'normal_indels_coverage_file',
    'tdna_snvs_coverage_file', 'tdna_indels_coverage_file',
    'trna_snvs_coverage_file', 'trna_indels_coverage_file']
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

    input_path = resolve_filepath(input_file)
    if not input_path:
        return (
            {
                'code':400,
                'message':'Unable to locate the given file: %s'%input_file,
                'fields':'input'
            },400
        )

    os.makedirs(current_path)
    additional_input_file_list = open(os.path.join(current_path, "additional_input_file_list.yml"), 'w')

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
            dump({"gene_expn_file": gene_expn_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"transcript_expn_file" :transcript_expn_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"normal_snvs_coverage_file" :normal_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"normal_indels_coverage_file" :normal_indels_coverage_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"tdna_snvs_coverage_file" :tdna_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"tdna_indels_coverage_file" :tdna_indels_coverage_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"trna_snvs_coverage_file" :trna_snvs_coverage_file_path}, additional_input_file_list, default_flow_style=False)

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
            dump({"trna_indels_coverage_file" :trna_indels_coverage_file_path}, additional_input_file_list)

    additional_input_file_list.flush()

    return start(input_path, samplename, alleles, epitope_lengths, prediction_algorithms, current_path,
              peptide_sequence_length, additional_input_file_list.name if additional_input_file_list.tell() else "", # check if any data written to file
              net_chop_method, len(netmhc_stab), len(top_result_per_mutation), top_score_metric,
              binding_threshold, minimum_fold_change,
              normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
              expn_val, net_chop_threshold,
              fasta_size, iedb_retries, downstream_sequence_length, len(keep_tmp_files))


def start(input, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, additional_input_file_list,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold,
          fasta_size, iedb_retries, downstream_sequence_length, keep_tmp_files):
    """Build the command for pVAC-Seq, then spawn a new process to run it"""
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

    # stdout and stderr from the child process will be directed to this file
    logfile = os.path.join(output, 'pVAC-Seq.log')
    current_app.config['storage']['synchronizer'].acquire()
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
            'status': "Task Started",
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
    current_app.config['storage']['synchronizer'].release()
    configObj = {
        'action':'run',
        'input_file': input,
        'sample_name':samplename,
        'alleles':alleles.split(','),
        'prediction_algorithms':prediction_algorithms.split(','),
        'output_directory':output,
        'epitope_lengths':epitope_lengths.split(','),
        'peptide_sequence_length':peptide_sequence_length,
        'additional_input_files':additional_input_file_list.split(','),
        'net_chop_method':net_chop_method,
        'netmhc_stab':netmhc_stab,
        'top_result_per_mutation':top_result_per_mutation,
        'top_score_metric':top_score_metric,
        'binding_threshold':binding_threshold,
        'minimum_fold_change':minimum_fold_change,
        'normal_coverage_cutoff':normal_cov,
        'tumor_dna_coverage_cutoff':tdna_cov,
        'tumor_rna_coverage_cutoff':trna_cov,
        'normal_vaf_cutoff':normal_vaf,
        'tumor_dna_vaf_cutoff':tdna_vaf,
        'tumor_rna_vaf_cutoff':trna_vaf,
        'expression_cutoff':expn_val,
        'netchop_threshold':net_chop_threshold,
        'fasta_size':fasta_size,
        'iedb_retries':iedb_retries,
        'downstream_sequence_length':downstream_sequence_length
    }

    writer = open(os.path.join(
        os.path.abspath(output),
        'config.json'
    ),'w')
    json.dump(configObj, writer, indent='\t')
    writer.close()
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
        allele_file = tempfile.TemporaryFile('w+')
        subprocess.call(['pvacseq', 'valid_alleles'], stdout=allele_file)
        current_app.config['allele_file'] = allele_file
    current_app.config['allele_file'].seek(0)
    for line in current_app.config['allele_file']:
        if line.strip() == allele:
            return True
    return False
