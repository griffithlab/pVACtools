import os
import csv
import sys
import itertools
from flask import session
import multiprocessing as mp
import subprocess
import re

spinner = re.compile(r'[\\\b\-/|]{2,}')

children = []
logs = []

descriptions = {
    'chop.tsv':"Processed and filtered data, with peptfileIDe cleavage data added",
    'combined.parsed.tsv':"Processed data from IEDB, but with no filtering or extra data",
    'filtered.binding.tsv':"Processed data filtered by binding strength",
    'filtered.coverage.tsv':"Processed data filtered by binding strength and coverage",
    'stab.tsv':"Processed and filtered data, with peptfileIDe stability data added",
    'final.tsv':"Final output data",
    'tsv':"Raw input data parsed out of the input vcf"
}

def column_filter(id, column):
    if 'columns' not in session['process-%d'%id]:
        session['process-%d'%id]['columns'] = {}
    result = column.replace(' ', '_').lower().strip()
    if column not in session['process-%d'%id]['columns']:
        session['process-%d'%id]['columns'][column]=result
    return result

def gen_files_list(id):
    if 'files' not in session['process-%d'%id]:
        session['process-%d'%id]['files'] = []
        base_dir = session['process-%d'%id]['output']
        for path in os.listdir(os.path.join(base_dir, 'class_i')):
            if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_i', path)):
                session['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_i', path))
        for path in os.listdir(os.path.join(base_dir, 'class_ii')):
            if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_ii', path)):
                session['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_ii', path))
        for path in os.listdir(base_dir):
            if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, path)):
                session['process-%d'%id]['files'].append(os.path.join(base_dir, path))

def results_get(id, count = None, page = None):
    processKey = 'process-%d'%id
    if processKey in session and children[id][1].is_alive():
        return []
    output = []
    gen_files_list(id)
    for fileID in range(len(session['process-%d'%id]['files'])):
        extension = '.'.join(os.path.basename(session['process-%d'%id]['files'][fileID]).split('.')[1:])
        display_name = os.path.join(
            os.path.dirname(session['process-%d'%id]['files'][fileID]),
            os.path.basename(session['process-%d'%id]['files'][fileID])
        )
        if 'class_i' not in display_name:
            display_name = os.path.basename(display_name)
        output.append({
            'fileID':fileID,
            'description':descriptions[extension] if extension in descriptions else "Unknown file",
            'display_name':display_name,
            'url':'/api/v1/processes/%d/results/%d'%(id, fileID),
            'size':os.path.getsize(os.path.join(session['process-%d'%id]['output'], session['process-%d'%id]['files'][fileID]))
        })
    return output[(page-1)*count:page*count]

def results_getfile(id, count = None, page = None, fileID = None):
    processKey = 'process-%d'%id
    if processKey in session and children[id][1].is_alive():
        return []
    gen_files_list(id)
    raw_reader = open(session['process-%d'%id]['files'][fileID])
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = [{column_filter(id, k):entry[k] for k in entry} for entry in itertools.islice(reader, (page-1)*count, page*count)]
    raw_reader.close()
    return output

def results_getcols(id, fileID):
    processKey = 'process-%d'%id
    if processKey in session and children[id][1].is_alive():
        return {}
    gen_files_list(id)
    raw_reader = open(session['process-%d'%id]['files'][fileID])
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = {column_filter(id, field):field for field in reader.fieldnames}
    raw_reader.close()
    return output

def start(input, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change, expn_val, net_chop_threshold,
          fasta_size, keep_tmp_files):
    command = [
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
        '--expn-val', str(expn_val),
        '-s', str(fasta_size)
    ]
    if len(gene_expn_file):
        command+= ['-g', gene_expn_file]
    if len(transcript_expn_file):
        command+=['-i', transcript_expn_file]
    if net_chop_method:
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
    (parent, child) = mp.Pipe()
    children.append([
        parent,
        mp.Process(target=_process_worker, args=(child, command))
    ])
    children[-1][1].start()
    session['process-%d'%(len(children)-1)] = {
        'command': "pvacseq run "+" ".join(command),
        'status': "Task Started",
        'log': "",
        'output':os.path.abspath(output)
    }
    logs.append("")
    return len(children)-1

def processes():
    return [i for i in range(len(children)) if 'process-%d'%i in session]

def process_info(id):
    intake = b''
    while children[id][0].poll():
        intake += os.read(children[id][0].fileno(), 512)
    intake = spinner.sub('', intake.decode()).strip().split("\n")
    logs[id]+='\n'.join(intake)
    session['process-%d'%id]['status'] = logs[id].split("\n")[-1]
    if not children[id][1].is_alive():
        session['process-%d'%id]['status'] = "Process complete: %d"%children[id][1].exitcode
    return {
        'pid':children[id][1].pid,
        'id':id,
        'command':session['process-%d'%id]['command'],
        'status':session['process-%d'%id]['status'],
        'log':logs[id],
        'output':session['process-%d'%id]['output'],
        'running':children[id][1].is_alive()
    }

def stop(id):
    status = process_info(id)
    children[id][1].terminate()
    children[id][0].close()
    return status

def shutdown():
    for i in range(len(children)):
        children[i][1].join(.1)
        if children[i][1].is_alive():
            children[i][1].terminate()
    sys.exit("Server closed")
    return {}

def _process_worker(pipe, command):
    from ...lib.main import main
    # setup this process' stdout and stderr to write to the pipe's file descriptor
    sys.stdout = open(pipe.fileno(), mode='w')
    sys.stderr = sys.stdout
    # run the pvacseq command
    main(command)
