import os
import csv
import sys
import itertools
from flask import session
import multiprocessing as mp
import subprocess
import re
import tempfile

spinner = re.compile(r'[\\\b\-/|]{2,}')

children = []
logs = []
staging_files = []

descriptions = {
    'chop.tsv':"Processed and filtered data, with peptide cleavage data added",
    'combined.parsed.tsv':"Processed data from IEDB, but with no filtering or extra data",
    'filtered.binding.tsv':"Processed data filtered by binding strength",
    'filtered.coverage.tsv':"Processed data filtered by binding strength and coverage",
    'stab.tsv':"Processed and filtered data, with peptide stability data added",
    'final.tsv':"Final output data",
    'tsv':"Raw input data parsed out of the input vcf"
}

def column_filter(id, column):
    return column.replace(' ', '_').lower().strip()

def gen_files_list(id):
    if 'files' not in session['process-%d'%id]:
        session['process-%d'%id]['files'] = []
        base_dir = session['process-%d'%id]['output']
        if os.path.isdir(os.path.join(base_dir, 'class_i')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'class_i'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_i', path)):
                    session['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_i', path))
        if os.path.isdir(os.path.join(base_dir, 'class_ii')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'class_ii'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_ii', path)):
                    session['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_ii', path))
        for path in sorted(os.listdir(base_dir)):
            if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, path)):
                session['process-%d'%id]['files'].append(os.path.join(base_dir, path))

def results_get(id):
    processKey = 'process-%d'%id
    if processKey in session and children[id][1].is_alive():
        return []
    output = []
    gen_files_list(id)
    for fileID in range(len(session['process-%d'%id]['files'])):
        extension = '.'.join(os.path.basename(session['process-%d'%id]['files'][fileID]).split('.')[1:])
        output.append({
            'fileID':fileID,
            'description':descriptions[extension] if extension in descriptions else "Unknown file",
            'display_name':os.path.relpath(session['process-%d'%id]['files'][fileID], session['process-%d'%id]['output']),
            'url':'/api/v1/processes/%d/results/%d'%(id, fileID),
            'size':"%0.3f KB"%(os.path.getsize(os.path.join(session['process-%d'%id]['output'], session['process-%d'%id]['files'][fileID]))/1024)
        })
    return output

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

def staging(input, samplename, alleles, epitope_lengths, prediction_algorithms,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change, expn_val, net_chop_threshold,
          fasta_size, keep_tmp_files):
    staged_input = tempfile.NamedTemporaryFile('wb')
    input.save(staged_input)
    staged_input.flush()
    staging_files.append(staged_input)

    staged_gene_expn_file = tempfile.NamedTemporaryFile('wb')
    gene_expn_file.save(staged_gene_expn_file)
    staged_gene_expn_file.flush()
    staging_files.append(staged_gene_expn_file)

    staged_transcript_expn_file = tempfile.NamedTemporaryFile('wb')
    transcript_expn_file.save(staged_transcript_expn_file)
    staged_transcript_expn_file.flush()
    staging_files.append(staged_transcript_expn_file)

    current_path = os.path.join(os.path.expanduser('~'), "Documents", "pVAC-Seq Output", samplename)
    if os.path.exists(current_path):
        i = 1
        while os.path.exists(current_path+"_"+str(i)):
            i+=1
        current_path += "_"+str(i)

    return start(staged_input.name, samplename, alleles, epitope_lengths, prediction_algorithms, current_path,
              peptide_sequence_length, staged_gene_expn_file.name if staged_gene_expn_file.tell() else "",
              staged_transcript_expn_file.name if staged_transcript_expn_file.tell() else "", #check if any data written to file
              net_chop_method, len(netmhc_stab), len(top_result_per_mutation), top_score_metric,
              binding_threshold, minimum_fold_change, expn_val, net_chop_threshold,
              fasta_size, len(keep_tmp_files))

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
    return [{
        'id':i,
        'running':children[i][1].is_alive()
    } for i in range(len(children)) if 'process-%d'%i in session]

def process_info(id):
    intake = b''
    while children[id][0].poll():
        intake += os.read(children[id][0].fileno(), 512)
    intake = spinner.sub('', intake.decode()).strip()
    logs[id]+=intake
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
    # children[id][0].close()
    return status

def shutdown():
    output = []
    for i in range(len(children)):
        if children[i][1].is_alive():
            output.append(i)
        children[i][1].join(.1)
        if children[i][1].is_alive():
            children[i][1].terminate()
    return output

def test():
    reader = open("/Users/agrauber/Documents/pVAC-Seq/pvacseq/server/test_start.html")
    data = reader.read()
    reader.close()
    return data

def _process_worker(pipe, command):
    from ...lib.main import main
    # setup this process' stdout and stderr to write to the pipe's file descriptor
    sys.stdout = open(pipe.fileno(), mode='w')
    sys.stderr = sys.stdout
    # run the pvacseq command
    main(command)
