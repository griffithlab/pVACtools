import os
import csv
import sys
import itertools
import subprocess
import re
import tempfile
import json
import shutil
from yaml import dump

spinner = re.compile(r'[\\\b\-/|]{2,}')
allele_file = None

children = {}  # Holds popen objects for processes spawned by this current session

# path to the configuration file where the server stores its data between runs
configfile = os.path.join(os.path.expanduser('~'), '.pvacseq_ui')

# check the last system reboot, because it means that the logged pid's can now
# be used by new processes
reboot = subprocess.check_output(['last', 'reboot']).decode().split("\n")[0]

# fetch data from the configuration file
if os.path.isfile(configfile):
    print("Resuming from saved state")
    data = json.load(open(configfile))
    if 'reboot' in data and data['reboot'] != reboot:
        print("A reboot has occurred since the server was first started")
        print(
            "pid's of old pVAC-Seq runs with id's",
            data['processid'],
            "and lower may be innacurate"
        )
else:
    print("No saved state found")
    data={
        'processid':-1
    }

# a mapping to provide a description of each result file based on its extension
descriptions = {
    'chop.tsv':"Processed and filtered data, with peptide cleavage data added",
    'combined.parsed.tsv':"Processed data from IEDB, but with no filtering or extra data",
    'filtered.binding.tsv':"Processed data filtered by binding strength",
    'filtered.coverage.tsv':"Processed data filtered by binding strength and coverage",
    'stab.tsv':"Processed and filtered data, with peptide stability data added",
    'final.tsv':"Final output data",
    'tsv':"Raw input data parsed out of the input vcf"
}

def column_filter(column):
    """standardize column names"""
    return column.replace(' ', '_').lower().strip()

def gen_files_list(id):
    """Generate the list of result files for a given process.  Stash them for later use"""
    #  if 'process-%d'%id not in data:
    #     raise KeyError("The requested process (%d) does not exist"%id)
    if 'files' not in data['process-%d'%id]:
        data['process-%d'%id]['files'] = []
        base_dir = data['process-%d'%id]['output']
        if os.path.isdir(os.path.join(base_dir, 'class_i')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'class_i'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_i', path)):
                    data['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_i', path))
        if os.path.isdir(os.path.join(base_dir, 'class_ii')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'class_ii'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'class_ii', path)):
                    data['process-%d'%id]['files'].append(os.path.join(base_dir, 'class_ii', path))
        for path in sorted(os.listdir(base_dir)):
            if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, path)):
                data['process-%d'%id]['files'].append(os.path.join(base_dir, path))
        savedata()

def fetch_process(id):
    """Fetch a tuple of available process info"""
    return (
        # This entry should always exist.  Any spawned process will have an entry in data
        data['process-%d'%id] if 'process-%d'%id in data else {},
        # This will only exist if spawned by the current run of the server
        # If the server crashed and restarted, the child is orphaned,
        # and the server will only be able to read its output files
        # (unable to access the return code or terminate it)
        children[id] if id in children else False
    )

def is_running(id):
    """Returns True if the requested process looks like it's still running"""
    process = fetch_process(id)
    if not process[0]:
        return False  # The process doesn't exist
    if process[1]:
        return process[1].poll() == None
    try:
        # check if the process is active by sending a dummy signal
        os.kill(process[0]['pid'], 0)
    except ProcessLookupError:
        return False
    return True

def savedata():
    """Saves the data object to the configfile"""
    writer = open(configfile, 'w')
    json.dump(data, writer)
    writer.close()


def results_get(id):
    """Get the list of result files from a specific pVAC-Seq run"""
    process = fetch_process(id)
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    if is_running(id):
        return []
    gen_files_list(id)
    output = []
    for fileID in range(len(process[0]['files'])):
        extension = '.'.join(os.path.basename(process[0]['files'][fileID]).split('.')[1:])
        output.append({
            'fileID':fileID,
            'description':descriptions[extension] if extension in descriptions else "Unknown file",
            'display_name':os.path.relpath(
                process[0]['files'][fileID],
                process[0]['output']
            ),
            'url':'/api/v1/processes/%d/results/%d'%(id, fileID),
            'size':"%0.3f KB"%(
                os.path.getsize(os.path.join(
                    process[0]['output'],
                    process[0]['files'][fileID]
                ))/1024
            )
        })
    return output


def results_getfile(id, count = None, page = None, fileID = None):
    """Read data directly from a specific output file"""
    process = fetch_process(id)
    if not process[0]:
        return (
            {
                "code": 400,
                "message": "The requested process (%d) does not exist"%id,
                "fields": "id"
            },400
        )
    if is_running(id):
        return []
    gen_files_list(id)
    if fileID not in range(len(process[0]['files'])):
        return (
            {
                "code": 400,
                "message": "The requested fileID (%d) does not exist for this process (%d)" %(fileID, id),
                "fields": "fileID"
            },400
        )
    raw_reader = open(process[0]['files'][fileID])
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = [
        {column_filter(k):entry[k] for k in entry}
        for entry in itertools.islice(reader, (page-1)*count, page*count)
    ]
    raw_reader.close()
    return output


def results_getcols(id, fileID):
    """Get a mapping of standardized column names -> original column names"""
    process = fetch_process(id)
    if not process[0]:
        return (
            {
                "code": 400,
                "message": "The requested process (%d) does not exist"%id,
                "fields": "id"
            },400
        )
    if is_running(id):
        return {}
    gen_files_list(id)
    if fileID not in range(len(process[0]['files'])):
        return (
            {
                "code":400,
                "message":"The requested fileID (%d) does not exist for this process (%d)" %(fileID, id),
                "fields":"fileID"
            }, 400
        )
    raw_reader = open(process[0]['files'][fileID])
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = {column_filter(field):field for field in reader.fieldnames}
    raw_reader.close()
    return output


def staging(input, samplename, alleles, epitope_lengths, prediction_algorithms,
          peptide_sequence_length, gene_expn_file, transcript_expn_file,
          normal_snvs_coverage_file, normal_indels_coverage_file,
          tdna_snvs_coverage_file, tdna_indels_coverage_file,
          trna_snvs_coverage_file, trna_indels_coverage_file,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold, fasta_size, downstream_sequence_length, keep_tmp_files):
    """Stage input for a new pVAC-Seq run.  Generate a unique output directory and \
    save uploaded files to temporary locations (and give pVAC-Seq the filepaths). \
    Then forward the command to start()"""

    current_path = os.path.join(os.path.expanduser('~'), "Documents", "pVAC-Seq Output", samplename)
    if os.path.exists(current_path):
        i = 1
        while os.path.exists(current_path+"_"+str(i)):
            i += 1
        current_path += "_"+str(i)

    os.makedirs(os.path.join(current_path, 'Staging'), exist_ok=True)

    staged_input = open(os.path.join(current_path, "Staging", "input.vcf"), 'wb')
    input.save(staged_input)
    staged_input.flush()

    staged_additional_input_file_list = open(os.path.join(current_path, "Staging", "additional_input_file_list.yml"), 'w')

    staged_gene_expn_file = open(os.path.join(current_path, "Staging", "genes.fpkm_tracking"), 'wb')
    gene_expn_file.save(staged_gene_expn_file)
    staged_gene_expn_file.flush()
    if staged_gene_expn_file.tell():
        dump({"gene_expn_file": staged_gene_expn_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_transcript_expn_file = open(os.path.join(current_path, "Staging", "transcript.fpkm_tracking"), 'wb')
    transcript_expn_file.save(staged_transcript_expn_file)
    staged_transcript_expn_file.flush()
    if staged_transcript_expn_file.tell():
        dump({"transcript_expn_file" : staged_transcript_expn_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_normal_snvs_coverage_file = open(os.path.join(current_path, "Staging", "normal_snvs.bam_readcount"), 'wb')
    normal_snvs_coverage_file.save(staged_normal_snvs_coverage_file)
    staged_normal_snvs_coverage_file.flush()
    if staged_normal_snvs_coverage_file.tell():
        dump({"normal_snvs_coverage_file" : staged_normal_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_normal_indels_coverage_file = open(os.path.join(current_path, "Staging", "normal_indels.bam_readcount"), 'wb')
    normal_indels_coverage_file.save(staged_normal_indels_coverage_file)
    staged_normal_indels_coverage_file.flush()
    if staged_normal_indels_coverage_file.tell():
        dump({"normal_indels_coverage_file" : staged_normal_indels_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_tdna_snvs_coverage_file = open(os.path.join(current_path, "Staging", "tdna_snvs.bam_readcount"), 'wb')
    tdna_snvs_coverage_file.save(staged_tdna_snvs_coverage_file)
    staged_tdna_snvs_coverage_file.flush()
    if staged_tdna_snvs_coverage_file.tell():
        dump({"tdna_snvs_coverage_file" : staged_tdna_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_tdna_indels_coverage_file = open(os.path.join(current_path, "Staging", "tdna_indels.bam_readcount"), 'wb')
    tdna_indels_coverage_file.save(staged_tdna_indels_coverage_file)
    staged_tdna_indels_coverage_file.flush()
    if staged_tdna_indels_coverage_file.tell():
        dump({"tdna_indels_coverage_file" : staged_tdna_indels_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_trna_snvs_coverage_file = open(os.path.join(current_path, "Staging", "trna_snvs.bam_readcount"), 'wb')
    trna_snvs_coverage_file.save(staged_trna_snvs_coverage_file)
    staged_trna_snvs_coverage_file.flush()
    if staged_trna_snvs_coverage_file.tell():
        dump({"trna_snvs_coverage_file" : staged_trna_snvs_coverage_file.name}, staged_additional_input_file_list, default_flow_style=False)

    staged_trna_indels_coverage_file = open(os.path.join(current_path, "Staging", "trna_indels.bam_readcount"), 'wb')
    trna_indels_coverage_file.save(staged_trna_indels_coverage_file)
    staged_trna_indels_coverage_file.flush()
    if staged_trna_indels_coverage_file.tell():
        dump({"trna_indels_coverage_file" : staged_trna_indels_coverage_file.name}, staged_additional_input_file_list)
    staged_additional_input_file_list.flush()

    return start(staged_input.name, samplename, alleles, epitope_lengths, prediction_algorithms, current_path,
              peptide_sequence_length, staged_additional_input_file_list.name if staged_additional_input_file_list.tell() else "", # check if any data written to file
              net_chop_method, len(netmhc_stab), len(top_result_per_mutation), top_score_metric,
              binding_threshold, minimum_fold_change,
              normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
              expn_val, net_chop_threshold,
              fasta_size, downstream_sequence_length, len(keep_tmp_files))


def start(input, samplename, alleles, epitope_lengths, prediction_algorithms, output,
          peptide_sequence_length, additional_input_file_list,
          net_chop_method, netmhc_stab, top_result_per_mutation, top_score_metric,
          binding_threshold, minimum_fold_change,
          normal_cov, tdna_cov, trna_cov, normal_vaf, tdna_vaf, trna_vaf,
          expn_val, net_chop_threshold,
          fasta_size, downstream_sequence_length, keep_tmp_files):
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

    data['processid']+=1
    os.makedirs(os.path.dirname(logfile), exist_ok = True)
    children[data['processid']] = subprocess.Popen(
        command,
        stdout=open(logfile, 'w'),  # capture stdout in the logfile
        stderr=subprocess.STDOUT,
        # isolate the child in a new process group
        # this way it will remainin running no matter what happens to this process
        preexec_fn=os.setpgrp
    )
    # Store some data about the child process
    data['process-%d'%(data['processid'])] = {
        # Do the replacement so that the displayed command is actually valid
        # The executed command is automatically escaped as part of Popen
        'command': " ".join(command),
        'logfile':logfile,
        'pid':children[data['processid']].pid,
        'status': "Task Started",
        'output':os.path.abspath(output)
    }
    if 'reboot' not in data:
        data['reboot'] = reboot
    savedata()
    return data['processid']


def processes():
    """Returns a list of processes, and whether or not each process is running"""
    return [{
        'id':i,
        'running':is_running(i)
    } for i in range(data['processid']+1) if 'process-%d'%i in data]


def process_info(id):
    """Returns more detailed information about a specific process"""
    process = fetch_process(id)
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    reader = open(process[0]['logfile'])
    log = spinner.sub('', reader.read()).strip()
    process[0]['status'] = log.split("\n")[-1]
    reader.close()
    if not is_running(id):
        if process[1]:
            process[0]['status'] = "Process Complete: %d"%process[1].returncode
        else:
            process[0]['status'] = "Process Complete"
        # If there is a staging directory, remove it
        if os.path.isdir(os.path.join(process[0]['output'], 'Staging')):
            shutil.rmtree(os.path.join(process[0]['output'], 'Staging'))
    savedata()
    return {
        'pid':process[0]['pid'],
        'id':id,
        'attached': bool(process[1]),
        'command':process[0]['command'],
        'status':process[0]['status'],
        'log':log,
        'output':process[0]['output'],
        'running':is_running(id)
    }


def stop(id):
    """Stops the requested process.  This is only allowed if the child is still attached"""
    status = process_info(id)
    if type(status) == dict:  # status could be an error object if the id is invalid
        if status['running'] and status['pid']>1:
            children[id].terminate()
    return status


def shutdown():
    """Stops all attached, running children"""
    output = []
    for i in range(data['processid']+1):
        if is_running(i) and i in children:
            output.append(i)
            children[i].wait(.1)
            if is_running(i):
                children[i].terminate()
    return output


def test():
    """Return the submission page (a stand-in until there is a proper ui for submission)"""
    reader = open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_start.html'))
    data = reader.read()
    reader.close()
    return data


def check_allele(allele):
    """Checks if the requested allele is supported by pVAC-Seq or not"""
    global allele_file
    if not allele_file:
        allele_file = tempfile.TemporaryFile('w+')
        subprocess.call(['pvacseq', 'valid_alleles'], stdout=allele_file)
    allele_file.seek(0)
    for line in allele_file:
        if line.strip() == allele:
            return True
    return False


def reset(clearall):
    """Clears out finished processes from the record"""
    output = []
    for i in range(data['processid']+1):
        if 'process-%d'%i in data and is_running(i):
            if i not in children:
                shutil.rmtree(data['process-%d'%i]['output'])
                del data['process-%d'%i]
                output.append(i)
            elif clearall:
                shutil.rmtree(data['process-%d'%i]['output'])
                del data['process-%d'%i]
                del children[i]
                output.append(i)
    # Set the processid to the highest child process from this session
    data['processid'] = max(i for i in range(data['processid']+1) if 'process-%d'%i in data)
    if clearall and 'reboot' in data:
        del data['reboot']
    savedata()
    return output
