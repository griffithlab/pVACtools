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
from flask import current_app
import watchdog.events

spinner = re.compile(r'[\\\b\-/|]{2,}')
queryfilters = re.compile(r'(.+)(<=?|>=?|!=|==)(.+)')
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

def initialize():
    """Setup anything that needs to be configured in the app, but within\
    the context of the controllers"""
    if not current_app.config['initialized']:
        current_app.config['initialized'] = True
        db = current_app.config['db']

        if 'filtertables' not in data:
            data['filtertables']={}

        #We need to setup listeners on the dropbox so that
        #fileID's are consistent between runs, and between changes to the directory
        #ie: file 'b' with id 0.  UI requests listing, gets {0: 'b'}
        #User adds file 'a'
        #UI requests id 0, expecting 'b', but gets 'a' since the listing has changed
        #The event listeners remove that race condition
        dbr = current_app.config['dropbox_dir']
        watcher = current_app.config['watcher']
        if 'dropbox' not in data:
            data['dropbox'] = {}
        fileID = 0
        current = set(os.listdir(current_app.config['dropbox_dir']))
        recorded = set(data['dropbox'].values())
        targets = {k for k in data['dropbox'] if data['dropbox'][k] in recorded-current}
        for fileID in targets:
            del data['dropbox'][fileID]
        for filename in current-recorded:
            while str(fileID) in data['dropbox']:
                fileID += 1
            print("Assigning file:", fileID,"-->",filename)
            data['dropbox'][str(fileID)] = filename
        savedata()

        def _create(event):
            filename = os.path.relpath(
                event.src_path,
                dbr
            )
            fileID = 0
            while str(fileID) in data['dropbox']:
                fileID += 1
            print("Creating file:", fileID, "-->",filename)
            data['dropbox'][str(fileID)] = filename
            savedata()
        watcher.subscribe(
            _create,
            watchdog.events.FileCreatedEvent
        )

        def _delete(event):
            filename = os.path.relpath(
                event.src_path,
                dbr
            )
            for key in data['dropbox']:
                if data['dropbox'][key] == filename:
                    del data['dropbox'][key]
                    print("Deleting file:",key,'-->', filename)
                    query = db.prepare("SELECT 1 FROM information_schema.tables WHERE table_name = $1")
                    if len(query('data_dropbox_'+str(key))):
                        db.execute("DROP TABLE data_dropbox_"+str(key))
                    savedata()
                    return
        watcher.subscribe(
            _delete,
            watchdog.events.FileDeletedEvent
        )

        def _move(event):
            filesrc = os.path.relpath(
                event.src_path,
                dbr
            )
            filedest = os.path.relpath(
                event.dest_path,
                dbr
            )
            for key in data['dropbox']:
                if data['dropbox'][key] == filesrc:
                    data['dropbox'][key] = filedest
                    print("Moving file:", key,'(',filesrc,'-->',filedest,')')
                    savedata()
                    return
        watcher.subscribe(
            _move,
            watchdog.events.FileMovedEvent
        )

# def assignFileId(filename):
#     if filename not in data['dropbox'].values():
#         watcher = current_app.config['watcher']
#         mover = watcher.subscribe(
#             lambda x:data['dropbox'][fileID]=x.dest_path,
#             watchdog.events.FileMovedEvent
#         )
#         def _deleter(event):
#             mover()
#             del data['dropbox'][fileID]
#         deleter = watcher.subscribe(
#             lambda x:data[]
#         )

# a mapping to provide a description of each result file based on its extension
descriptions = {
    'json':"Metadata regarding a specific run of pVAC-Seq",
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
    return column.replace(' ', '_').replace('-', '_').lower().strip()

def gen_files_list(id):
    """Generate the list of result files for a given process.  Stash them for later use"""
    #  if 'process-%d'%id not in data:
    #     raise KeyError("The requested process (%d) does not exist"%id)
    if 'files' not in data['process-%d'%id]:
        data['process-%d'%id]['files'] = []
        base_dir = data['process-%d'%id]['output']
        if os.path.isdir(os.path.join(base_dir, 'MHC_Class_I')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'MHC_Class_I'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'MHC_Class_I', path)):
                    data['process-%d'%id]['files'].append(os.path.join(base_dir, 'MHC_Class_I', path))
        if os.path.isdir(os.path.join(base_dir, 'MHC_Class_II')):
            for path in sorted(os.listdir(os.path.join(base_dir, 'MHC_Class_II'))):
                if path.endswith('.tsv') and os.path.isfile(os.path.join(base_dir, 'MHC_Class_II', path)):
                    data['process-%d'%id]['files'].append(os.path.join(base_dir, 'MHC_Class_II', path))
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
    initialize()
    if id == -1:
        return list_dropbox()
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


def results_getfile(id, fileID, count, page, filters, sort, direction):
    """(deprecated) Read data directly from a specific output file"""
    return filterfile(
        parentID = id,
        fileID = fileID,
        count = count,
        page = page,
        filters = filters,
        sort = sort,
        direction = direction
    )
    # initialize()
    # process = fetch_process(id)
    # if not process[0]:
    #     return (
    #         {
    #             "code": 400,
    #             "message": "The requested process (%d) does not exist"%id,
    #             "fields": "id"
    #         },400
    #     )
    # if is_running(id):
    #     return []
    # gen_files_list(id)
    # if fileID not in range(len(process[0]['files'])):
    #     return (
    #         {
    #             "code": 400,
    #             "message": "The requested fileID (%d) does not exist for this process (%d)" %(fileID, id),
    #             "fields": "fileID"
    #         },400
    #     )
    # raw_reader = open(process[0]['files'][fileID])
    # reader = csv.DictReader(raw_reader, delimiter='\t')
    # output = [
    #     {column_filter(k):entry[k] for k in entry}
    #     for entry in itertools.islice(reader, (page-1)*count, page*count)
    # ]
    # raw_reader.close()
    # return output


def results_getcols(id, fileID):
    """Get a mapping of standardized column names -> original column names"""
    initialize()
    if id==-1:
        if str(fileID) not in data['dropbox']:
            return {
                'code':400,
                'message':'The requested file (%d) does not exist'%fileID,
                'fields':'fileID'
            }
        raw_reader = open(
            os.path.join(
                os.path.abspath(current_app.config['dropbox_dir']),
                data['dropbox'][str(fileID)]
            )
        )
        reader = csv.DictReader(raw_reader, delimiter='\t')
        output = {
            column_filter(field):field for field in reader.fieldnames
        }
        raw_reader.close()
        return output
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
    initialize()

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

    #Regexes to build the config object:
    #>>(\w+)\|(.+)\|(.+) -> if $1!=$2:\n\t\tconfigObj['$3']=$1
    #if (\w+)!=False -> if $1
    #configObj\['CST/(.+)/(.+)'\]=.* -> configObj['$2']=$1
    configObj = {
        'action':'run',
        'input_file': input,
        'sample_name':samplename,
        'alleles':','.split(alleles),
        'prediction_algorithms':','.split(prediction_algorithms),
        'output_directory':output
    }
    if epitope_lengths!=10:
        configObj['epitope_lengths']=','.split(epitope_lengths)
    if peptide_sequence_length!=21:
        configObj['peptide_sequence_length']=peptide_sequence_length
    if additional_input_file_list!='':
        configObj['additional_input_files']=','.split(additional_input_file_list)
    if net_chop_method!='':
        configObj['net_chop_method']=net_chop_method
    if netmhc_stab:
        configObj['netmhc_stab']=True
    if top_result_per_mutation:
        configObj['top_result_per_mutation']=True
    if top_score_metric!='median':
        configObj['top_score_metric']=top_score_metric
    if binding_threshold!=500:
        configObj['binding_threshold']=binding_threshold
    if minimum_fold_change!=0:
        configObj['minimum_fold_change']=minimum_fold_change
    if normal_cov!=5:
        configObj['normal_coverage_cutoff']=normal_cov
    if tdna_cov!=10:
        configObj['tumor_dna_coverage_cutoff']=tdna_cov
    if trna_cov!=10:
        configObj['tumor_rna_coverage_cutoff']=trna_cov
    if normal_vaf!=2:
        configObj['normal_vaf_cutoff']=normal_vaf
    if tdna_vaf!=40:
        configObj['tumor_dna_vaf_cutoff']=tdna_vaf
    if trna_vaf!=40:
        configObj['tumor_rna_vaf_cutoff']=trna_vaf
    if expn_val!=1:
        configObj['expression_cutoff']=expn_val
    if net_chop_threshold!=0.5:
        configObj['netchop_threshold']=net_chop_threshold
    if fasta_size!=200:
        configObj['fasta_size']=fasta_size
    if downstream_sequence_length!="1000":
        configObj['downstream_sequence_length']=downstream_sequence_length

    writer = open(os.path.join(
        os.path.abspath(output),
        'config.json'
    ),'w')
    json.dump(configObj)
    writer.close()
    return data['processid']


def processes():
    """Returns a list of processes, and whether or not each process is running"""
    initialize()
    return [{
        'id':i,
        'running':is_running(i)
    } for i in range(data['processid']+1) if 'process-%d'%i in data]


def process_info(id):
    """Returns more detailed information about a specific process"""
    initialize()
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
    initialize()
    status = process_info(id)
    if type(status) == dict:  # status could be an error object if the id is invalid
        if status['running'] and status['pid']>1:
            children[id].terminate()
    return status


def shutdown():
    """Stops all attached, running children"""
    initialize()
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
    initialize()
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
    initialize()
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

def list_dropbox():
    initialize()
    return [
        {
            'fileID':key,
            'description':descriptions[ext] if ext in descriptions else "Unknown file",
            'display_name':os.path.relpath(
                os.path.join(current_app.config['dropbox_dir'], data['dropbox'][key])
            ),
            'url':'/api/v1/dropbox/files/%d'%(int(key)),
            'size':"%0.3f KB"%(
                os.path.getsize(os.path.join(
                    current_app.config['dropbox_dir'],
                    data['dropbox'][key]
                ))/1024
            )
        } for (key, ext) in zip(
            data['dropbox'].keys(),
            ('.'.join(data['dropbox'][key].split('.')[1:]) for key in data['dropbox'])
        )
    ]

def init_column_mapping(row):
    mapping = {column_filter(col):str for col in row}
    defs = {column_filter(col):'text' for col in row}
    for (col, val) in row.items():
        col = column_filter(col)
        if mapping[col] != float:
            try:
                int(val)
                if mapping[col] == str:
                    print("Assigning int to",col,"based on",val)
                    mapping[col]=int
                    defs[col] = 'integer'
            except ValueError:
                try:
                    float(val)
                    print("Assigning float to",col,"based on",val)
                    mapping[col]=float
                    defs[col] = 'decimal'
                except ValueError:
                    pass
    defs['start'] = 'bigint'
    defs['stop'] = 'bigint'
    return (mapping, defs)

def column_mapping(row, mapping):
    #we have to read the whole file, in case there's just a NA
    #in a normally numerical field on the first row
    #ALTER TABLE ? ALTER COLUMN ? SET DATA TYPE ? USING null
    output = {}
    changes = {}
    for (col, val) in row.items():
        col = column_filter(col)
        if mapping[col] == str:
            try:
                int(val)
                print("Assigning int to",col,"based on",val)
                mapping[col]=int
                changes[col] = int
            except ValueError:
                try:
                    float(val)
                    print("Assigning float to",col,"based on",val)
                    mapping[col]=float
                    changes[col] = float
                except ValueError:
                    pass
        try:
            output[col] = mapping[col](val)
        except ValueError:
            output[col] = None
    return (mapping, output, changes)

def filterfile(parentID, fileID, count, page, filters, sort, direction):
    """Gets the file ID belonging to the parent.\
    For result files, the parentID is the process ID that spawned them.\
    For dropbox files, the parentID is -1"""
    print(parentID)
    initialize()

    #first, generate the key
    tablekey = "data_%s_%s"%(
        (parentID if parentID >=0 else 'dropbox'),
        fileID
    )

    #check if the table exists:
    db = current_app.config['db']
    query = db.prepare("SELECT 1 FROM information_schema.tables WHERE table_name = $1")
    if not len(query(tablekey)): #table does not exist
        #Open a reader to cache the file in the database
        if parentID != -1:
            process = fetch_process(parentID)
            if not process[0]:
                return (
                    {
                        "code": 400,
                        "message": "The requested process (%d) does not exist"%parentID,
                        "fields": "parentID"
                    },400
                )
            if is_running(parentID):
                return []
            gen_files_list(parentID)
            if fileID not in range(len(process[0]['files'])):
                return (
                    {
                        "code": 400,
                        "message": "The requested fileID (%d) does not exist for this process (%d)" %(fileID, parentID),
                        "fields": "fileID"
                    },400
                )
            raw_reader = open(process[0]['files'][fileID])
        else:
            if str(fileID) not in data['dropbox']:
                return (
                    {
                        "code": 400,
                        "message": "The requested fileID (%d) does not exist in the dropbox"%fileID,
                        "fields":"fileID"
                    }
                )
            raw_reader = open(os.path.join(
                os.path.abspath(current_app.config['dropbox_dir']),
                data['dropbox'][str(fileID)]
            ))
        reader = csv.DictReader(raw_reader, delimiter='\t')

        tmp_reader = open(raw_reader.name)
        tmp = csv.DictReader(tmp_reader, delimiter='\t')
        init = next(tmp)
        tmp_reader.close()

        #Get an initial estimate of column datatypes from the first row
        (mapping, column_names) = init_column_mapping(init)
        tablecolumns = "\n".join( #use the estimated types to create the table
            "%s %s,"%(colname, column_names[colname])
            for colname in column_names
        )[:-1]
        CREATE_TABLE = "CREATE TABLE %s (\
            rowid SERIAL PRIMARY KEY NOT NULL,\
            %s\
        )"%(tablekey, tablecolumns)
        db.execute(CREATE_TABLE)
        #mark the table for deletion when the server shuts down
        if 'db-clean' not in current_app.config:
            current_app.config['db-clean'] = [tablekey]
        else:
            current_app.config['db-clean'].append(tablekey)
        #prepare the insertion query
        insert = db.prepare("INSERT INTO %s (%s) VALUES (%s)" %(
            tablekey,
            ','.join(column_names),
            ','.join('$%d'%i for (_,i) in zip(column_names, range(1,sys.maxsize)))
        ))
        update = "ALTER TABLE %s "%tablekey
        for row in reader:
            #process each row
            #We format the data in the row and update column data types, if necessary
            (mapping, formatted, changes) = column_mapping(row, mapping)
            alter_cols = []
            for (k,v) in changes.items():
                #if there were any changes to the data type, update the table
                #since we only ever update a text column to int/decimal, then
                #it's okay to nullify the data
                typ = ''
                if v == int:
                    typ = 'bigint' if k in {'start', 'stop'} else 'integer'
                elif v == float:
                    typ = 'decimal'
                alter_cols.append(
                    "ALTER COLUMN %s SET DATA TYPE %s USING null"%(
                        k,
                        typ
                    )
                )
            if len(changes):
                #Re-generate the insert statement since the data types changed
                print("Alter:",update+','.join(alter_cols))
                db.execute(update+','.join(alter_cols))
                insert = db.prepare("INSERT INTO %s (%s) VALUES (%s)" %(
                    tablekey,
                    ','.join(column_names),
                    ','.join('$%d'%i for (_,i) in zip(column_names, range(1,sys.maxsize)))
                ))
            #insert the row
            insert(*[formatted[column] for column in column_names])
        raw_reader.close()
    typequery = db.prepare("SELECT column_name, data_type FROM information_schema.columns WHERE table_name = $1")
    column_defs = typequery(tablekey)
    column_maps = {}
    for (col, typ) in column_defs:
        if 'int' in typ:
            column_maps[col] = int
        elif typ == 'numeric'or typ == 'decimal':
            column_maps[col] = float
        else:
            column_maps[col] = str
    formatted_filters = []
    for i in range(len(filters)):
        f = filters[i].strip()
        if not len(f):
            continue
        result = queryfilters.match(f)
        if not result:
            return {
                "code":400,
                "message":"Encountered an invalid filter (%s)"%f,
                "fields":"filters"
            }
        colname = column_filter(result.group(1))
        if colname not in column_maps:
            return {
                "code":400,
                "message":"Unknown column name %s"%result.group(1),
                "fields":"filters"
            }
        op = result.group(2)
        typ = column_maps[colname]
        val = None
        try:
            val = column_maps[colname](result.group(3))
        except ValueError:
            return {
                "code":400,
                "message":"Value %s cannot be formatted to match the type of column %s (%s)"%(
                    result.group(3),
                    result.group(1),
                    typ
                )
            }
        if typ == str and (op in {'==', '!='}):
            formatted_filters.append(
                json.dumps(colname) + (' not ' if '!' in op else ' ') + "LIKE '%s'"%(
                    json.dumps(val)[1:-1]
                )
            )
        else: #type is numerical
            op = op.replace('==', '=')
            formatted_filters.append(
                '%s %s %s'%(
                    json.dumps(colname),
                    op,
                    json.dumps(val)
                )
            )
    raw_query = "SELECT %s FROM %s"%(
        ','.join([k[0] for k in column_defs]),
        tablekey
    )
    if len(formatted_filters):
        raw_query += " WHERE "+" AND ".join(formatted_filters)
    if sort:
        if column_filter(sort) not in column_maps:
            return {
                'code':400,
                'message':'Invalid column name %s'%sort,
                'fields':'sort'
            }
        raw_query += " ORDER BY %s"%(column_filter(sort))
        if direction:
            raw_query += " "+direction
    if count:
        raw_query += " LIMIT %d"%count
    if page:
        raw_query += " OFFSET %d"%(page*count)
    print("Query:",raw_query)
    query = db.prepare(raw_query)
    import decimal
    decimalizer = lambda x:(float(x) if type(x) == decimal.Decimal else x)
    return [
        {
            colname:decimalizer(value) for (colname, value) in zip(
                [k[0] for k in column_defs],
                [val for val in row]
            )
        } for row in query.rows()
    ]
