import os
import re
import shutil
from flask import current_app
from .utils import initialize, savedata

spinner = re.compile(r'[\\\b\-/|]{2,}')

def gen_files_list(id, data):
    """Generate the list of result files for a given process.  Stash them for later use"""
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
        savedata(data)
    return data

def fetch_process(id,data,children):
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

def is_running(process):
    """Returns True if the requested process looks like it's still running"""
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


def processes():
    """Returns a list of processes, and whether or not each process is running"""
    data = initialize()
    return [{
        'id':i,
        'running':is_running(fetch_process(i, data, current_app.config['children']))
    } for i in range(data['processid']+1) if 'process-%d'%i in data]


def process_info(id):
    """Returns more detailed information about a specific process"""
    data = initialize()
    process = fetch_process(id, data, current_app.config['children'])
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
    if not is_running(process):
        if process[1]:
            process[0]['status'] = "Process Complete: %d"%process[1].returncode
        else:
            process[0]['status'] = "Process Complete"
        # If there is a staging directory, remove it
        if os.path.isdir(os.path.join(process[0]['output'], 'Staging')):
            shutil.rmtree(os.path.join(process[0]['output'], 'Staging'))
    savedata(data)
    return {
        'pid':process[0]['pid'],
        'id':id,
        'attached': bool(process[1]),
        'command':process[0]['command'],
        'status':process[0]['status'],
        'log':log,
        'output':process[0]['output'],
        'running':is_running(process)
    }


def stop(id):
    """Stops the requested process.  This is only allowed if the child is still attached"""
    data = initialize()
    status = process_info(id)
    if type(status) == dict:  # status could be an error object if the id is invalid
        if status['running'] and status['pid']>1:
            current_app.config['children'][id].terminate()
    return status


def shutdown():
    """Stops all attached, running children"""
    data = initialize()
    output = []
    for i in range(data['processid']+1):
        proc = fetch_process(i, data, current_app.config['children'])
        if is_running(proc) and i in current_app.config['children']:
            output.append(i)
            current_app.config['children'][i].wait(.1)
            if is_running(proc):
                current_app.config['children'][i].terminate()
    return output


def reset(clearall):
    """Clears out finished processes from the record"""
    data = initialize()
    output = []
    for i in range(data['processid']+1):
        proc = fetch_process(i, data, current_app.config['children'])
        if 'process-%d'%i in data and is_running(proc):
            if i not in current_app.config['children']:
                shutil.rmtree(data['process-%d'%i]['output'])
                del data['process-%d'%i]
                output.append(i)
            elif clearall:
                shutil.rmtree(data['process-%d'%i]['output'])
                del data['process-%d'%i]
                del current_app.config['children'][i]
                output.append(i)
    # Set the processid to the highest child process from this session
    data['processid'] = max(i for i in range(data['processid']+1) if 'process-%d'%i in data)
    if clearall and 'reboot' in data:
        del data['reboot']
    savedata(data)
    return output
