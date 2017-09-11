import os
import re
import shutil
from flask import current_app, redirect
import json
import sys
import subprocess
from shlex import split
from .utils import descriptions, filterdata
from shutil import move as movetree

spinner = re.compile(r'[\\\b\-/|]{2,}')

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

def fixpath(src, root, *keys):
    if type(src) == 'str':
        return os.path.relpath(
            src,
            root
        )
    for key in keys:
        src[key] = os.path.relpath(
            src[key],
            root,
        )
    return src

def processes(filters, sorting, page, count):
    """Returns a list of processes, and whether or not each process is running"""
    data = current_app.config['storage']['loader']()
    #Python comprehensions are great!
    # data['process-%d'%id]['returncode'] = process[1].returncode
    # data['process-%d'%id]['status'] = 1 if process[1].returncode == 0 else -1
    return filterdata([
         {
             'id':proc[0],
             'running':is_running(proc[1]),
             'returncode':proc[1][1].returncode if (not is_running(proc[1])) and proc[1][1] else (
                 proc[1][0]['returncode'] if 'returncode' in proc[1][0] else 0
             ),
             'status':0 if is_running(proc[1]) else (
                -1 if (
                    (proc[1][1] and proc[1][1].returncode != 0) or
                    ('returncode' in proc[1][0] and proc[1][0]['returncode'] != 0)
                ) else 1
             ),
             'url':'/api/v1/processes/%d'%proc[0],
             'results_url':'/api/v1/processes/%d/results'%proc[0],
             'attached':bool(proc[1][1]),
             'output':os.path.relpath(
                 proc[1][0]['output'],
                 current_app.config['files']['data-dir']
             ),
             'pid':proc[1][0]['pid'],
             'command':proc[1][0]['command'],
             'files':[
                 {
                     'fileID':fileID,
                     'url':'/api/v1/processes/%d/results/%s'%(
                         proc[0],
                         fileID
                     ),
                     'display_name':filedata['display_name'],
                     'description':filedata['description']
                 }
                 for (fileID, filedata) in data['process-%d'%proc[0]]['files'].items()
             ],
             'parameters':(
                 fixpath(json.load(open(os.path.join(
                     proc[1][0]['output'],
                     'config.json'
                 ))), current_app.config['files']['data-dir'], 'output')
                 if os.path.isfile(os.path.join(
                     proc[1][0]['output'],
                     'config.json'
                 )) else {}
             )
         } for proc in
            map(
                lambda x: (x,fetch_process(
                    x,
                    data,
                    current_app.config['storage']['children']
                )),
                range(data['processid']+1)
            ) if 'process-%d'%(proc[0]) in data
    ], filters, sorting, page, count)


def process_info(id):
    """Returns more detailed information about a specific process"""
    data = current_app.config['storage']['loader']()
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    reader = open(process[0]['logfile'])
    log = spinner.sub('', reader.read()).strip().split(os.linesep)
    process[0]['last_message'] = log[-1]
    reader.close()
    if not is_running(process):
        if process[1]:
            data['process-%d'%id]['returncode'] = process[1].returncode
            data['process-%d'%id]['status'] = 1 if process[1].returncode == 0 else -1
        # If there is a staging directory, remove it
        if os.path.isdir(os.path.join(process[0]['output'], 'Staging')):
            shutil.rmtree(os.path.join(process[0]['output'], 'Staging'))
    data.save()
    configfile = os.path.join(
        process[0]['output'],
        'config.json'
    )
    params = {}
    if os.path.isfile(configfile):
        params = json.load(open(configfile))
        params['output'] = os.path.relpath(
            params['output'],
            current_app.config['files']['data-dir']
        )
    return {
        'pid':process[0]['pid'],#
        'id':id,#
        'results_url':'/api/v1/processes/%d/results'%id,#
        'attached': bool(process[1]),#
        'command':process[0]['command'],#
        'returncode':process[0]['returncode'] if 'returncode' in process[0] else 0,
        'status':0 if is_running(process) else (process[0]['status'] if 'status' in process[0] else 0),
        'last_message':process[0]['last_message'],
        'log':log,
        'log_updated_at':(
            int(os.stat(process[0]['logfile']).st_mtime)
            if os.path.isfile(process[0]['logfile'])
            else 0
        ),
        'output':os.path.relpath(
            process[0]['output'],
            current_app.config['files']['data-dir']
        ),#
        'running':is_running(process),#
        'files':[
            {
                'fileID':fileID,
                'url':'/api/v1/processes/%d/results/%s'%(
                    id,
                    fileID
                ),
                'display_name':filedata['display_name'],
                'description':filedata['description']
            }
            for (fileID, filedata) in data['process-%d'%id]['files'].items()
        ],
        'parameters':params
    }


def stop(id):
    """Stops the requested process.  This is only allowed if the child is still attached"""
    data = current_app.config['storage']['loader']()
    status = process_info(id)
    if type(status) == dict:  # status could be an error object if the id is invalid
        if status['running'] and status['pid']>1:
            current_app.config['storage']['children'][id].terminate()
    return status


def shutdown():
    """Stops all attached, running children"""
    data = current_app.config['storage']['loader']()
    output = []
    for i in range(data['processid']+1):
        proc = fetch_process(i, data, current_app.config['storage']['children'])
        if is_running(proc) and i in current_app.config['storage']['children']:
            output.append(i)
            try:
                current_app.config['storage']['children'][i].wait(.1)
            except subprocess.TimeoutExpired:
                current_app.config['storage']['children'][i].terminate()
    return output


def reset(clearall):
    """Clears out and archives finished processes from the record"""
    with current_app.config['storage']['synchronizer']:
        data = current_app.config['storage']['loader']()
        output = []
        for i in range(data['processid']+1):
            proc = fetch_process(i, data, current_app.config['storage']['children'])
            if 'process-%d'%i in data and not is_running(proc):
                if i not in current_app.config['storage']['children']:
                    result = archive(i)
                    if type(result) == tuple:
                        return result
                    output.append(i)
                elif clearall:
                    stop(i)
                    result = archive(i)
                    if type(result) == tuple:
                        return result
                    output.append(i)
        if clearall and 'reboot' in data:
            del data['reboot']
        data.save()
        return output

def archive(processID):
    """Archives the results from the given process"""
    with current_app.config['storage']['synchronizer']:
        data = current_app.config['storage']['loader']()
        if 'process-%d'%processID not in data:
            return (
                {
                    'code':400,
                    'message': "The requested process (%d) does not exist"%processID,
                    'fields':"processID"
                },
                400
            )
        proc = fetch_process(processID, data, current_app.config['storage']['children'])
        if is_running(proc):
            return (
                {
                    'code':400,
                    'message': "The requested process (%d) is still running.\
                    Stop the process or wait for it to finish before archiving"%processID,
                    'fields':"processID"
                },
                400
            )
        dirname = os.path.basename(proc[0]['output'])
        archive = os.path.join(
            current_app.config['files']['data-dir'],
            'archive'
        )
        destpath = os.path.join(archive, dirname)
        if os.path.exists:
            i = 1
            while os.path.exists(destpath+'_%d'%i):
                i+=1
            destpath += '_%d'%i
        movetree(
            proc[0]['output'],
            destpath
        )
        del data['process-%d'%processID]
        if processID in current_app.config['storage']['children']:
            del current_app.config['storage']['children'][processID]
        # Set the processid to the highest child process from this session
        data['processid'] = max([0]+[i for i in range(data['processid']+1) if 'process-%d'%i in data])
        data.save()
        return "OK"

def restart(processID):
    """Restart a previously started process"""
    data = current_app.config['storage']['loader']()
    key = 'process-%d' % processID
    if key not in data:
        return (
            {
                'code':400,
                'message':"Invalid ProcessID: %d"%processID,
                'fields':"processID"
            },
            400
        )
    proc = fetch_process(processID, data, current_app.config['storage']['children'])
    if is_running(proc):
        return (
            {
                'code':400,
                'message':"The requested process (%d) is still running"%processID,
                'fields':"processID"
            },
            400
        )
    with current_app.config['storage']['synchronizer']:
        data = current_app.config['storage']['loader']() #refresh the data
        logfile = data[key]['logfile']
        os.makedirs(os.path.dirname(logfile), exist_ok = True)
        current_app.config['storage']['children'][processID] = subprocess.Popen(
            split(data[key]['command']),
            stdout=open(logfile, 'w'),  # capture stdout in the logfile
            stderr=subprocess.STDOUT,
            # isolate the child in a new process group
            # this way it will remainin running no matter what happens to this process
            preexec_fn=os.setpgrp
        )
        # Store some data about the child process
        data[key]['files']={}
        data[key]['status']=0
        data[key]['pid']=current_app.config['storage']['children'][processID].pid
        if 'reboot' not in data:
            data.addKey(
                'reboot',
                current_app.config['reboot'],
                current_app.config['files']['processes']
            )
        data.save()
        return redirect("/api/v1/processes/%d"%processID, 302)
