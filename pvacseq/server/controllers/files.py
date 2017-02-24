import os
import csv
from flask import current_app
import subprocess
from .processes import fetch_process, is_running, gen_files_list
from .database import filterfile
from .utils import initialize, savedata, descriptions, column_filter


def results_get(id):
    """Get the list of result files from a specific pVAC-Seq run"""
    data = initialize()
    if id == -1:
        return list_dropbox()
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    if is_running(process):
        return []
    data = gen_files_list(id, data)
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
            ),
            'rows':int(subprocess.check_output(['wc', '-l', os.path.join(
                process[0]['output'],
                process[0]['files'][fileID]
            )]).decode().split()[0])-1,
        })
    return output


def results_getfile(id, fileID, count, page, filters, sort, direction):
    """(DEPRECATED) Read data directly from a specific output file"""
    return filterfile(
        parentID = id,
        fileID = fileID,
        count = count,
        page = page,
        filters = filters,
        sort = sort,
        direction = direction
    )



def results_getcols(id, fileID):
    """Get a mapping of standardized column names -> original column names"""
    data = initialize()
    if id==-1:
        if str(fileID) not in data['dropbox']:
            return {
                'code':400,
                'message':'The requested file (%d) does not exist'%fileID,
                'fields':'fileID'
            }
        raw_reader = open(
            os.path.join(
                os.path.abspath(current_app.config['files']['dropbox-dir']),
                data['dropbox'][str(fileID)]
            )
        )
        reader = csv.DictReader(raw_reader, delimiter='\t')
        output = {
            column_filter(field):field for field in reader.fieldnames
        }
        raw_reader.close()
        return output
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code": 400,
                "message": "The requested process (%d) does not exist"%id,
                "fields": "id"
            },400
        )
    if is_running(process):
        return {}
    data = gen_files_list(id, data)
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

def list_dropbox():
    data = initialize()
    return [
        {
            'fileID':key,
            'description':descriptions[ext] if ext in descriptions else "Unknown file",
            'display_name':os.path.relpath(
                os.path.join(current_app.config['files']['dropbox-dir'], data['dropbox'][key])
            ),
            'url':'/api/v1/processes/-1/results/%d'%(int(key)),
            'size':"%0.3f KB"%(
                os.path.getsize(os.path.join(
                    current_app.config['files']['dropbox-dir'],
                    data['dropbox'][key]
                ))/1024
            ),
            'rows':int(subprocess.check_output(['wc', '-l', os.path.join(
                current_app.config['files']['dropbox-dir'],
                data['dropbox'][key]
            )]).decode().split()[0])-1,
        } for (key, ext) in zip(
            data['dropbox'].keys(),
            ('.'.join(data['dropbox'][key].split('.')[1:]) for key in data['dropbox'])
        )
    ]
