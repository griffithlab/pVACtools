import os
import csv
import re
from flask import current_app
import subprocess
from .processes import fetch_process, is_running
from .database import filterfile
from .utils import descriptions, is_visualizable, visualization_type, column_filter, filterdata, sort, fullresponse, nav_to_dir

# details for each file to be appended to the output of results_get
def resultfile(id, process, file):
    fileID = file['fileID']
    return({
        'fileID':int(fileID),
        'description':file['description'],
        'display_name':file['display_name'],
        'is_visualizable': file['is_visualizable'],
        'visualization_type': file['visualization_type'],
        'type':'file',
        'url':'/api/v1/processes/%d/results/%s'%(id, fileID),
        'size':os.path.getsize(process[0]['files'][fileID]['fullname']),
        'rows':int(subprocess.check_output([
            'wc',
            '-l',
            process[0]['files'][fileID]['fullname']
        ]).decode().split()[0])-1,
    })

def prelim_res(id, type, process, data):
    output = []
    for file in data:
        if file['type'] == 'file':
            if "default" in type:
                if not (re.search("/tmp/", file['display_name'])):
                    output.append(resultfile(id, process, file))
            elif "all" in type:
                output.append(resultfile(id, process, file))
            else:
                for filter in type:
                    if (re.search('%s.tsv'%filter, file['display_name'])):
                        output.append(resultfile(id, process, file))
        elif file['type'] == 'directory':
            output.append({
                'display_name':file['display_name'],
                'type':'directory',
                'contents':prelim_res(id, type, process, file['contents'])
            })
    return output

def results_get(id, type, filters, sorting, page, count):
    """Get the list of result files from a specific pVAC-Seq run"""
    data = current_app.config['storage']['loader']()
    if id == -1:
        return list_visualize(type, filters, sorting, page, count)
    process = fetch_process(id, data, current_app.config['storage']['children'])
    if not process[0]:
        return (
            {
                "code":400,
                "message":"The requested process (%d) does not exist"%id,
                "fields":"id"
            },400
        )
    res_dir = os.path.join(current_app.config['files']['data-dir'], 'results')
    #nav_to_dir expects a file for arg1, so add a fake file to end of path
    tree = nav_to_dir(
        os.path.join(process[0]['output'],'fake_file'), res_dir, current_app.config['storage']['manifest']['results']
    )
    output = prelim_res(id, type, process, tree)
    return filterdata(output, filters, sorting, page, count)

def inputfile(file):
    return({
        'fileID':int(file['fileID']),
        'description':file['description'],
        'display_name':file['display_name'],
        'is_input': file['display_name'].endswith('.vcf') or file['display_name'].endswith('.vcf.gz'),
        'is_visualizable': file['is_visualizable'],
        'visualization_type': file['visualization_type'],
        'type':'file',
    })

def prelim_inp(data):
    output = []
    for file in data:
        if file['type'] == 'file':
            output.append(inputfile(file))
        elif file['type'] == 'directory':
            output.append({
                'display_name':file['display_name'],
                'type':'directory',
                'contents':prelim_inp(file['contents'])
            })
    return output

def list_input(filters = [], sorting = [], page = 1, count = -1):
    """Fetches a list of input files from the input directory"""
    output = prelim_inp(current_app.config['storage']['manifest']['input'])
    return filterdata(output, filters, sorting, page, count)

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
    data = current_app.config['storage']['loader']()
    if id==-1:
        if str(fileID) not in data['visualize']:
            return ({
                'code':400,
                'message':'The requested file (%d) does not exist'%fileID,
                'fields':'fileID'
            },400)
        raw_reader = open(
            os.path.join(
                os.path.abspath(current_app.config['files']['data-dir']),
                'visualize',
                data['visualize'][str(fileID)]['display_name']
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
    if str(fileID) not in process[0]['files']:
        return (
            {
                "code":400,
                "message":"The requested fileID (%d) does not exist for this process (%d)" %(fileID, id),
                "fields":"fileID"
            }, 400
        )
    raw_reader = open(os.path.join(
        process[0]['output'],
        process[0]['files'][str(fileID)]['display_name']
    ))
    reader = csv.DictReader(raw_reader, delimiter='\t')
    output = {column_filter(field):field for field in reader.fieldnames}
    raw_reader.close()
    return output

def visualizefile(file):
    data = current_app.config['storage']['loader']()
    return {
            'fileID':file['fileID'],
            'description':file['description'],
            'display_name':file['display_name'],
            'is_visualizable': file['is_visualizable'],
            'visualization_type': file['visualization_type'],
            'type':'file',
            'url':'/api/v1/processes/-1/results/%s'%(file['fileID']),
            'size':(
                os.path.getsize(
                    data['visualize'][file['fileID']]['fullname']
                )
            ),
            'rows':int(subprocess.check_output([
                'wc', '-l', data['visualize'][file['fileID']]['fullname']
            ]).decode().split()[0])-1,
    }

def prelim_db(type, data):
    output = []
    for file in data:
        if file['type'] == 'file':
            if "default" in type:
                if not (re.search("/tmp/", file['display_name'])):
                    output.append(visualizefile(file))
            elif "all" in type:
                output.append(visualizefile(file))
            else:
                for filter in type:
                    if (re.search('%s.tsv'%filter, file['display_name'])):
                        output.append(visualizefile(file))
        elif file['type'] == 'directory':
            output.append({
                'display_name':file['display_name'],
                'type':'directory',
                'contents':prelim_db(type, file['contents'])
            })
    return output

def list_visualize(type = 'all', filters = [], sorting = [], page = 1, count = -1):
    data = current_app.config['storage']['loader']()
    db_dir = os.path.join(current_app.config['files']['data-dir'], 'visualize')
    tree = current_app.config['storage']['manifest']['visualize']
    output = prelim_db(type, tree)
    return filterdata(output, filters, sorting, page, count)
