#common utils for all controllers
import os
import json
import subprocess
from flask import current_app
import watchdog.events

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

def loaddata():
    configfile = os.path.join(os.path.expanduser('~'), '.pvacseq_ui')
    if os.path.isfile(configfile):
        data = json.load(open(configfile))
    else:
        print("No saved state found")
        data={
            'processid':-1
        }
    data['configfile'] = configfile
    return data

def savedata(data):
    """Saves the data object to the configfile"""
    writer = open(data['configfile'], 'w')
    json.dump(data, writer)
    writer.close()


def initialize():
    """Setup anything that needs to be configured in the app, but within\
    the context of the controllers"""
    data = loaddata()
    if not current_app.config['initialized']:
        current_app.config['initialized'] = True
        db = current_app.config['db']
        current_app.config['children'] = {}
        if 'filtertables' not in data:
            data['filtertables']={}
        reboot = subprocess.check_output(['last', 'reboot']).decode().split("\n")[0]
        if 'reboot' in data and data['reboot'] != reboot:
            print("A reboot has occurred since the server was first started")
            print(
                "pid's of old pVAC-Seq runs with id's",
                data['processid'],
                "and lower may be innacurate"
            )

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
        savedata(data)

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
            savedata(data)
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
                    savedata(data)
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
                    savedata(data)
                    return
        watcher.subscribe(
            _move,
            watchdog.events.FileMovedEvent
        )
    return data
