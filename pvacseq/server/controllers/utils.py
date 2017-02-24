#common utils for all controllers
import os
import json
import subprocess
from flask import current_app
import watchdog.events
import postgresql as psql
from postgresql.exceptions import Exception as psqlException
from .watchdir import Observe

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

def loaddata(datafiles):
    data = {
        '_datafiles':{datafiles[datafile]:[] for datafile in datafiles if not datafile.endswith('-dir')}
    }
    for datafile in data['_datafiles']:
        if os.path.isfile(datafile):
            current = json.load(open(datafile))
            for (key, value) in current.items():
                data['_datafiles'][datafile].append(key)
                data[key] = value
    return data

def savedata(data):
    """Saves the data object to the configfile"""
    for datafile in data['_datafiles']:
        writer = open(datafile, 'w')
        json.dump(
            {
                key:data[key] for key in data['_datafiles'][datafile] if key in data
            },
            writer,
            indent='\t'
        )
        writer.close()


def initialize():
    """Setup anything that needs to be configured in the app, but within\
    the context of the controllers"""
    if not current_app.config['initialized']:
        #This section is run once, when the API responds to the first request
        current_app.config['initialized'] = True

        config_dir = os.path.join(
            os.path.dirname(__file__),
            '..',
            'config'
        )
        #read the index of config files from root.json
        #these cannot be overridden
        reader = open(os.path.join(config_dir, 'root.json'))
        config = {
            "_configroot":json.load(reader),
            "storage": {}
        }
        reader.close()
        user_config_dir = os.path.expanduser("~/.pvacseq")
        if not os.path.isdir(user_config_dir):
            os.makedirs(user_config_dir)
        #For ever config file listed, check if the user has defined overrides in
        #user_config_dir/file.json
        for key in config['_configroot']:
            reader = open(os.path.expanduser(os.path.join(config_dir, config['_configroot'][key])))
            config[key] = json.load(reader)
            reader.close()
            try:
                reader = open(os.path.expanduser(os.path.join(user_config_dir, config['_configroot'][key])))
                if key == 'schema':
                    config[key].update({
                        column_filter(key):value for (key,value) in json.load(reader).items()
                    })
                else:
                    config[key].update(json.load(reader))
                reader.close()
            except FileNotFoundError:
                pass
        for key in config['files']:
            config['files'][key] = os.path.abspath(os.path.expanduser(config['files'][key]))
        current_app.config.update(config) #save to the app configuration object
        data = loaddata(current_app.config['files'])
        os.makedirs(current_app.config['files']['dropbox-dir'], exist_ok=True)
        reboot = subprocess.check_output(['last', 'reboot']).decode().split("\n")[0]
        current_app.config['reboot'] = reboot
        if 'reboot' in data and data['reboot'] != reboot:
            print("A reboot has occurred since the server was first started")
            print(
                "pid's of old pVAC-Seq runs with id's",
                data['processid'],
                "and lower may be innacurate"
            )
        current_app.config['storage'] = {'children':{}}

        watcher = Observe(current_app.config['files']['dropbox-dir'])
        current_app.config['storage']['watcher'] = watcher
        watcher.subscribe(lambda x:print("FS Event:", x))

        try:
            tmp = psql.open("localhost/postgres")
        except psqlException as e:
            raise SystemExit("Unable to connect to your Postgres server.\
                             The pVAC-Seq API requires a running local Postgres server") from e
        if not len(tmp.prepare("SELECT 1 FROM pg_database WHERE datname = 'pvacseq'")()):
            tmp.execute("CREATE DATABASE pvacseq")
        tmp.close()
        db = psql.open("localhost/pvacseq")
        current_app.config['storage']['db'] = db
        #We need to setup listeners on the dropbox so that
        #fileID's are consistent between runs, and between changes to the directory
        #ie: file 'b' with id 0.  UI requests listing, gets {0: 'b'}
        #User adds file 'a'
        #UI requests id 0, expecting 'b', but gets 'a' since the listing has changed
        #The event listeners remove that race condition
        dbr = current_app.config['files']['dropbox-dir']
        fileID = 0
        current = set(os.listdir(current_app.config['files']['dropbox-dir']))
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

    else:
        data = loaddata(current_app.config['files'])

    return data
