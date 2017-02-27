#common utils for all controllers
import os
from glob import iglob
import json
import subprocess
from flask import current_app
import watchdog.events
import postgresql as psql
from postgresql.exceptions import Exception as psqlException
from .watchdir import Observe

class dataObj(dict):
    def __init__(self, datafiles):
        super().__init__()
        super().__setitem__(
            '_datafiles',
            {datafile:[] for datafile in datafiles}
        )

    def __setitem__(self, key, value):
        if key not in self and key not in {k for parent in self['_datafiles'] for k in parent}:
            raise KeyError("Key %s has no associated file.  Use addKey() first"%key)
        super().__setitem__(key, value)

    def addKey(self, key, value, dest):
        """Adds a new root key to the app data storage object."""
        if dest not in self['_datafiles']:
            self['_datafiles'][dest] = [key]
        else:
            self['_datafiles'][dest].append(key)
        super().__setitem__(key, value)

    def save(self):
        """Saves the data object to the various data files"""
        for datafile in self['_datafiles']:
            os.makedirs(os.path.dirname(datafile), exist_ok=True)
            writer = open(datafile, 'w')
            json.dump(
                {
                    key:self[key] for key in self['_datafiles'][datafile] if key in self
                },
                writer,
                indent='\t'
            )
            writer.close()

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
    data = dataObj({datafiles[datafile] for datafile in datafiles if not datafile.endswith('-dir')})
    for datafile in data['_datafiles']:
        if os.path.isfile(datafile):
            current = json.load(open(datafile))
            for (key, value) in current.items():
                data.addKey(key, value, datafile)
    return data

def initialize():
    """Setup anything that needs to be configured in the app, but within\
    the context of the controllers"""
    if not current_app.config['initialized']:
        #This section is run once, when the API responds to the first request
        current_app.config['initialized'] = True

        #First, read all the json config files to load app configuration
        config = {'storage': {}}
        config_dir = os.path.join(
            os.path.dirname(__file__),
            '..',
            'config'
        )
        user_config_dir = os.path.expanduser("~/.pvacseq")
        if not os.path.isdir(user_config_dir):
            os.makedirs(user_config_dir)
        #For every config file predefined in the config directory,
        #first read and load the file, then
        #check the user config directory for an override
        for configfile in iglob(os.path.join(config_dir, '*.json')):
            reader = open(configfile)
            key = os.path.splitext(os.path.basename(configfile))[0]
            config[key] = json.load(reader)
            reader.close()
            try:
                reader = open(os.path.join(user_config_dir, os.path.basename(configfile)))
                if key == 'schema':
                    config[key].update({
                        column_filter(k):v for (k,v) in json.load(reader).items()
                    })
                else:
                    config[key].update(json.load(reader))
                reader.close()
            except FileNotFoundError:
                pass
        for key in config['files']:
            config['files'][key] = os.path.abspath(os.path.expanduser(config['files'][key]))
        current_app.config.update(config) #save to the app configuration object

        #Now load the data object from the files specified in the configuration
        data = loaddata(current_app.config['files'])
        if 'processid' not in data:
            data.addKey('processid', 0, current_app.config['files']['processes'])
        if 'dropbox' not in data:
            data.addKey('dropbox', {}, current_app.config['files']['dropbox'])
        os.makedirs(current_app.config['files']['dropbox-dir'], exist_ok=True)

        #Check the last reboot (because pid's won't remain valid after a reboot)
        reboot = subprocess.check_output(['last', 'reboot']).decode().split("\n")[0]
        current_app.config['reboot'] = reboot
        if 'reboot' in data and data['reboot'] != reboot:
            print("A reboot has occurred since the server was first started")
            print(
                "pid's of old pVAC-Seq runs with id's",
                data['processid'],
                "and lower may be innacurate"
            )
        current_app.config['storage']['children']={}

        #Setup the watcher to observe the dropbox
        watcher = Observe(current_app.config['files']['dropbox-dir'])
        current_app.config['storage']['watcher'] = watcher
        watcher.subscribe(lambda x:print("FS Event:", x))

        #Establish a connection to the local postgres database
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

        #Now we set up event handlers for the dropbox
        #This ensures that file ids are held consistent
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
        data.save()

        data_path = current_app.config['files']
        def _create(event):
            data = loaddata(data_path)
            filename = os.path.relpath(
                event.src_path,
                dbr
            )
            fileID = 0
            while str(fileID) in data['dropbox']:
                fileID += 1
            print("Creating file:", fileID, "-->",filename)
            data['dropbox'][str(fileID)] = filename
            data.save()
        watcher.subscribe(
            _create,
            watchdog.events.FileCreatedEvent
        )

        def _delete(event):
            data = loaddata(data_path)
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
                    data.save()
                    return
        watcher.subscribe(
            _delete,
            watchdog.events.FileDeletedEvent
        )

        def _move(event):
            data = loaddata(data_path)
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
                    data.save()
                    return
        watcher.subscribe(
            _move,
            watchdog.events.FileMovedEvent
        )

    else:
        #if the app was already initialized, just load the data and return
        data = loaddata(current_app.config['files'])
    return data
