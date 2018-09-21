#common utils for all controllers
import os
from glob import iglob
import time
import json
import csv
import re
import sys
import subprocess
import watchdog.events
from shlex import quote
import postgresql as psql
from postgresql.exceptions import Exception as psqlException
from .watchdir import Observe
import atexit
import site
import webbrowser
import threading
from postgresql.exceptions import UndefinedTableError
from math import ceil
import operator

class dataObj(dict):
    def __init__(self, datafiles, sync):
        super().__init__()
        super().__setitem__(
            '_datafiles',
            {datafile:[] for datafile in datafiles}
        )
        self.sync = sync

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
        self.sync.acquire()
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
        self.sync.release()

_file_info = {
    'json': {
        'description': "Metadata regarding a specific run of pVAC-Seq",
        'visualizable': False,
    },
    'yml': {
        'description': "Manifest of auxiliary files supplied to pVAC-Seq",
        'visualizable': False,
    },
    'yaml': {
        'description': "Manifest of auxiliary files supplied to pVAC-Seq",
        'visualizable': False,
    },
    'log': {
        'description': "Transcript of messages produced by pVAC-Seq",
        'visualizable': False,
    },
    'chop.tsv': {
        'description': "Processed and filtered data, with peptide cleavage data added",
        'visualizable': True,
        'visualization_type': 'full',
    },
    'all_epitopes.tsv': {
        'description': "Processed data from IEDB, but with no filtering or extra data",
        'visualizable': True,
        'visualization_type': 'full',
    },
    'filtered.tsv': {
        'description': "Processed data with all filters applied",
        'visualizable': True,
        'visualization_type': 'full',
    },
    'stab.tsv': {
        'description': "Processed and filtered data, with peptide stability data added",
        'visualizable': True,
        'visualization_type': 'full',
    },
    'filtered.condensed.ranked.tsv': {
        'description': "A condensed report of the processed and filtered data, with ranking score added",
        'visualizable': True,
        'visualization_type': 'condensed',
    },
    'tsv': {
        'description': "Raw input data parsed out of the input vcf",
        'visualizable': False,
    },
    'vcf': {
        'description': "Unprocessed input VCF",
        'visualizable': False,
    },
    'vcf.gz': {
        'description': "Unprocessed input VCF",
        'visualizable': False,
    },
}

def descriptions(ext):
    if ext in _file_info:
        return _file_info[ext]['description']
    elif re.search(r'(split|tsv)_\d+-\d+$', ext):
        return "A temporary file to cache a subset of the data when working with IEDB"
    elif re.search(r'key$', ext):
        return "Data used by pVAC-Seq to parse results from IEDB"
    return "Unknown File"

def is_visualizable(ext):
    if ext in _file_info:
        return _file_info[ext]['visualizable']
    else:
        return False

def visualization_type(ext):
    if ext in _file_info and 'visualization_type' in _file_info[ext]:
        return _file_info[ext]['visualization_type']
    else:
        return None

def check_size(file, max_tries = 5):
    tries = 1
    while tries <= max_tries:
        try:
            with open(file) as f:
                read = csv.DictReader(f, delimiter='\t')
                try:
                    next(read)
                except StopIteration:
                    if tries >= max_tries:
                        return False
        except IOError as e:
            if tries >= max_tries:
                print('Could not check size of file',file)
        time.sleep(1)
        tries += 1
    return True

def column_filter(column):
    """standardize column names"""
    return column.replace(' ', '_').replace('-', '_').lower().strip()

def loaddata(datafiles, sync):
    sync.acquire()
    data = dataObj({datafiles[datafile] for datafile in datafiles if not datafile.endswith('-dir')}, sync)
    for datafile in data['_datafiles']:
        if os.path.isfile(datafile):
            try:
                current = json.load(open(datafile))
            except BaseException as e:
                #got to make sure that lock is released so we don't stall the app
                sync.release()
                raise e
            for (key, value) in current.items():
                data.addKey(key, value, datafile)
    sync.release()
    return data

def initialize(current_app, args):
    """Setup anything that needs to be configured before the app start"""
    #This section is run once, when the API spins up
    print("Initializing app configuration")
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
    synchronizer = threading.RLock()
    data = loaddata(current_app.config['files'], synchronizer)
    if 'processid' not in data:
        data.addKey('processid', 0, current_app.config['files']['processes'])
    if 'visualize' not in data:
        data.addKey('visualize', {}, current_app.config['files']['visualize'])
    if 'input' not in data:
        data.addKey('input', {}, current_app.config['files']['input'])
    #Check the last reboot (because pid's won't remain valid after a reboot)
    current_app.config['storage']['data'] = data
    import weakref
    current_app.config['storage']['loader'] = weakref.ref(current_app.config['storage']['data'])
    loader = current_app.config['storage']['loader']
    reboot = subprocess.check_output(['last', 'reboot']).decode().split("\n")[0]
    current_app.config['reboot'] = reboot
    if 'reboot' in data and data['reboot'] != reboot:
        print("A reboot has occurred since the server was first started")
        print(
            "pid's of old pVAC-Seq runs with id's",
            data['processid'],
            "and lower may be inaccurate"
        )
    current_app.config['storage']['children']={}
    current_app.config['storage']['manifest']={}

    visapp_path = os.path.relpath(
        os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'visualizations.py'
        )
    )
    #Check if the bokeh port is already in use.  Attempt to reconnect?
    current_app.config['storage']['bokeh']=subprocess.Popen(
        'bokeh serve %s --allow-websocket-origin=localhost:8080'%(
            quote(visapp_path)
        ),
        shell=True,
        stdout=subprocess.DEVNULL
    )
    print(
        "Visualization server started on PID",
        current_app.config['storage']['bokeh'].pid
    )

    @atexit.register
    def cleanup_bokeh():
        print("Cleaning up visualization server")
        import signal
        current_app.config['storage']['bokeh'].send_signal(signal.SIGINT)
        try:
            current_app.config['storage']['bokeh'].wait(1)
        except subprocess.TimeoutExpired:
            current_app.config['storage']['bokeh'].terminate()

    #Establish a connection to the local postgres database
    try:
        tmp = psql.open("localhost/postgres")
    except psqlException as e:
        raise SystemExit("Unable to connect to your Postgres server.\
                         The pVAC-Seq API requires a running local Postgres server") from e
    if not len(tmp.prepare("SELECT 1 FROM pg_database WHERE datname = $1")('pvacseq')):
        tmp.execute("CREATE DATABASE pvacseq")
    tmp.close()
    db = psql.open("localhost/pvacseq")
    db.synchronizer = threading.RLock()
    current_app.config['storage']['db'] = db

    @atexit.register
    def cleanup_database():
        print("Cleaning up database connections")
        if 'db-clean' in current_app.config:
            with db.synchronizer:
                for table in current_app.config['db-clean']:
                    try:
                        current_app.config['storage']['db'].execute("DROP TABLE %s"%table)
                    except UndefinedTableError:
                        pass
        current_app.config['storage']['db'].close()

    #setup directory structure:
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'input'),
        exist_ok=True
    )
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'.processes'),
        exist_ok=True
    )
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'archive'),
        exist_ok=True
    )
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'visualize'),
        exist_ok=True
    )
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'export'),
        exist_ok=True
    )
    os.makedirs(
        os.path.join(current_app.config['files']['data-dir'],'.tmp'),
        exist_ok=True
    )

    def make_config():
        import yaml
        base = os.path.join(current_app.config['files']['data-dir'],'visualize')
        runs = [d for d in os.listdir(base) if os.path.isdir(os.path.join(base, d))]
        for run in runs:
            config_path = os.path.join(base, run, 'config.json')
            MHCI = os.path.join(base, run, 'MHC_Class_I', 'log', 'inputs.yml')
            MHCII = os.path.join(base, run, 'MHC_Class_II', 'log', 'inputs.yml')
            if os.path.exists(MHCI):
                with open(MHCI, 'r') as MHCI_input:
                    MHC_dict = yaml.load(MHCI_input)
                    if MHC_dict:
                        if os.path.exists(MHCII):
                            with open(MHCII, 'r') as MHCII_input:
                                temp_dict = yaml.load(MHCII_input)
                                if temp_dict:
                                    MHC_dict.update({k:v for k,v in temp_dict.items() if k not in MHC_dict})
                                    MHC_dict['alleles'].extend(temp_dict['alleles'])
                                    MHC_dict['prediction_algorithms'].extend(temp_dict['prediction_algorithms'])
                        del MHC_dict['tmp_dir']
                        MHC_dict['output'] = MHC_dict['output_dir']
                        del MHC_dict['output_dir']
                        if 'MHC_Class' in os.path.basename(MHC_dict['output']):
                            MHC_dict['output'] = MHC_dict['output'][:MHC_dict['output'].rfind('/')]
                        if os.path.exists(config_path):
                            old_dict = json.load(open(config_path))
                            if old_dict and MHC_dict != old_dict:
                                with open(config_path, 'w') as config_file:
                                    json.dump(MHC_dict, config_file, indent='\t')
                        else:
                            with open(config_path, 'w') as config_file:
                                json.dump(MHC_dict, config_file, indent='\t')
            elif os.path.exists(MHCII):
                with open(MHCII, 'r') as MHCII_input:
                    MHC_dict = yaml.load(MHCII_input)
                    if MHC_dict:
                        del MHC_dict['tmp_dir']
                        MHC_dict['output'] = MHC_dict['output_dir']
                        del MHC_dict['output_dir']
                        if 'MHC_Class' in os.path.basename(MHC_dict['output']):
                            MHC_dict['output'] = MHC_dict['output'][:MHC_dict['output'].rfind('/')]
                        if os.path.exists(config_path):
                            old_dict = json.load(open(config_path))
                            if old_dict and MHC_dict != old_dict:
                                with open(config_path, 'w') as config_file:
                                    json.dump(MHC_dict, config_file, indent='\t')
                        else:
                            with open(config_path, 'w') as config_file:
                                json.dump(MHC_dict, config_file, indent='\t')

    #checks if any previous runs results are already provided and creates subsequent config files if so
    if os.listdir(os.path.join(current_app.config['files']['data-dir'],'visualize')): make_config()

    #Setup the watchers to observe the files
    current_app.config['storage']['watchers'] = []

    inputdir = os.path.join(current_app.config['files']['data-dir'],'input')
    manifest_data = current_app.config['storage']['manifest']
    input_watcher = Observe(inputdir)
    input_watcher.subscribe(lambda x:print("Input Event:", x))

    manifest_data['input'] = []
    hier_inp = manifest_data['input']

    current = {
        os.path.join(path, filename)
        for (path, _, files) in os.walk(inputdir)
        for filename in files
    }
    for (key, filename) in data['input'].items():
        if type(data['input'][key])==str:
            ext = '.'.join(os.path.basename(filename).split('.')[1:])
            print("Updating input entry",key,"to new format")
            fullname = os.path.join(inputdir, filename)
            viz = is_visualizable(ext)
            size = check_size(fullname, 1) if viz else None
            data['input'][key] = {
                'fullname':fullname,
                'display_name':os.path.relpath(
                    filename,
                    inputdir
                ),
                'description':descriptions(ext),
                'is_visualizable': viz and size,
                'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
            }
    recorded = {item['fullname'] for item in data['input'].values()}
    targets = {k for k in data['input'] if data['input'][k]['fullname'] in recorded-current}
    for file_id in targets:
        del data['input'][file_id]
    file_id = 0
    for filename in current-recorded:
        while str(file_id) in data['input']:
            file_id += 1
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        print("Assigning file:", file_id,"-->",filename)
        fullname = os.path.abspath(os.path.join(inputdir, filename))
        viz = is_visualizable(ext)
        size = check_size(fullname, 1) if viz else None
        data['input'][str(file_id)] = {
            'fullname':fullname,
            'display_name':os.path.relpath(
                filename,
                inputdir
            ),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        }
    for filename in current:
        file_path = os.path.abspath(os.path.join(inputdir, filename))
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        file_id = str([k for k,v in data['input'].items() if v['fullname'] == file_path][0])
        viz = is_visualizable(ext)
        size = check_size(data['input'][file_id]['fullname'], 1) if viz else None
        nav_to_dir(file_path, inputdir, hier_inp).append({
            'display_name':filename[filename.rfind('/')+1:],
            'type':'file',
            'fileID':file_id,
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })

    def _create(event):
        data = loader()
        filename = os.path.relpath(
            event.src_path,
            inputdir
        )
        file_id = 0
        while str(file_id) in data['input']:
            file_id += 1
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        print("Creating file:", file_id, "-->",filename)
        viz = is_visualizable(ext)
        size = check_size(event.src_path) if viz else None
        data['input'][str(file_id)] = {
            'fullname':os.path.abspath(os.path.join(
                inputdir,
                filename
            )),
            'display_name':filename,
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        }
        nav_to_dir(event.src_path, inputdir, hier_inp).append({
            'display_name':filename[filename.rfind('/')+1:],
            'type':'file',
            'fileID':str(file_id),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })
        data.save()
    input_watcher.subscribe(
        _create,
        watchdog.events.FileCreatedEvent
    )

    def _delete(event):
        data = loader()
        filename = os.path.relpath(
            event.src_path,
            inputdir
        )
        current = nav_to_dir(event.src_path, inputdir, hier_inp)
        for entity in current:
            if entity['display_name'] == filename[filename.rfind('/')+1:]:
                current.remove(entity)
        clean_tree(hier_inp)
        for key in list(data['input']):
            if data['input'][key]['display_name'] == filename:
                del data['input'][key]
                print("Deleting file:",key,'-->', filename)
                data.save()
                return
    input_watcher.subscribe(
        _delete,
        watchdog.events.FileDeletedEvent
    )

    def _move(event):
        data = loader()
        filesrc = os.path.relpath(
            event.src_path,
            inputdir
        )
        filedest = os.path.relpath(
            event.dest_path,
            inputdir
        ) 
        ext = '.'.join(os.path.basename(filedest).split('.')[0b1:])
        viz = is_visualizable(ext)
        size = check_size(event.dest_path) if viz else None
        #This accounts for how Watchdog records duplicate symlinks (i.e. symlinks of the same file) 
        #as File Moved Events from the previously added duplicate symlink, resulting in said symlinks not being 
        #properly recorded and causing situations where the source file of such events may not also be recorded.
        current = {
            filename
            for (_, _, files) in os.walk(inputdir)
            for filename in files
        }
        if filesrc in [data['input'][k]['display_name'] for k in data['input']] and os.path.basename(filesrc) not in current:
            file_id = [k for k in data['input'] if data['input'][k]['display_name'] == filesrc][0]
            current_src = nav_to_dir(event.src_path, inputdir, hier_inp)
            current_src.remove([
                entity for entity in current_src if entity['type'] == 'file' and entity['fileID'] == str(file_id)
            ][0])
        else:
            file_id = 0
            while str(file_id) in data['input']:
                file_id += 1
                
        nav_to_dir(event.dest_path, inputdir, hier_inp).append({
            'display_name':filedest[filedest.rfind('/')+1:],
            'type':'file',
            'fileID':str(file_id),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })
        clean_tree(hier_inp)

        data['input'][str(file_id)] = {
            'fullname':os.path.abspath(os.path.join(
                inputdir,
                filedest
            )),
            'display_name':filedest,
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        }
        print("Moving file:", key,'(',filesrc,'-->',filedest,')')
        data.save()

    input_watcher.subscribe(
        _move,
        watchdog.events.FileMovedEvent
    )
    current_app.config['storage']['watchers'].append(input_watcher)

    vsz = os.path.join(current_app.config['files']['data-dir'],'visualize')
    visualize_watcher = Observe(vsz)
    visualize_watcher.subscribe(lambda x:print("visualize Event:", x))

    manifest_data['visualize'] = []
    hier_vz = manifest_data['visualize']
    #Now we set up event handlers for the visualize
    #This ensures that file ids are held consistent
    current = {
        os.path.join(path, filename)
        for (path, _, files) in os.walk(vsz)
        for filename in files
    }
    for (key, filename) in data['visualize'].items():
        if type(data['visualize'][key])==str:
            ext = '.'.join(os.path.basename(filename).split('.')[1:])
            print("Updating visualize entry",key,"to new format")
            fullname = os.path.join(vsz, filename)
            viz = is_visualizable(ext)
            size = check_size(fullname, 1) if viz else None
            data['visualize'][key] = {
                'fullname':fullname,
                'display_name':os.path.relpath(
                    filename,
                    vsz
                ),
                'description':descriptions(ext),
                'is_visualizable': viz and size,
                'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
            }
    recorded = {item['fullname'] for item in data['visualize'].values()}
    targets = {k for k in data['visualize'] if data['visualize'][k]['fullname'] in recorded-current}
    for file_id in targets:
        del data['visualize'][file_id]
    file_id = 0
    for filename in current-recorded:
        while str(file_id) in data['visualize']:
            file_id += 1
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        print("Assigning file:", file_id,"-->",filename)
        fullname = os.path.abspath(os.path.join(vsz, filename))
        viz = is_visualizable(ext)
        size = check_size(fullname, 1) if viz else None
        data['visualize'][str(file_id)] = {
            'fullname':fullname,
            'display_name':os.path.relpath(
                filename,
                vsz
            ),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        }
    for filename in current:
        file_path = os.path.abspath(os.path.join(vsz, filename))
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        file_id = str([k for k,v in data['visualize'].items() if v['fullname'] == file_path][0])
        viz = is_visualizable(ext)
        size = check_size(data['visualize'][file_id]['fullname'], 1) if viz else None
        nav_to_dir(file_path, vsz, hier_vz).append({
            'display_name':filename[filename.rfind('/')+1:],
            'type':'file',
            'fileID':file_id,
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })

    data_path = current_app.config['files']
    def _create(event):
        data = loader()
        make_config()
        filename = os.path.relpath(
            event.src_path,
            vsz
        )
        file_id = 0
        while str(file_id) in data['visualize']:
            file_id += 1
        ext = '.'.join(os.path.basename(filename).split('.')[0b1:])
        print("Creating file:", file_id, "-->",filename)
        viz = is_visualizable(ext)
        size = check_size(event.src_path) if viz else None
        data['visualize'][str(file_id)] = {
            'fullname':os.path.abspath(os.path.join(
                vsz,
                filename
            )),
            'display_name':filename,
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        }
        nav_to_dir(event.src_path, vsz, hier_vz).append({
            'display_name':filename[filename.rfind('/')+1:],
            'type':'file',
            'fileID':str(file_id),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })
        data.save()
    visualize_watcher.subscribe(
        _create,
        watchdog.events.FileCreatedEvent
    )

    def _delete(event):
        data = loader()
        filename = os.path.relpath(
            event.src_path,
            vsz
        )
        current = nav_to_dir(event.src_path, vsz, hier_vz)
        for entity in current:
            if entity['display_name'] == filename[filename.rfind('/')+1:]:
                current.remove(entity)
        clean_tree(hier_vz)
        for key in list(data['visualize']):
            if data['visualize'][key]['display_name'] == filename:
                del data['visualize'][key]
                print("Deleting file:",key,'-->', filename)
                with db.synchronizer:
                    query = db.prepare("SELECT 1 FROM information_schema.tables WHERE table_name = $1")
                    if len(query('data_visualize_'+str(key))):
                        db.execute("DROP TABLE data_visualize_"+str(key))
                data.save()
                return
    visualize_watcher.subscribe(
        _delete,
        watchdog.events.FileDeletedEvent
    )

    def _move(event):
        data = loader()
        filesrc = os.path.relpath(
            event.src_path,
            vsz
        )
        filedest = os.path.relpath(
            event.dest_path,
            vsz
        )
        file_id = [k for k in data['visualize'] if data['visualize'][k]['display_name'] == filesrc][0]
        ext = '.'.join(os.path.basename(filedest).split('.')[0b1:])
        viz = is_visualizable(ext)
        size = check_size(event.dest_path) if viz else None
        current_src = nav_to_dir(event.src_path, vsz, hier_vz)
        nav_to_dir(event.dest_path, vsz, hier_vz).append({
            'display_name':filedest[filedest.rfind('/')+1:],
            'type':'file',
            'fileID':str(file_id),
            'description':descriptions(ext),
            'is_visualizable': viz and size,
            'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
        })
        current_src.remove([
            entity for entity in current_src if entity['type'] == 'file' and entity['fileID'] == str(file_id)
        ][0])
        clean_tree(hier_vz)
        for key in data['visualize']:
            if key == file_id:
                data['visualize'][key] = {
                    'fullname':os.path.abspath(os.path.join(
                        vsz,
                        filedest
                    )),
                    'display_name':filedest,
                    'description':descriptions(ext),
                    'is_visualizable': viz and size,
                    'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                }
                print("Moving file:", key,'(',filesrc,'-->',filedest,')')
                data.save()
                return
    visualize_watcher.subscribe(
        _move,
        watchdog.events.FileMovedEvent
    )
    current_app.config['storage']['watchers'].append(visualize_watcher)

    manifest_data['results'] = []
    hier_res = manifest_data['results']

    resultdir = os.path.join(current_app.config['files']['data-dir'], '.processes')
    results_watcher = Observe(resultdir)
    results_watcher.subscribe(lambda x:print("Results Event:", x))
    for processID in range(data['processid']+1):
        processkey = 'process-%d'%processID
        if processkey in data:
            print("Checking files for process", processID)
            if 'files' in data[processkey]:
                if type(data[processkey]['files']) == list:
                    print("Updating file manifest of process",processID,"to new format")
                    for (filename, file_id) in zip(data[processkey]['files'], range(sys.maxsize)):
                        ext = '.'.join(os.path.basename(filename).split('.')[1:])
                        viz = is_visualizable(ext)
                        size = check_size(filename, 1) if viz else None
                        data[processkey]['files']={
                            file_id:{
                                'fullname':filename,
                                'display_name':os.path.relpath(
                                    filename,
                                    data[processkey]['output']
                                ),
                                'description':descriptions(
                                    '.'.join(os.path.basename(filename).split('.')[1:])
                                ),
                                'is_visualizable': viz and size,
                                'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                            }
                        }
            else:
                data[processkey]['files'] = {}
            current = {
                os.path.join(path, filename)
                for (path, _, files) in os.walk(data[processkey]['output'])
                for filename in files
            }
            recorded = {entry['fullname']:k for k,entry in data[processkey]['files'].items()}
            for file_id in recorded.keys()-current:
                print("Deleting file",file_id,"from manifest")
                file_id = recorded[file_id]
                del data[processkey]['files'][file_id]
            for filename in current-recorded.keys():
                file_id = len(data[processkey]['files'])
                while str(file_id) in data[processkey]['files']:
                    file_id += 1
                file_id = str(file_id)
                ext = '.'.join(os.path.basename(filename).split('.')[1:])
                print("Assigning file:",file_id,"-->",filename)
                viz = is_visualizable(ext)
                size = check_size(filename, 1) if viz else None
                data[processkey]['files'][file_id] = {
                    'fullname':filename,
                    'display_name':os.path.relpath(
                        filename,
                        data[processkey]['output']
                    ),
                    'description':descriptions(ext),
                    'is_visualizable': viz and size,
                    'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                }
            for filename in current:
                file_path = os.path.abspath(os.path.join(data[processkey]['output'], filename))
                ext = '.'.join(os.path.basename(filename).split('.')[1:])
                file_id = str([k for k,v in data[processkey]['files'].items() if v['fullname'] == file_path][0])
                viz = is_visualizable(ext)
                size = check_size(data[processkey]['files'][file_id]['fullname'], 1) if viz else None
                nav_to_dir(file_path, resultdir, hier_res).append({
                    'display_name':filename[filename.rfind('/')+1:],
                    'type':'file',
                    'fileID':file_id,
                    'description':descriptions(ext),
                    'is_visualizable': viz and size,
                    'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                })

    def _create(event):
        data = loader()
        parentpaths = {
            (data['process-%d'%i]['output'], i)
            for i in range(data['processid']+1)
            if 'process-%d'%i in data
        }
        filepath = event.src_path
        for (parentpath, parentID) in parentpaths:
            if os.path.commonpath([filepath, parentpath])==parentpath:
                print("New output from process",parentID)
                processkey = 'process-%d'%parentID
                file_id = len(data[processkey]['files'])
                while str(file_id) in data[processkey]['files']:
                    file_id+=1
                file_id = str(file_id)
                display_name = os.path.relpath(
                    filepath,
                    data[processkey]['output']
                )
                ext = '.'.join(os.path.basename(filepath).split('.')[1:])
                print("Assigning id",file_id,'-->',display_name)
                viz = is_visualizable(ext)
                size = check_size(filepath) if viz else None
                data[processkey]['files'][file_id] = {
                    'fullname':filepath,
                    'display_name':display_name,
                    'description':descriptions(ext),
                    'is_visualizable': viz and size,
                    'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                }
                nav_to_dir(filepath, resultdir, hier_res).append({
                    'display_name':filepath[filepath.rfind('/')+1:],
                    'type':'file',
                    'fileID':file_id,
                    'description':descriptions(ext),
                    'is_visualizable': viz and size,
                    'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                })
                data.save()
                return
    results_watcher.subscribe(
        _create,
        watchdog.events.FileCreatedEvent
    )

    def _delete(event):
        data = loader()
        parentpaths = {
            (data['process-%d'%i]['output'], i)
            for i in range(data['processid']+1)
            if 'process-%d'%i in data
        }
        filepath = event.src_path
        current = nav_to_dir(filepath, resultdir, hier_res)
        for entity in current:
            if entity['display_name'] == filepath[filepath.rfind('/')+1:]:
                current.remove(entity)
            clean_tree(hier_res)
        for (parentpath, parentID) in parentpaths:
            if os.path.commonpath([filepath, parentpath])==parentpath:
                print("Deleted output from process",parentID)
                processkey = 'process-%d'%parentID
                for (file_id, filedata) in list(data[processkey]['files'].items()):
                    if filedata['fullname'] == filepath:
                        del data[processkey]['files'][file_id]
                        print("Deleted file:", file_id,'-->',filepath)
                        with db.synchronizer:
                            query = db.prepare("SELECT 1 FROM information_schema.tables WHERE table_name = $1")
                            if len(query('data_%d_%s'%(parentID, file_id))):
                                db.execute("DROP TABLE data_%d_%s"%(parentID, file_id))
                data.save()
                return
    results_watcher.subscribe(
        _delete,
        watchdog.events.FileDeletedEvent
    )

    def _move(event):
        data = loader()
        filesrc = event.src_path
        filedest = event.dest_path
        parentpaths = {
            (data['process-%d'%i]['output'], i)
            for i in range(data['processid']+1)
            if 'process-%d'%i in data
        }
        srckey = ''
        destkey = ''
        for (parentpath, parentID) in parentpaths:
            if os.path.commonpath([filesrc, parentpath])==parentpath:
                srckey = 'process-%d'%parentID
            elif os.path.commonpath([filedest, parentpath]) == parentpath:
                destkey = 'process-%d'%parentID

        ext = '.'.join(os.path.basename(filedest).split('.')[1:])
        viz = is_visualizable(ext)
        size = check_size(filedest) if viz else None
        if srckey == destkey:
            for (file_id, filedata) in data[srckey]['files'].items():
                if filedata['fullname'] == filesrc:
                    nav_to_dir(filedest, resultdir, hier_res).append({
                        'display_name':filedest[filedest.rfind('/')+1:],
                        'type':'file',
                        'fileID':file_id,
                        'description':descriptions(ext),
                        'is_visualizable': viz and size,
                        'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                    })
                    current_src = nav_to_dir(filesrc, resultdir, hier_res)
                    current_src.remove([
                        entity for entity in current_src if entity['type'] == 'file' and entity['fileID'] == str(file_id)
                    ][0])
                    data[srckey]['files'][file_id] = {
                        'fullname':filedest,
                        'display_name':os.path.relpath(
                            filedest,
                            data[srckey]['output']
                        ),
                        'description':descriptions(ext),
                        'is_visualizable': viz and size,
                        'visualization_type': 'File contains no data' if viz and not size else visualization_type(ext),
                    }
        else:
            _delete(event)
            evt = lambda:None
            evt.src_path = event.dest_path
            _create(evt)
        clean_tree(hier_res)
    results_watcher.subscribe(
        _move,
        watchdog.events.FileMovedEvent
    )
    current_app.config['storage']['watchers'].append(results_watcher)


    @atexit.register
    def cleanup_watchers():
        print("Cleaning up observers")
        for watcher in current_app.config['storage']['watchers']:
            watcher.stop()
            watcher.join()

    current_app.config['storage']['synchronizer'] = synchronizer
    data.save()

    print("Initialization complete.  Booting API")


### filtering, sorting, and paging functions shared by multiple files ###
queryfilters = re.compile(r'(.+)(<=?|>=?|!=|==)(.+)')

ops = {
    '<': operator.lt,
    '<=': operator.le,
    '==': operator.eq,
    '!=': operator.ne,
    '>=': operator.ge,
    '>': operator.gt
}

#return the associated list within data that represents the directory of file_dir, 
#creating dictionaries for each directory in the path that doesn't already exist.
#Note: assumes file_dir is a path to a file and navigates to it's directory.
def nav_to_dir(file_dir, home_dir, data):
    current = data
    for d in os.path.relpath(file_dir, home_dir).split('/'):
        if not d == '.' and not d == file_dir[file_dir.rfind('/')+1:]:
            if d not in [f['display_name'] for f in current if f['type'] == 'directory']:
                current.append({
                    'display_name':d,
                    'type':'directory',
                    'contents': []
                })
            for item in current:
                if item['display_name'] == d and item['type'] == 'directory':
                    current = item['contents']
                    break
    return current

#recursively remove all directory dictionaries in data that are empty from bottom up
def clean_tree(data):
    for entity in data:
        if entity['type'] == 'directory':
            if not clean_tree(entity['contents']):
                data.remove(entity)
    return data

# see if string is a number
def is_number(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def cmp(arg1, op, arg2):
    operation = ops.get(op)
    return operation(arg1,arg2)

def fullresponse(data, page, count):
    if count == -1:
        count = len(data)
    if count == 0:
        total_pages = 0
    else:
        total_pages = ceil(len(data)/count)
    return ({
        "_meta": {
            "page":page,
            "count":count,
            "total_pages":total_pages,
            "total_count":len(data)
        },
        "result": data[(count*(page-1)):((count*page)) if (count*page)<len(data) else len(data)]
    })

def value_type(data,col):
    for entity in data:
        if entity['type'] == 'file':
            return type(entity[col])()
        elif entity['type'] == 'directory' and tree_size(entity['contents']):
            return value_type(entity['contents'], col)

def sort_tree(data,col):
    col_type = value_type(data,col[1:])
    #accounts for visualization_type which has NoneType and String objects
    col_type = str() if col_type == None else col_type
    data.sort(key=lambda x: x[col[1:]] if col[1:] in x and x[col[1:]] != None else col_type, reverse=True if col.startswith('-') else False)
    for file in data:
        if file['type'] == 'directory':
            sort_tree(file['contents'], col)

#sort filterdata
def sort_data(data, sorting, page, count, columns):
    if not len(sorting) or sorting[0]=="none":
        sorting = ['+display_name']
    if not [item for item in sorting if re.search(r'[-\s\+]type$', item) != None]:
        sorting.insert(0,'-type')
    i = len(sorting)-1
    while i > -1:
        col = sorting[i] if not sorting[i].startswith(' ') else '+' + sorting[i][1:]
        if not col.startswith('-') and not col.startswith('+'):
            return ({
                "code": 400,
                "message": "Please indicate which direction you'd like to sort by by putting a + or - in front of the column name",
                "fields": "sorting"
            }, 400)
        if col[1:] not in columns:
            return ({
                "code": 400,
                "message": "Unknown column name %s" % col[1:],
                "fields": "sorting"
            }, 400)
        sort_tree(data,col)
        i-=1
    return data#fullresponse(data, page, count)

#sort filterprocess
def sort(data, sorting, page, count, columns):
    if not len(sorting) or sorting[0]=="none":
        return fullresponse(data, page, count)
    i = len(sorting)-1
    while i > -1:
        col = sorting[i]
        if not col.startswith('-') and not col.startswith('+'):
            return ({
                "code": 400,
                "message": "Please indicate which direction you'd like to sort by by putting a + or - in front of the column name",
                "fields": "sorting"
            }, 400)
        if col[1:] not in columns:
            return ({
                "code": 400,
                "message": "Unknown column name %s" % col[1:],
                "fields": "sorting"
            }, 400)
        data = sorted(data, key=operator.itemgetter(col[1:]), reverse=True if col.startswith('-') else False)
        i-=1
    return fullresponse(data, page, count)

def rem_inv(data):
    for file in [f for f in data]:
        if file['type'] == 'file' and file['display_name'].startswith('.'):
            data.remove(file)
        elif file['type'] == 'directory':
            rem_inv(file['contents'])

def tree_size(data):
    count = 0
    for file in data:
        if file['type'] == 'file':
            count += 1
        elif file['type'] == 'directory':
            count += tree_size(file['contents'])
    return count

def col_names(data):
    col = []
    for file in data:
        if file['type'] == 'file':
            col = [name for name in file]
            break;
        elif file['type'] == 'directory':
            col = col_names(file['contents'])
    return col

def filter_tree(data, filters, columns):
    for file in [f for f in data]: #remove from data without the issues of removing while iterating over it
        if file['type'] == 'file':
            keep = True
            for j in range(len(filters)):
                f = filters[j].strip()
                if not len(f):
                    continue
                result = queryfilters.match(f)
                if not result:
                    return ({
                        "code":400,
                        "message": "Encountered an invalid filter (%s)" % f,
                        "fields": "filtering"
                    }, 400)
                colname = result.group(1)
                if colname not in columns:
                    return ({
                        "code": 400,
                        "message": "Unknown column name %s" % result.group(1),
                        "fields": "filtering"
                    }, 400)
                op = result.group(2)
                val = result.group(3)
                comp = file[colname]
                if type(comp) == int:
                    val = int(val)
                elif is_number(comp):
                    val = float(val)
                if not cmp(comp, op, val):
                    keep = False
            if not keep:
                data.remove(file)
        elif file['type'] == 'directory':
            filter_tree(file['contents'],filters,columns)

#a version of filterprocess for tree structured data
def filterdata(data, filters, sorting, page, count):
    rem_inv(data) #remove invisible files (if this becomes a filter option, remove)
    if not tree_size(data):
        return data
    columns = col_names(data)
    if not len(filters) or filters[0]=="none":
        return sort_data(data, sorting, page, count, columns)
    filter_tree(data, filters, columns)
    return sort_data(data, sorting, page, count, columns)

def filterprocess(data, filters, sorting, page, count):
    if not len(data):
        return fullresponse(data, page, count)
    columns = [name for name in data[0]]
    if not len(filters) or filters[0]=="none":
        return sort(data, sorting, page, count, columns)
    filteredlist = []
    for i in range(len(data)):
        comparisons = []
        for j in range(len(filters)):
            f = filters[j].strip()
            if not len(f):
                continue
            result = queryfilters.match(f)
            if not result:
                return ({
                    "code":400,
                    "message": "Encountered an invalid filter (%s)" % f,
                    "fields": "filtering"
                }, 400)
            colname = result.group(1)
            if colname not in columns:
                return ({
                    "code": 400,
                    "message": "Unknown column name %s" % result.group(1),
                    "fields": "filtering"
                }, 400)
            op = result.group(2)
            val = result.group(3)
            comp = data[i][colname]
            if type(comp) == int:
                val = int(val)
            # see if string is actually a number for accurate number comparisons,
            # avoiding string comparisons of numbers in cmp() function
            elif is_number(comp):
                data[i]
                val = float(val)
            if not cmp(comp, op, val):
                break
            if j == len(filters)-1:
                filteredlist.append(data[i])
    return sort(filteredlist, sorting, page, count, columns)
