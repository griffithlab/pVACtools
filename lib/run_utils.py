import sys
import os
import csv
from lib.prediction_class import *

def change_permissions_recursive(path, dir_mode, file_mode):
    for root, dirs, files in os.walk(path, topdown=False):
        for dir in [os.path.join(root,d) for d in dirs]:
            os.chmod(dir, dir_mode)
        for file in [os.path.join(root, f) for f in files]:
            os.chmod(file, file_mode)
