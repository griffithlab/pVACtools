import argparse
import os
import sys

def main():
    print("Clearing the api's process cache file...")
    cache_dir = os.path.expanduser(os.path.join('~', '.pvacseq'))
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)
    cache_file = os.path.join(cache_dir, 'processes.json')
    try:
        with open(cache_file, 'w+') as cache:
            cache.write('{\n}'); #empty dictionary
        print("Successfully cleared the cache file.")
    except IOError as e:
        print("Could not clear the cache file.")

if __name__ == '__main__':
    main()
