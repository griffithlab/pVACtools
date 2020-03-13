import inspect
import re
import os
import unittest

def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result

def pvac_directory():
    return os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

mock_fhs = []
def mock_ncbiwww_qblast(algorithm, reference, peptide, entrez_query):
    for stack in inspect.stack():
        if '/tests/' in stack.filename and not 'test_utils.py' in stack.filename:
            tool = re.compile('.+/test_(.+).py').match(stack.filename).group(1)
            break
    base_dir      = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    test_data_dir = os.path.join(base_dir, "tests", "test_data", tool)
    fh = open(os.path.join(test_data_dir, 'response_{}.xml'.format(peptide[0:100])), 'r')
    mock_fhs.append(fh)
    return fh

def close_mock_fhs():
    for fh in mock_fhs:
        fh.close()

def make_response(data, files, path):
    if not files:
        if 'length' in data:
            filename = 'response_%s_%s_%s.tsv' % (data['allele'], data['length'], data['method'])
        else:
            filename = 'response_%s_%s.tsv' % (data['allele'], data['method'])
        reader = open(os.path.join(
            path,
            filename
        ), mode='r')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.text = reader.read()
        reader.close()
        return response_obj
    else:
        basefile = os.path.basename(data['configfile'])
        reader = open(os.path.join(
            path,
            'net_chop.html' if basefile == 'NetChop.cf' else 'Netmhcstab.html'
        ), mode='rb')
        response_obj = lambda :None
        response_obj.status_code = 200
        response_obj.content = reader.read()
        reader.close()
        return response_obj

def generate_class_i_call(method, allele, length, input_file):
    reader = open(input_file, mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhci/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'length':        length,
        'user_tool':     'pVac-seq',
    })

def generate_class_ii_call(method, allele, input_file):
    reader = open(input_file, mode='r')
    text = reader.read()
    reader.close()
    return unittest.mock.call('http://tools-cluster-interface.iedb.org/tools_api/mhcii/', data={
        'sequence_text': ""+text,
        'method':        method,
        'allele':        allele,
        'user_tool':     'pVac-seq',
    })
