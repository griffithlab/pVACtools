import inspect
import re
import os

def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result

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
