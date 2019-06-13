import unittest.mock
import shutil
import threading
import subprocess
import sys
import os

def mock_process(args, *others, **kwargs):
    global popen_mock # -infinity style points
    if args[0] == 'pvacseq' and args[1] == 'run':
        print("FAKE")
        src = os.path.join(
            os.path.dirname(__file__),
            'test_data',
            'api',
            'sample_%s'%args[3].replace(' ', ' ')
        )
        dest = args[args.index('-e')-1]
        return popen_mock.__class__(
            [
                sys.executable,
                '-c',
                "import shutil; shutil.rmtree('%s', ignore_errors=True); shutil.copytree('%s', '%s')"%(
                    dest,
                    src,
                    dest
                )
            ]
        )
    else:
        print("REAL")
        return popen_mock.__class__(
            args,
            *others,
            **kwargs
        )

popen_mock = unittest.mock.create_autospec(subprocess.Popen, side_effect=mock_process)

def mock_api():
    print("starting mock api...")
    pvac_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    sys.path.append(pvac_dir)
    import utils.pvacapi.controllers.staging
    utils.pvacapi.controllers.staging.subprocess.Popen = popen_mock
    utils.pvacapi.controllers.staging.hashfile = unittest.mock.create_autospec(
        utils.pvacapi.controllers.staging.hashfile,
        return_value = b'Probably going to cause some collisions'
    )
    import utils.pvacapi
    import utils.pvacapi.app
    utils.pvacapi.app.sys.argv.append('--ip-address')
    utils.pvacapi.app.sys.argv.append('localhost')
    utils.pvacapi.app.main()

if __name__ == '__main__':
    mock_api()
