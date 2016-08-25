import sys
import os

if __name__ == '__main__':
    try:
        from ...lib.main import main
    except SystemError:
        sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        from lib.main import main

    main(sys.argv[1:])
