# Installed when IPYTHONDIR points to tests/ipython (see tests/gaia_replay.py)

import os
import sys

sys.path.insert(0, os.path.dirname(os.environ['IPYTHONDIR']))

import gaia_replay

gaia_replay.install()
del gaia_replay, os, sys
