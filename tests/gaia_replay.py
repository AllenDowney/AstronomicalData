"""Record and replay Gaia archive calls, so the notebook tests don't depend on the Gaia service.

The IPython startup file in tests/ipython installs this in each notebook kernel.
GAIA_REPLAY selects the mode:

    replay (default): return the saved result of each call; fail if there is none.
    record: make the real calls and save their results in tests/gaia_cache.

Results are keyed by the method name and arguments (whitespace in queries is
collapsed), so if a notebook changes a query, record again.
"""

import hashlib
import json
import os
from pathlib import Path

from astropy.table import Table
from astroquery.gaia.core import GaiaClass

CACHE = Path(__file__).parent / 'gaia_cache'
MODE = os.environ.get('GAIA_REPLAY', 'replay')

_real = {name: getattr(GaiaClass, name)
         for name in ['launch_job', 'launch_job_async', 'load_tables', 'load_table']}


def _key(name, *args, **kwargs):
    query = kwargs.get('query', args[0] if args else '')
    text = ' '.join(str(query).split())
    flags = sorted((k, v) for k, v in kwargs.items() if k != 'query')
    digest = hashlib.sha1(f'{name} {text} {flags}'.encode()).hexdigest()[:12]
    return f'{name}-{digest}'


class ReplayJob:
    def __init__(self, text, results):
        self.text = text
        self.results = results

    def __str__(self):
        return self.text

    def get_results(self):
        return self.results.copy()


class ReplayColumn:
    def __init__(self, name):
        self.name = name


class ReplayTable:
    def __init__(self, name, text='', columns=()):
        self.name = name
        self.text = text
        self.columns = [ReplayColumn(c) for c in columns]

    def get_qualified_name(self):
        return self.name

    def __str__(self):
        return self.text


def _missing(key):
    raise RuntimeError(f'No saved Gaia result for this call ({key}). '
                       'Run the notebooks with GAIA_REPLAY=record to save it.')


def _job_method(name):
    def method(self, *args, **kwargs):
        key = _key(name, *args, **kwargs)
        info, ecsv = CACHE / f'{key}.json', CACHE / f'{key}.ecsv'
        if MODE == 'record':
            job = _real[name](self, *args, **kwargs)
            CACHE.mkdir(exist_ok=True)
            job.get_results().write(ecsv, format='ascii.ecsv', overwrite=True)
            info.write_text(json.dumps({'text': str(job)}, indent=1))
            return job
        if not info.exists():
            _missing(key)
        text = json.loads(info.read_text())['text']
        return ReplayJob(text, Table.read(ecsv, format='ascii.ecsv'))
    return method


def _load_tables(self, *args, **kwargs):
    key = _key('load_tables', *args, **kwargs)
    info = CACHE / f'{key}.json'
    if MODE == 'record':
        tables = _real['load_tables'](self, *args, **kwargs)
        CACHE.mkdir(exist_ok=True)
        info.write_text(json.dumps([t.name for t in tables], indent=1))
        return tables
    if not info.exists():
        _missing(key)
    return [ReplayTable(name) for name in json.loads(info.read_text())]


def _load_table(self, *args, **kwargs):
    key = _key('load_table', *args, **kwargs)
    info = CACHE / f'{key}.json'
    if MODE == 'record':
        table = _real['load_table'](self, *args, **kwargs)
        CACHE.mkdir(exist_ok=True)
        info.write_text(json.dumps({'name': table.name, 'text': str(table),
                                    'columns': [c.name for c in table.columns]},
                                   indent=1))
        return table
    if not info.exists():
        _missing(key)
    return ReplayTable(**json.loads(info.read_text()))


def install():
    GaiaClass.launch_job = _job_method('launch_job')
    GaiaClass.launch_job_async = _job_method('launch_job_async')
    GaiaClass.load_tables = _load_tables
    GaiaClass.load_table = _load_table
