import gzip
from logging import getLogger
from pathlib import Path
import pickle
import re

from pymatgen.core import Structure
from pymatgen.io.cif import CifParser

from ..util.struc_validation import validate_loaded_structures


logger = getLogger('cryspy')


def load_structures(rin):
    # ---------- input files
    files = _get_input_files()
    if not files:
        # An empty pool leaves all initial structures for random generation.
        return {}, []

    # ---------- log
    logger.info('# ---------- Load initial structure data')
    logger.info('Load initial structures from ./data/load_struc/')

    # ---------- readers
    readers = {
        'pkl': _read_pickle,
        'cif': _read_cif,
        'poscars': _read_poscars,
        'structure': _read_structure,
    }

    # ---------- read input files
    init_struc_data = {}
    for path in files:
        try:
            file_type = _get_file_type(path)
            items = readers[file_type](path)
        except Exception as e:
            raise ValueError(
                f'{path}: failed to load structures: {e}'
            ) from e

        # ------ source identifiers
        source_label = (
            'original ID'
            if file_type in ('pkl', 'poscars')
            else 'structure number'
        )

        # ------ assign new IDs
        for source_id, struc in items:
            cid = len(init_struc_data) + 1
            init_struc_data[cid] = struc
            logger.info(
                f'Loaded structure ID {cid}: '
                f'{path}, {source_label} {source_id}'
            )

    # ---------- number of valid structures
    if not init_struc_data:
        raise ValueError('No valid structures found in ./data/load_struc')

    # ---------- validate structures
    validate_loaded_structures(rin, init_struc_data)
    logger.info(f'Validated {len(init_struc_data)} loaded structures')

    # ---------- return
    return init_struc_data, files


def move_loaded_files(files):
    # ---------- no loaded files
    if not files:
        return

    # ---------- destination directory
    done_dir = Path('./data/load_struc/done')
    done_dir.mkdir(exist_ok=True)

    # ---------- move files
    for path in files:
        destination = done_dir / path.name
        path.replace(destination)
        logger.info(f'Moved loaded file: {path} -> {destination}')


def _get_input_files():
    # ---------- input directory
    load_dir = Path('./data/load_struc')
    if not load_dir.exists():
        # No input directory leaves all structures for random generation.
        return []
    if not load_dir.is_dir():
        raise ValueError(f'{load_dir} is not a directory')

    # ---------- input files
    files = sorted(
        path for path in load_dir.iterdir()
        if path.is_file() and not path.name.startswith('.')
    )
    # ---------- return
    return files


def _get_file_type(path):
    # ---------- filename without compression suffix
    name = path.name.lower()
    if name.endswith('.gz'):
        name = name[:-3]

    # ---------- file type
    if name.endswith('.pkl'):
        return 'pkl'
    if name.endswith('.cif'):
        return 'cif'
    if 'init_poscars' in name or 'opt_poscars' in name:
        return 'poscars'

    # ---------- other formats supported by pymatgen
    return 'structure'


def _open_file(path, mode):
    # ---------- compression
    open_func = gzip.open if path.suffix.lower() == '.gz' else open

    # ---------- binary or text
    if 'b' in mode:
        return open_func(path, mode)
    return open_func(path, mode, encoding='utf-8')


def _read_pickle(path):
    # ---------- load data
    with _open_file(path, 'rb') as f:
        data = pickle.load(f)
    if not isinstance(data, dict):
        raise ValueError(f'{path}: pickle data must be a dictionary')

    # ---------- valid structures
    items = [
        (cid, struc)
        for cid, struc in data.items()
        if isinstance(cid, int)
        and not isinstance(cid, bool)
        and isinstance(struc, Structure)
    ]

    # ---------- skipped entries
    n_skipped = len(data) - len(items)
    if n_skipped:
        logger.warning(
            f'{path}: skipped {n_skipped} entries with invalid IDs '
            'or non-Structure values'
        )

    # ---------- return
    return sorted(items, key=lambda item: item[0])


def _read_cif(path):
    # ---------- read all structures
    parser = CifParser(path)
    structures = parser.parse_structures(
        primitive=False,
        on_error='raise',
    )

    # ---------- return
    return list(enumerate(structures, start=1))


def _read_structure(path):
    # ---------- read a single structure
    struc = Structure.from_file(
        path,
        primitive=False,
        sort=False,
    )

    # ---------- return
    return [(1, struc)]


def _read_poscars(path):
    # ---------- read text
    with _open_file(path, 'rt') as f:
        text = f.read()

    # ---------- structure headers
    headers = list(re.finditer(
        r'^ID_([+-]?[0-9]+)[ \t]*$',
        text,
        flags=re.MULTILINE,
    ))
    if not headers:
        raise ValueError(f'{path}: no ID headers found')
    if text[:headers[0].start()].strip():
        raise ValueError(f'{path}: unexpected content before the first ID header')

    # ---------- block boundaries
    ends = [header.start() for header in headers[1:]]
    ends.append(len(text))

    # ---------- read structures
    items = []
    for header, end in zip(headers, ends):
        block = text[header.start():end]
        struc = Structure.from_str(
            block,
            fmt='poscar',
            primitive=False,
            sort=False,
        )
        items.append((int(header.group(1)), struc))

    # ---------- return
    return items
