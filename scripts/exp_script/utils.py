import bz2
import gzip
import hashlib
import io
import json
import lzma
import os
import re
import shutil
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from enum import Enum, auto
# from importlib import metadata
from pathlib import Path
from typing import List

import numpy as np


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    zstd = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    """
    Checks if a file is compressed (gzip, bzip2 or xz).

    Parameters
    ----------
    filepath : Path
        Path object pointing to a file.

    Returns
    -------
    Compression
        Returns a `Compression` enum with one of the following attributes:
        `bzip2`, `gzip`, `xz`, or `uncompressed`
    """
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        elif tuple(signature[:4]) == (0x28, 0xB5, 0x2F, 0xFD):
            return Compression.zstd
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath):
    filepath_compression = is_compressed(filepath)
    if filepath_compression is Compression.gzip:
        fin = gzip.open(filepath, "rt")
    elif filepath_compression is Compression.bzip2:
        fin = bz2.open(filepath, "rt")
    elif filepath_compression is Compression.xz:
        fin = lzma.open(filepath, "rt")
    else:
        fin = open(filepath, "r")
    try:
        yield fin
    finally:
        fin.close()


def read_file(filepath: Path, skip_header: bool = False) -> str:
    with open_file(filepath) as fin:
        if skip_header:
            next(fin)
        yield from fin


@contextmanager
def suppress_stdout(suppress=True):
    std_ref = sys.stdout
    if suppress:
        sys.stdout = open("/dev/null", "w")
        yield
    sys.stdout = std_ref


def natsort(iterable):
    return sorted(
        iterable,
        key=lambda s: [
            int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)
        ],
    )


def check_executables(executables: List[str]) -> List[str]:
    """
    Checks if the executables from the input list are available in the PATH.

    Parameters
    ----------
    executables : list of str
        Names of the executables

    Returns
    -------
    list of str
        Returns a list containing the executables that are not in the PATH.
    """
    return [executable for executable in executables if not shutil.which(executable)]


def get_md5(filepath, size=io.DEFAULT_BUFFER_SIZE):
    m = hashlib.md5()
    with open(filepath, "rb") as fin:
        b = fin.read(size)
        while len(b) > 0:
            m.update(b)
            b = fin.read(size)
        return str(m.hexdigest())


def get_n_available_cpus():
    try:
        n_cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        n_cpus = os.cpu_count()
    return n_cpus


def write_execution_info(
    module_name: str, input_file: Path, parameters: dict, output_file: Path
) -> None:
    input_file_hash = get_md5(input_file)
    start_time = datetime.now(timezone.utc).astimezone().isoformat()
    json_dump = json.dumps(
        {
            "module": module_name,
            "input": input_file.name,
            "input_md5": input_file_hash,
            "start_time": start_time,
            "parameters": parameters,
        },
        indent=4,
    )
    with open(output_file, "w") as fout: fout.write(f"{json_dump}\n")


def get_execution_info(input_file: Path):
    with open(input_file) as fin:
        execution_info = json.load(fin)
    input_md5 = execution_info["input_md5"]
    module_name = execution_info["module"]
    parameters = execution_info["parameters"]
    return input_md5, module_name, parameters


def compare_executions(
    input_file: Path,
    parameters: dict,
    execution_info_file: Path,
    only_md5: bool = False,
) -> bool:
    input_md5 = get_md5(input_file)
    previous_input_md5, _, previous_parameters = get_execution_info(execution_info_file)
    if only_md5:
        return input_md5 == previous_input_md5
    else:
        return (parameters == previous_parameters) and (input_md5 == previous_input_md5)


def logistic(x, temperature=1.0):
    return 1 / (1 + np.exp(-x / temperature))


def softmax(x, temperature=1.0, axis=1):
    x_max = np.max(x, axis=axis, keepdims=True)
    e_x = np.exp(x - x_max)
    return e_x / np.sum(e_x, axis=axis, keepdims=True)


def entropy(x):
    x = np.array(x)
    n = len(x)
    if not np.any(x):
        return np.log2(n)
    p = x / np.sum(x)
    p = p[p != 0]
    return -1 * np.dot(p, np.log2(p))


def specificity(x):
    x = np.array(x)
    if not np.any(x):
        return 0.0
    n = len(x)
    if n == 1:
        return 0.0
    ss = np.log2(n) - entropy(x)
    return ss / np.log2(n)


def rle_encode(array):
    """
    Run-length encoding of an input array
    """
    value_array = []
    count_array = []
    i = 0
    while i <= len(array) - 1:
        c = 1
        v = array[i]
        j = i
        while j < len(array) - 1 and array[j] == array[j + 1]:
            c += 1
            j += 1
        count_array.append(c)
        value_array.append(v)
        i = j + 1
    return count_array, value_array


def rle_decode(count_array, value_array):
    decoded_array = []
    for c, v in zip(count_array, value_array):
        decoded_array += [v] * c
    return decoded_array
