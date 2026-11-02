#!/usr/bin/env python3

"""
Script to download a compressed tar archive, verify its checksum, and unpack
its contents without the top-level folder of the archive into the given folder.
Used for the MathJax files that are needed by the HTML version of the manual.
"""

import hashlib
import os
import shutil
import sys
import tarfile
import tempfile
from argparse import ArgumentParser
from urllib.request import Request, urlopen

parser = ArgumentParser(prog='fetch-archive.py',
                        description="Download, verify, and unpack a compressed tar archive")
parser.add_argument("-u", "--url", action='append', required=True,
                    help="URL of the archive. Repeat option to add fallback URLs")
parser.add_argument("-s", "--sha256", required=True, help="SHA256 checksum of the archive")
parser.add_argument("-d", "--dest", required=True, help="Folder to unpack the archive into")
parser.add_argument("-k", "--keep", action='append', default=[],
                    help="Only unpack this file or folder. Repeat option to keep more")

args = parser.parse_args()
dest = os.path.abspath(args.dest)
parent = os.path.dirname(dest)
if not os.path.isdir(parent):
    sys.exit(f"Folder {parent} does not exist")

def download(url, fname):
    """Download archive from URL to file and return its SHA256 checksum"""
    checksum = hashlib.sha256()
    request = Request(url, headers={'User-Agent': 'LAMMPS-doc-build'})
    with urlopen(request, timeout=60) as response, open(fname, 'wb') as fh:
        while True:
            data = response.read(81920)
            if not data:
                break
            checksum.update(data)
            fh.write(data)
    return checksum.hexdigest()

def wanted(path):
    """Check if the path without the top-level folder is selected for unpacking"""
    if not path or '..' in path or '' in path:
        return False
    if not args.keep:
        return True
    name = '/'.join(path)
    return any(name == keep or name.startswith(keep + '/') for keep in args.keep)

# unpack into a temporary folder and rename it at the end,
# so that an incomplete folder is never taken for the real thing
tmpdir = tempfile.mkdtemp(prefix='fetch-archive-', dir=parent)
archive = os.path.join(tmpdir, 'archive.tar.gz')
unpacked = os.path.join(tmpdir, 'unpacked')
try:
    success = False
    for url in args.url:
        print(f"Downloading {url}")
        try:
            checksum = download(url, archive)
        except Exception as e:
            print(f"Download failed: {e}")
            continue
        if checksum == args.sha256.lower():
            success = True
            break
        print(f"Checksum mismatch: expected {args.sha256} but got {checksum}")
    if not success:
        sys.exit(f"Could not download the archive for {dest}")

    numfiles = 0
    with tarfile.open(archive, 'r:gz') as tar:
        for member in tar:
            path = member.name.split('/')[1:]
            if not member.isfile() or not wanted(path):
                continue
            target = os.path.join(unpacked, *path)
            os.makedirs(os.path.dirname(target), exist_ok=True)
            with tar.extractfile(member) as src, open(target, 'wb') as fh:
                shutil.copyfileobj(src, fh)
            numfiles += 1

    if numfiles == 0:
        sys.exit(f"The archive for {dest} does not contain the expected files")
    shutil.rmtree(dest, ignore_errors=True)
    os.rename(unpacked, dest)
finally:
    shutil.rmtree(tmpdir, ignore_errors=True)
