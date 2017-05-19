"""
Copyright (C) 2015 - 2016 OpenEye Scientific Software, Inc.
"""
import os
import sys
import json
from invoke import task, run


@task
def update_manifest(ctx):
    """
    Updates manifest.json with the correct version of chunkdb_modules
    """
    spec = json.load(open('manifest.json', 'r'))
    sys.path.append(os.path.dirname(__file__))
    spec['version'] = __import__("chunkdb_modules").__version__
    json.dump(spec, open('manifest.json', 'w'))


@task
def package(ctx):
    """
    Create Floe Package
    """
    update_manifest(ctx)
    sys.path.append(os.path.dirname(__file__))
    run("python setup.py sdist")
    run("python setup.py bdist_wheel")
