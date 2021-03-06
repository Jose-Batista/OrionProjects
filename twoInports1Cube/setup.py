#!/usr/bin/env python
import re
import ast

from pip.req import parse_requirements

from setuptools import setup, find_packages


def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]

try:
    install_reqs = get_reqs(parse_requirements("requirements.txt"))
except TypeError:
    from pip.download import PipSession
    install_reqs = get_reqs(
        parse_requirements("requirements.txt", session=PipSession())
    )


def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('hwfloe_mod/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version

setup(
    name="OpenEye-hwfloe",
    version=get_version(),
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    author="Jose Batista",
    author_email="batista@eyesopen.com",
    description=open("README.md", "r").read(),
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ]
)
