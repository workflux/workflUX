#!/usr/bin/env python3
import os
from setuptools import find_packages, setup

SETUP_DIR = os.path.dirname(__file__)
README = os.path.join(SETUP_DIR, 'README.rst')

setup(
    name='cwlab',
    version='0.1.2',    
    description='A platform-agnostic, cloud-ready framework for simplified deployment of the Common Workflow Language using a graphical web interface',
    long_description=open(README).read(),
    long_description_content_type="text/x-rst",
    url='https://github.com/CompEpigen/CWLab',
    download_url="https://github.com/CompEpigen/CWLab",
    author='Kersten Henrik Breuer',
    author_email='k.breuer@dkfz.de',
    license='Apache 2.0',
    include_package_data=True,
    packages=find_packages(exclude=("scratchs")),
    entry_points={
        "console_scripts": [
            "cwlab=cwlab.__main__:main",
        ]
    },
    install_requires=['flask',
                      'flask_wtf',
                      'flask-login',
                      'flask-sqlalchemy',
                      'pyexcel',
                      'pyexcel-io',
                      'pyexcel-ods',
                      'pyexcel-ods3',
                      'pyexcel-xls',
                      'pyexcel-xlsx',
                      'PyYAML',
                      'unidecode',
                      'cwltool',
                      'psutil'               
                      ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Operating System :: POSIX', 
        'Operating System :: POSIX :: Linux',    
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: OS Independent',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: Microsoft :: Windows :: Windows 10',
        'Operating System :: Microsoft :: Windows :: Windows 8.1', 
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: System :: Distributed Computing',
    ]
)