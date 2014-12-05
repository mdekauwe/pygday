#!/usr/bin/env python
"""
Build the G'Day model

that's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (13.08.2011)"
__email__ = "mdekauwe@gmail.com"

import os, subprocess, re
from distutils.core import setup, Command
from distutils.command.sdist import sdist as _sdist

import subprocess

def info(x):
    gitArg = ["git"]
    gitArg.extend(x)
    git_info = subprocess.check_output(gitArg)
    return(git_info.split('\n'))

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# version'

__version__ = '%s'
"""

def update_version_py():
    if not os.path.isdir(".git"):
        print "This does not appear to be a Git repository."
        return
    try:
        ver = info(["rev-parse","HEAD"])[0]
    except EnvironmentError:
        print "unable to run git, leaving src/_version.py alone"
        return
    f = open("src/_version.py", "w")
    f.write(VERSION_PY % ver)
    f.close()
    print "set src/_version.py to '%s'" % ver

def get_version():
    try:
        f = open("src/_version.py")
    except EnvironmentError:
        return None
    for line in f.readlines():
        mo = re.match("__version__ = '([^']+)'", line)
        if mo:
            ver = mo.group(1)
            return ver
    return None

class Version(Command):
    description = "update _version.py from Git repo"
    user_options = []
    boolean_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        update_version_py()
        print "Version is now", get_version()

class sdist(_sdist):
    def run(self):
        update_version_py()
        # unless we update this, the sdist command will keep using the old
        # version
        self.distribution.metadata.version = get_version()
        return _sdist.run(self)

setup(name="pygday",
    version=get_version(),
    description="A Python implementation of the G'DAY (Generic Decomposition And Yield) model",
    long_description="GDAY model, simulates the cycling of C, N and water in plants and soil",
    author="Martin De Kauwe",
    author_email='mdekauwe@gmail.com',
    platforms = ['any'],
    package_dir = {'gday': 'src'},
    packages = ['gday'],
    cmdclass={"version": Version, "sdist": sdist },
)