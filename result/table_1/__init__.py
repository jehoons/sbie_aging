# -*- coding: utf-8 -*-
#*************************************************************************
# Author: {Je-Hoon Song, <song.je-hoon@kaist.ac.kr>
#
# This file is part of {sbie_aging}.
#*************************************************************************
__all__ = []

from os.path import join,dirname
import json
import pandas as pd

def get_workdir(targetdir):
    return join(dirname(__file__), targetdir)

def get_config():
    return config

config = {
    'parameter': None,
    'table': {
        'a': join(get_workdir('.'), 'Table-1A-Aging-network.csv'),
        }
    }
    