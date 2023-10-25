#!/usr/bin/env python3

from ast import parse
import sys
import json
import os
import argparse

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='config.json', help='config file name')
parser.add_argument('-runData', help="Run over data", action="store_true")
parser.add_argument('-runMC', help="Run over MC", action="store_true")
parser.add_argument('--arg', help='Configuration argument')
parser.add_argument("--aod-writer-json", help = "Name of the json configuration file", action = "store", type = str)

extrargs = parser.parse_args()

# Make some checks on provided arguments
if len(sys.argv) < 3:
  print("ERROR: Invalid syntax! The command line should look like this:")
  print("  ./runAnalysis.py <yourConfig.json> <runData|runMC> [task:param:value] ...")
  sys.exit()

# Load the configuration file provided as the first parameter
config = {}
with open(extrargs.cfgFileName) as configFile:
  config = json.load(configFile)

# Check whether we run over data or MC
if not (extrargs.runMC or extrargs.runData):
  print("ERROR: You have to specify either runMC or runData !")
  sys.exit()
"LMeePbPb/scripts/runAnalysis.py" 73L, 2869C                                                                                                                   1,1           Top
