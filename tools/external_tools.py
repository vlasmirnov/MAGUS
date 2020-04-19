'''
Created on Apr 14, 2020

@author: Vlad
'''

import subprocess
from configuration import Configs

def runMafft(fastaPath, subtablePath, workingDir, outputPath, threads = 1):
    args = [Configs.mafftPath, "--localpair", "--maxiterate", "1000", "--ep", "0.123", 
            "--quiet", "--thread", str(threads), "--anysymbol"]
    if subtablePath is not None:
        args.extend(["--merge", subtablePath])
    args.extend([fastaPath, ">", outputPath])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command

def runMcl(matrixPath, inflation, workingDir, outputPath):
    args = [Configs.mclPath, matrixPath, "--abc", "-o", outputPath]
    if inflation is not None:
        args.extend(["-I", str(inflation)])
    command = subprocess.list2cmdline(args)
    Configs.log("Running command: {}".format(command))
    subprocess.run(command, shell=True, cwd = workingDir)
    return command