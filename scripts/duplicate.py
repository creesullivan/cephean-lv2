#!/usr/bin/python

# Python shell script for Duplicating Cephean Audio Library LV2 projects,
#       modified from original script by Bibha Tripathi http://code.activestate.com/recipes/435904-sedawk-python-script-to-rename-subdirectories-of-a/
#       modified from WDL-OL script by Oli Larkin 2012-2019 http://www.olilarkin.co.uk
#       current version modified by Chris Sullivan 2025
#
# USAGE:
# duplicate.py [inputprojectname] [outputprojectname]
#
# Designed to work with Python 2 and 3
#
# FEATURES:
# Copies a plugin project folder located at ..\projects\inputprojectname relative to this script
# into a new plugin project folder located at ..\projects\outputprojectname with only source code
# files, find/replacing all occurrences of inputprojectname with outputprojectname in all file
# names and within all code/text files. This is not guaranteed to be safe in general, but when used
# with the plugintemplate project, this operation is always safe. 


from __future__ import generators

import fileinput, glob, string, sys, os, re, uuid, pprint, random
from shutil import copy, copytree, ignore_patterns, rmtree
from os.path import join

scriptpath = os.path.dirname(os.path.realpath(__file__))

VERSION = "1.0"

# binary files that we don't want to do find and replace inside
FILTERED_FILE_EXTENSIONS = [".ico",".icns", ".pdf", ".png", ".zip", ".exe", ".wav", ".aif", ".dll"]

# files/folders that we don't want to duplicate
DONT_COPY = (".vs", "x64", "*.exe", "*.dll", "*.dmg", "*.pkg", "*.mpkg", "*.svn", "*.ncb", "*.suo", "*sdf", "ipch", "build-*", "*.layout", "*.depend", ".DS_Store", "xcuserdata", "*.aps")

def checkdirname(name, searchproject):
  "check if directory name matches with the given pattern"
  print("")
  if name == searchproject:
    return True
  else:
    return False

def replacestrs(filename, s, r):
  files = glob.glob(filename)

  for line in fileinput.input(files,inplace=1):
    line.find(s)
    line = line.replace(s, r)
    sys.stdout.write(line)

def dirwalk(dir, searchproject, replaceproject):
  for f in os.listdir(dir):
    fullpath = os.path.join(dir, f)

    #for each folder, recurse to find each file
    if os.path.isdir(fullpath) and not os.path.islink(fullpath):
        foldername = os.path.basename(fullpath)
        newfoldername = foldername.replace(searchproject, replaceproject)
      
        print('recursing in ' + f + ' directory: ')
        for x in dirwalk(fullpath, searchproject, replaceproject):
          yield x
          
        if foldername != newfoldername:
            print("Renaming folder " + foldername + " to " + newfoldername)
            os.rename(fullpath, os.path.join(dir, newfoldername))

    #for each file, replace the old project name with the new one
    if os.path.isfile(fullpath):
      filename = os.path.basename(fullpath)
      newfilename = filename.replace(searchproject, replaceproject)
      base, extension = os.path.splitext(filename)

      if (not(extension in FILTERED_FILE_EXTENSIONS)):

        print("Replacing project name strings in file " + filename)
        replacestrs(fullpath, searchproject, replaceproject)

      else:
        print("NOT replacing name strings in file " + filename)

      if filename != newfilename:
        print("Renaming file " + filename + " to " + newfilename)
        os.rename(fullpath, os.path.join(dir, newfilename))

      yield f, fullpath
      
    else:
      yield f, fullpath

def main():
  global VERSION
  print("\nCephean LV2 Project Duplicator v" + VERSION + " by Chris Sullivan ------------------------------\n")

  numargs = len(sys.argv) - 1

  if not (numargs == 2):
    print("Usage: duplicate.py inputprojectname outputprojectname")
    sys.exit(1)
    
  inputprojectname=sys.argv[1]
  outputprojectname=sys.argv[2]

  if ' ' in inputprojectname:
    print("error: input project name has spaces")
    sys.exit(1)

  if ' ' in outputprojectname:
    print("error: output project name has spaces")
    sys.exit(1)

  # remove a trailing slash if it exists
  if inputprojectname[-1:] == "/":
    inputprojectname = inputprojectname[0:-1]

  if outputprojectname[-1:] == "/":
    outputprojectname = outputprojectname[0:-1]
    
  scriptpath=os.getcwd()
  projectpath = os.path.join(scriptpath, "..", "projects")
  inputpath = os.path.join(projectpath, inputprojectname)
  outputpath = os.path.join(projectpath, outputprojectname)

  #check that the folders are OK
  if os.path.isdir(inputpath) == False:
    print("error: input project not found at " + inputpath)
    sys.exit(1)

  if os.path.isdir(outputpath):
    print("error: output project already exists")
    sys.exit(1)

  print("copying " + inputprojectname + " folder to " + outputpath)
  copytree(inputpath, outputpath, ignore=ignore_patterns(*DONT_COPY))

  #replace strings
  for dir in dirwalk(outputpath, inputprojectname, outputprojectname):
    pass

  print("\ndone.")

if __name__ == '__main__':
  main()
