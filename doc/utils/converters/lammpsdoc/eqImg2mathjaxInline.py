#! /usr/bin/env python3
# LAMMPS Documentation Utilities
#
# Scan for duplicate anchor labels in documentation files
#
# Copyright (C) 2019 E. Anne Gunn
# Based largely on doc_anchor_check.py by Richard Berger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import argparse
import os
import re
import shutil
import sys


# We only want to replace image lines where image is
# pulled from Eqs subfolder
image_pattern = re.compile(r'.*image:: Eqs/(.*)\.jpg')
tex_eq_pattern = re.compile(r'\$\$')
latex_begin_eq_pattern = re.compile(r'\\begin{equation}')
latex_end_eq_pattern = re.compile(r'\\end{equation}')
latex_begin_eqArray_pattern = re.compile(r'\\begin{eqnarray\*}')
latex_end_eqArray_pattern = re.compile(r'\\end{eqnarray\*}')

imageMarker = ">>>image was here"
image_marker_pattern = re.compile(r'>>>image was here')
align_pattern = re.compile(r'.*:align: center')

leading_space_pattern = re.compile(r'^ *(.*)')  # first group will match rest of line after leading spaces

modifiedRstFolder = "src/modifiedRst/"
safeRstFolder = "src/safeRst/"
# Since this is a proof of concept implementation,
# skip any rst files that are known to cause problems
skipFileList = ["pair_tersoff_zbl.rst"]

runReport = {
}


def checkForEquationStart(texLine):
    eqType = None
    texMatch = tex_eq_pattern.match(texLine)
    if texMatch:
        eqType = "texMatch"
    else:
        eqMatch = latex_begin_eq_pattern.match(texLine)
        if eqMatch:
            eqType = "eqMatch"
        else:
            eqArrayMatch = latex_begin_eqArray_pattern.match(texLine)
            if eqArrayMatch:
                eqType = "eqArrayMatch"
    return eqType


def checkForEquationEnd(texLine, eqType):
    endPattern = tex_eq_pattern
    if eqType == "texMatch":
        endPattern = tex_eq_pattern
    elif eqType == "eqMatch":
        endPattern = latex_end_eq_pattern
    elif eqType == "eqArrayMatch":
        endPattern = latex_end_eqArray_pattern
    else:
        print("***error: unexpected eqType %s, will look for tex delimiter" % eqType)

    endMatch = endPattern.match(texLine)
    endFound = endMatch is not None
    if endFound:
        print("found pattern end, line: %s" % texLine)
    return endFound


def startMathjax():
    mathjaxLines = []
    mathjaxLines.append(".. math::\n\n")
    return mathjaxLines


def endMathjax(mathjaxLines):
    mathjaxLines.append("\n")
    mathjaxLines.append("%s\n" % imageMarker)
    return mathjaxLines


def adjustIndentation(texLine):
    adjustedLine = texLine
    m = leading_space_pattern.match(texLine)
    if m:
        lineContent = m.group(1)
        adjustedLine = "   %s\n" % lineContent

    return adjustedLine


def processFile(filename):
    print("in processFile for filename: %s" % filename)
    imageCount = 0

    modifiedFileLines = []
    doWriteModifiedFile = False
    with open(filename, 'rt') as f:
        for line_number, line in enumerate(f):
            m = image_pattern.match(line)
            if m:
                fileroot = m.group(1)
                print("fileroot: {0}".format(fileroot))
                imageCount += 1
                texFilename = "src/Eqs/{0}.tex".format(fileroot)
                print("will try to open %s" % texFilename)
                eqType = None
                eqLines = []
                try:
                    with open(texFilename, 'rt', encoding='utf-8') as t:
                        print("%s file opened ok" % texFilename)
                        eqLines = startMathjax()
                        try:
                            for dummy, texLine in enumerate(t):
                                #print(texLine)
                                if eqType == None:
                                    eqType = checkForEquationStart(texLine)
                                    if eqType != None:
                                        print("equation type: {0}".format(eqType))
                                else:
                                    endFound = checkForEquationEnd(texLine, eqType)
                                    if endFound != True:
                                        # rst syntax rules require consistent indentation for block
                                        # the txt2rst converter appears to have used 3 spaces
                                        # which is a bit unusual but does mean that content lines up under the
                                        # directive word for common markup prefixes like:
                                        # .. parsed-literal:
                                        # .. math:
                                        # So, we double-clutch
                                        # - delete any leading spaces on the line as we find it
                                        # - add 3 spaces to the front of the line
                                        newEqLine = adjustIndentation(texLine)
                                        eqLines.append(newEqLine)
                                    else:
                                        eqType = None
                                        eqLines = endMathjax(eqLines)
                                        print("Equation lines will be:")
                                        print("-----------------------------")
                                        print(*eqLines, sep="\n")
                                        print("-----------------------------")
                        except UnicodeDecodeError:
                            print("UnicodeDecodeError reading file file %s, image markup will be left in place" % texFilename)
                            break
                except EnvironmentError:
                    error = "could not open source tex file {0}, line: {1}".format(texFilename, line)
                    print(error)
                    print("image markup will be left in place")
                    if filename not in runReport:
                        runReport[filename] = []
                    runReport[filename].append(error)
                    # put the image line we could not replace back into the output
                    eqLines.append(line)
                if len(eqLines) > 0:
                    modifiedFileLines.extend(eqLines)
                    doWriteModifiedFile = True
                    eqLines = []
            else:
                # not an equation line, so simply queue it up for output as is
                modifiedFileLines.append(line)
    if doWriteModifiedFile:
        # We're going to write out a modified file, so first copy the original rst
        # file into the original file folder.
        nameParts = filename.split("/")
        filenamePos = len(nameParts) - 1
        safeFilePath = "{0}{1}".format(safeRstFolder, nameParts[filenamePos])
        shutil.copyfile(filename, safeFilePath)

        print("modifiedFileLines has %d lines before align center cleanup" % len(modifiedFileLines))
        # First, go through the file and pull out the lines where there is
        # now an image file marker followed by an align center directive
        deleteLines = []
        for lineNumber, line in enumerate(modifiedFileLines):
            m = image_marker_pattern.match(line)
            if m:
                print("found image marker in line %d" % lineNumber)
                n = align_pattern.match(modifiedFileLines[lineNumber+1])
                if n:
                    print("found align center")
                    deleteLines.append(lineNumber)
                    deleteLines.append(lineNumber+1)
        #When deleting, always work from the back of the list to the front
        for lineNumber in reversed(deleteLines):
            print(lineNumber)
            del modifiedFileLines[lineNumber]
        print("modifiedFileLines has %d lines after align center cleanup" % len(modifiedFileLines))
        # Now we can actually write out the new contents
        try:
            modFilePath = "{0}{1}".format(modifiedRstFolder, nameParts[filenamePos])
            modRst = open(modFilePath, "w")
            for rstLine in modifiedFileLines:
                modRst.write(rstLine)
            modRst.close()
        except OSError:
            print('Error: Creating directory. ' + modifiedRstFolder)
    return imageCount


def main():
    fileCount = 0
    totalImageCount = 0

    parser = argparse.ArgumentParser(description='replace image markup in rst files with inline mathjax markup from .txt source of images')
    parser.add_argument('files',  metavar='file', nargs='+', help='one or more files to scan')
    parsed_args = parser.parse_args()
    print(parsed_args)

    if not os.path.exists(safeRstFolder):
        os.makedirs(safeRstFolder)
    if not os.path.exists(modifiedRstFolder):
        os.makedirs(modifiedRstFolder)

    # Because we may decide to add files to the skip list between runs,
    # if we have more than one file to process,
    # files from both original and modified folders
    # zombie modifications
    if len(parsed_args.files) > 1:
        for outputFile in os.listdir(modifiedRstFolder):
            filePath = os.path.join(modifiedRstFolder, outputFile)
            try:
                if os.path.isfile(filePath):
                    os.unlink(filePath)
            except Exception as e:
                print(e)
                sys.exit(1)
        for safeFile in os.listdir(safeRstFolder):
            filePath = os.path.join(safeRstFolder, safeFile)
            try:
                if os.path.isfile(filePath):
                    os.unlink(filePath)
            except Exception as e:
                print(e)
                sys.exit(1)

    for filename in parsed_args.files:
        print("filename: %s" % filename)
        doSkip = False
        for skipName in skipFileList:
            if filename.find(skipName) != -1:
                print("skipping file: %s" % filename)
                doSkip = True
                runReport[filename] = ["skipped based on skipFileList"]
                break
        if not doSkip:
            fileCount += 1
            ic = processFile(filename)
            totalImageCount += ic

    print("============================================")
    print("Processed %d rst files." % fileCount)
    print("Found %d image lines." % totalImageCount)

    for fileKey in runReport:
        print("--------------------------------------------")
        print("run report for %s:" % fileKey)
        print(*runReport[fileKey], sep="\n")

    print("============================================")

if __name__ == "__main__":
    main()
