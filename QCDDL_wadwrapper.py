#!/usr/bin/env python
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
#
#
# This code is an analysis module for WAD-QC 2.0: a server for automated 
# analysis of medical images for quality control.
#
# The WAD-QC Software can be found on 
# https://bitbucket.org/MedPhysNL/wadqc/wiki/Home
# 
#
# Changelog:
#   20190426: Fix for matplotlib>3
#   20170825: fixed misinterpretation of auto_suffix
#   20170621: added auto_suffix param
#   20161220: removed testing stuff; removed class variabled
#   20160802: sync with pywad1.0
#   20160622: removed adding limits (now part of analyzer)
#   20160620: remove quantity and units
#
# ./QCDDL_wadwrapper.py -c Config/rf_philips_omni.json -d TestSet/StudyDRE -r results_dre.json
from __future__ import print_function

__version__ = '20190426'
__author__ = 'aschilham'

import os
# this will fail unless wad_qc is already installed
from wad_qc.module import pyWADinput
from wad_qc.modulelibs import wadwrapper_lib

if not 'MPLCONFIGDIR' in os.environ:
    import pkg_resources
    try:
        #only for matplotlib < 3 should we use the tmp work around, but it should be applied before importing matplotlib
        matplotlib_version = [int(v) for v in pkg_resources.get_distribution("matplotlib").version.split('.')]
        if matplotlib_version[0]<3:
            os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 
    except:
        os.environ['MPLCONFIGDIR'] = "/tmp/.matplotlib" # if this folder already exists it must be accessible by the owner of WAD_Processor 

import matplotlib
matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.

try:
    import pydicom as dicom
except ImportError:
    import dicom
import QCDDL_lib

# MODULE EXPECTS PYQTGRAPH DATA: X AND Y ARE TRANSPOSED!

def logTag():
    return "[QCDDL_wadwrapper] "

##### Series wrappers
def qc_series(data, results, action):
    """
    QCDDL_UMCU checks:
        Horizontal uniformity
        LowContrast
        DynamicRange

    Workflow:
        2. Check data format
        3. Build and populate qcstructure
        4. Run tests
        5. Build xml output
        6. Build artefact picture thumbnail
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    # do we want a suffix added to the results, based on table/wall or detectorname?
    try:
        auto_suffix = (str(params['auto_suffix']).lower() == 'true')
    except:
        auto_suffix = False
    print(logTag()+' auto_suffix set to ', auto_suffix)

    inputfile = data.series_filelist[0]  # give me a filename

    ## 2. Check data format
    dcmInfile, pixeldataIn, dicomMode = wadwrapper_lib.prepareInput(inputfile, headers_only=False, logTag=logTag())

    ## 3. Build and populate qcstructure
    remark = ""
    qclib = QCDDL_lib.DDLQC()
    cs = QCDDL_lib.DDLStruct(dcmInfile,pixeldataIn,dicomMode)
    cs.verbose = False # do not produce detailed logging

    ## 4. Run tests
    error,msg = qclib.QC(cs)

    ## 5. Build xml output
    ## Struct now contains all the results and we can write these to the WAD IQ database
    stand = qclib.HorizontalOrVertical(cs)
    if auto_suffix:
        idname = '_'+stand
    else:
        idname = ''
        
    labvals = qclib.ReportEntries(cs)
    tmpdict={}
    for key, val in labvals:
        varname = key+str(idname)
        results.addFloat(varname, val)

    ## 6. Build artefact picture thumbnail
    filename = 'test'+idname+'.jpg' # Use jpg if a thumbnail is desired

    qclib.saveAnnotatedImage(cs, filename)
    varname = 'AnnotatedImage'+idname
    results.addObject(varname, filename)

def acqdatetime_series(data, results, action):
    """
    Read acqdatetime from dicomheaders and write to IQC database

    Workflow:
        1. Read only headers
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    ## 1. read only headers
    dcmInfile = dicom.read_file(data.series_filelist[0][0], stop_before_pixels=True)

    dt = wadwrapper_lib.acqdatetime_series(dcmInfile)

    results.addDateTime('AcquisitionDateTime', dt) 

def header_series(data, results, action):
    """
    Read selected dicomfields and write to IQC database

    Workflow:
        1. Read only headers
        2. Run tests
        3. Build xml output
    """
    try:
        params = action['params']
    except KeyError:
        params = {}

    # do we want a suffix added to the results, based on table/wall or detectorname?
    try:
        auto_suffix = bool(params['auto_suffix'])
    except:
        auto_suffix = False
    print(logTag()+' auto_suffix set to ', auto_suffix)

    info = 'qc'

    ## 1. read only headers
    dcmInfile, pixeldataIn, dicomMode = wadwrapper_lib.prepareInput(data.series_filelist[0], headers_only=True, logTag=logTag())

    ## 2. Run tests
    qclib = QCDDL_lib.DDLQC()

    ## Table or Wall? from distances and sensitivity; for well defined protocols to be defined in DESCRIPTION field
    cs = QCDDL_lib.DDLStruct(dcmInfile,None,dicomMode)
    cs.verbose = False # do not produce detailed logging
    dicominfo = qclib.DICOMInfo(cs,info)
    if auto_suffix:
        idname = '_'+qclib.HorizontalOrVertical(cs)
    else:
        idname = ''

    ## 3. Build xml output
    varname = 'pluginversion'+idname
    results.addString(varname, str(qclib.qcversion))
    floatlist = []
    for di in dicominfo:
        varname = di[0]+idname
        if di[0] in floatlist:
            results.addFloat(varname, di[1]) 
        else:
            results.addString(varname, str(di[1])[:min(len(str(di[1])),100)])

    varname = 'room'+idname
    results.addString(varname, cs.guessroom.name) 
    varname = 'stand'+idname
    results.addString(varname, qclib.HorizontalOrVertical(cs))

if __name__ == "__main__":
    data, results, config = pyWADinput()

    # read runtime parameters for module
    for name,action in config['actions'].items():
        if name == 'acqdatetime':
            acqdatetime_series(data, results, action)

        elif name == 'header_series':
            header_series(data, results, action)
        
        elif name == 'qc_series':
            qc_series(data, results, action)

    results.write()
