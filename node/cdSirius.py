#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:04:19 2024

@author: pleeferguson
"""

import pyeds
import sys
import os
import json
import pandas as pd
from CdScriptingNodeHelper import ScriptingResponse
from submitJob import startSirius,makeProjectSpace,importCDfeatures,configureJob,executeSirius,retrieveSiriusResults,shutdownSirius
import string
import random

def print_error(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def main():
    print('cdSirius node')
                
    # start in development mode where nodeargs are given explicitely rather than reading it as command line argument
    if sys.argv[1] == '-devel':    
        print(f'Development mode: Current Dir is {os.getcwd()}')
        nodeargs_path = 'node_args.json'        
    else:
        nodeargs_path = sys.argv[1]        

    # parse node args from Compound Discoverer
    try:
        with open(nodeargs_path, 'r') as rf:
            nodeargs = json.load(rf) 
            compounds_path = ''
            response_path = nodeargs['ExpectedResponsePath']
            cdResult_path = nodeargs['ResultFilePath']
            parameters = nodeargs['NodeParameters']
            tables = nodeargs['Tables']
            for table in tables:
                if table['TableName'] == 'Compounds':                    
                    compounds_path = table['DataFile']                                    
                    if table['DataFormat'] != 'CSV':
                        print_error(f"Unknown Data Format {table['DataFormat']}")    
                        exit(1)
    except Exception as e: 
        print_error('Could not read Compound Discoverer node args')
        print_error(str(e))
        exit(1)

    if not compounds_path:
        print_error('Compounds file not defined in node args.')
        exit(1)
    
    # Parse Sirius start & login variables from parameters
    siriusPath = str(parameters['Sirius Program Path'])
    siriusUser = str(parameters['Sirius Username'])
    siriusPW = str(parameters['Sirius Password'])

    # Parse CD result file path, scratch folder path, generate Sirius project space name
    cdResult = cdResult_path
    projectSpaceName = ''.join(random.choice(string.ascii_lowercase) for x in range(10))
    (workdir, _ ) = os.path.split(response_path)
    projectSpacePath = workdir
    
    # Parse CD compound settings
    CheckedOnly = parameters['Checked Feature Status Handling'] == "Checked"
    MinPeakRating = float(parameters['Peak Quality Threshold'])
    MaxMass = float(parameters['Maximum MW'])
    Limit = 2000 # Arbitrary limit switch for debugging

    # Sirius job parameters
    # Formula prediction params
    doSirius = parameters['Predict Formulas'] == "True"
    profile = "ORBITRAP"
    formulaCandidates = int(parameters['Maximum Formula Candidates'])
    MS2accuracy_ppm = float(parameters['MS1 Mass Accuracy [ppm]'])
    MS1accuracy_ppm = float(parameters['MS2 Mass Accuracy [ppm]'])
    filterByIsotopes = parameters['Filter by Isotope Pattern'] == "True"
    enforceLipidFormula = parameters['Enforce Lipid Detection Filtering'] == "True"
    performBottomUpSearch = parameters['Perform Bottom-Up Formula Search'] == "True"
    deNovoBelowMz = float(parameters['De novo Formula Generation Threshold'])
    formulaConstraints = str(parameters['Formula Elemental Constraints'])
    detectableElements = ['B', 'S', 'Cl', 'Se', 'Br']
    formulaSearchDBs = None
    timeOuts = {"numberOfSecondsPerDecomposition": 0,
                "numberOfSecondsPerInstance": 0}
    # CSI-FingerID params
    doCSIFID = parameters['Predict Structures']== "True"
    structureDBs = ['DSSTOX', 'PUBCHEM']
    # ClassyFire params
    doClassyFire = parameters['Predict Compound Classes']== "True"
    # msNovelist params
    doMsNovelist = parameters['Predict de Novo Structures']== "True"
    msNovelistCandidates = int(parameters['De Novo Structure Candidates Limit'])    
            
    
    # Import Compounds table from CD export   
    try:
        compoundsImport = pd.read_csv(compounds_path, 
                                      header=0,
                                      sep = '\t')

    except Exception as e:
        print_error('Could not process data')
        print_error(e)        
        exit(1)
    
    # Start the Sirius API
    try:
        api = startSirius(siriusPath, siriusUser, siriusPW)
        
    except Exception as e:
        print_error('Could not start the Sirius engine')
        print_error(e)
        exit(1)

    # Create and connect to Sirius project space
    try:
        ps_info = makeProjectSpace(api, projectSpaceName, projectSpacePath)
        
    except Exception as e:
        print_error('Could not create Sirius project space')
        print_error(e)
        exit(1)
    
    # Import features from Compound Discoverer to Sirius project space
    try:
        importCDfeatures(cdResult, 
                         CheckedOnly, 
                         MinPeakRating, 
                         MaxMass, 
                         Limit, 
                         ps_info, 
                         api)
        
    except Exception as e:
        print_error('Could not create import Compound Discoverer features to Sirius')
        print_error(e)
        exit(1)
        
    # Configure job for Sirius processing
    try:
        jobSub = configureJob(doSirius, profile, formulaCandidates, MS2accuracy_ppm,
                         MS1accuracy_ppm, filterByIsotopes, enforceLipidFormula,
                         performBottomUpSearch, deNovoBelowMz, formulaConstraints,
                         detectableElements, formulaSearchDBs, timeOuts,
                         doCSIFID, structureDBs, doClassyFire, doMsNovelist,
                         msNovelistCandidates, api)
        
    except Exception as e:
        print_error('Could not configure Sirius job input')
        print_error(e)
        exit(1)
        
    # Initiate Sirius job processing
    try:
        executeSirius(api, ps_info, jobSub)
        
    except Exception as e:
        print_error('Sirius processing failed')
        print_error(e)
        exit(1)
        
    # Import results from Sirius job completion
    try:
        results_dict = retrieveSiriusResults(api, ps_info, jobSub)
         
    except Exception as e:
        print_error('Could not retrieve Sirius processing results')
        print_error(e)
        exit(1)
        
    # Shut down Sirius API
    shutdownSirius()
         

     

    outTable = compoundsImport.merge(scoredCompounds,
                                     how = 'left',
                                     on = 'Compounds ID')
    outTable['Formula Carbons'] = outTable['Formula Carbons'].astype('Int64')
    
    # write data file
    outfilename = 'ScoredCompounds.txt'
    (workdir, _ ) = os.path.split(response_path)
    outfile_path = os.path.join(workdir, outfilename)
    outTable.to_csv(outfile_path, 
                    sep='\t', 
                    encoding='utf-8', 
                    index=False, 
                    header=True, 
                    quoting=1,
                    na_rep ='')
                
    # Assemble response JSON file
    response = ScriptingResponse()
    response.add_table('Compounds', outfile_path)
    response.add_column('Compounds', 'Compounds ID', 'Int', 'ID')
    response.add_column('Compounds', 'Background', 'Boolean')
    response.add_column('Compounds', 'Formula Carbons', 'Int')
    response.set_column_option('Compounds', 'Formula Carbons', 'RelativePosition', '122')
    response.add_column('Compounds', 'Estimated Carbons', 'Float')
    response.set_column_option('Compounds', 'Estimated Carbons', 'RelativePosition', '123')
    response.set_column_option('Compounds', 'Estimated Carbons', 'FormatString', 'F1')
    response.add_column('Compounds', 'Mass Defect Calc', 'Float')
    response.set_column_option('Compounds', 'Mass Defect Calc', 'RelativePosition', '501')
    response.set_column_option('Compounds', 'Mass Defect Calc', 'FormatString', 'F3')
    response.add_column('Compounds', 'MD/C', 'Float')
    response.set_column_option('Compounds', 'MD/C', 'RelativePosition', '502')
    response.set_column_option('Compounds', 'MD/C', 'FormatString', 'F3')
    response.add_column('Compounds', 'm/C', 'Float')
    response.set_column_option('Compounds', 'm/C', 'RelativePosition', '503')
    response.set_column_option('Compounds', 'm/C', 'FormatString', 'F1')
    response.add_column('Compounds', 'MD/Cm', 'Float')
    response.set_column_option('Compounds', 'MD/Cm', 'RelativePosition', '504')
    response.set_column_option('Compounds', 'MD/Cm', 'FormatString', 'F3')
    response.add_column('Compounds', 'm/Cm', 'Float')
    response.set_column_option('Compounds', 'm/Cm', 'RelativePosition', '505')
    response.set_column_option('Compounds', 'm/Cm', 'FormatString', 'F1')
    response.add_column('Compounds', 'rCF2', 'Float')
    response.set_column_option('Compounds', 'rCF2', 'RelativePosition', '506')
    response.set_column_option('Compounds', 'rCF2', 'FormatString', 'F3')
    
    # Save response file to disk
    response.save(response_path)




if __name__== "__main__" :
    main()