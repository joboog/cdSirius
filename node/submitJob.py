#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:18:10 2024

@author: pleeferguson
"""

import PySirius
from PySirius import SiriusSDK
from PySirius.models.feature_import import FeatureImport
from PySirius.models.formula_candidate import FormulaCandidate
from PySirius.models.structure_candidate import StructureCandidate
import os
import time
import pandas as pd
from formatSpectra import makeFeatures

"""
Set variables and run functions here
"""
# Sirius start & login variables
siriusPath = "/Users/pleeferguson/miniconda3/bin/sirius"
siriusUser = "lee.ferguson@duke.edu"
siriusPW = "121#Hudson"

# CD result file, scratch folder, sirius project space, and paths
cdResult = '/Users/pleeferguson/Documents/Projects/CleanStreak 2022/CleanStreak_extracts_FISh.cdResult'
projectSpaceName = 'test1'
projectSpacePath = '/Users/pleeferguson/scratch/'

# Sirius job parameters
# Formula prediction params
doSirius = True
profile = "ORBITRAP"
formulaCandidates = 10
MS2accuracy_ppm = 5
MS1accuracy_ppm = 2
filterByIsotopes = True
enforceLipidFormula = True
performBottomUpSearch = True
deNovoBelowMz = 400
formulaConstraints = "HCNOP[4]F[40]"
detectableElements = ['B', 'S', 'Cl', 'Se', 'Br']
formulaSearchDBs = None
timeOuts = {"numberOfSecondsPerDecomposition": 0,
            "numberOfSecondsPerInstance": 0}
# CSI-FingerID params
doCSIFID = True
structureDBs = ['DSSTOX', 'PUBCHEM']
# ClassyFire params
doClassyFire = True
# msNovelist params
doMsNovelist = False
msNovelistCandidates = 128



"""
Function blocks below
"""
# Start a local Sirius instance and connect to it
def startSirius(siriusPath, siriusUser, siriusPW):
   sdk = SiriusSDK()
   api = sdk.start_sirius(sirius_path=os.path.abspath(siriusPath), 
                          port=8080,
                          headless=True)
   api = sdk.connect("http://localhost:8080")

   time.sleep(10)

   # Set login information
   loginCreds = {'username': 'lee.ferguson@duke.edu',
                 'password': '121#Hudson',
                 'refreshToken': None}
   loginCreds = PySirius.AccountCredentials().from_dict(loginCreds)
   # Check if user is logged-in
   isLoggedIn = api.account().is_logged_in()
   if not isLoggedIn:
       try:
           # Login into SIRIUS web services and activate default subscription if available.
           api_response = api.account().login(True, 
                                              loginCreds, 
                                              fail_when_logged_in=False, 
                                              include_subs=False)
       except Exception as e:
           print("Exception when calling LoginAndAccountApi->login: %s\n" % e) 
   return(api)


# Create project space
def makeProjectSpace(api, projectSpaceName, projectSpacePath):
    path = os.path.abspath(projectSpacePath)+"/"+projectSpaceName+".sirius"
    ps_info = api.projects().create_project_space(projectSpaceName, path)
    return(ps_info)
    

# Make features from a cdResult file and import into Sirius
# Creates a list of PySirius feature import objects from a 
# list of formatted feature dicts (siriusCompounds) created
# in formatSpectra.py
def importCDfeatures(cdResult, ps_info, api):
    siriusCompounds = makeFeatures(os.path.abspath(cdResult))
    feature_import =[]
    for feat in siriusCompounds:
        feature_import_from_dict = FeatureImport.from_dict(feat)
        feature_import.append(feature_import_from_dict)
    api.features().add_aligned_features(ps_info.project_id,
                                        feature_import,
                                        profile=PySirius.InstrumentProfile("ORBITRAP"), 
                                        opt_fields = ["msData"]) 


# Configure job for Sirius processing
def configureJob(doSirius, profile, formulaCandidates, MS2accuracy_ppm,
                 MS1accuracy_ppm, filterByIsotopes, enforceLipidFormula,
                 performBottomUpSearch, deNovoBelowMz, formulaConstraints,
                 detectableElements, formulaSearchDBs, timeOuts,
                 doCSIFID, structureDBs, doClassyFire, doMsNovelist,
                 msNovelistCandidates, api):
    # Load default job submission template
    jobSub = api.jobs().get_default_job_config()
    # Set parameters for job
    jobSub.spectra_search_params.enabled = False
    jobSub.formula_id_params.enabled = doSirius
    jobSub.formula_id_params.profile = profile
    jobSub.formula_id_params.mass_accuracy_ms2ppm = float(MS2accuracy_ppm)
    jobSub.formula_id_params.filter_by_isotope_pattern = filterByIsotopes
    jobSub.formula_id_params.enforce_el_gordo_formula = enforceLipidFormula
    jobSub.formula_id_params.perform_bottom_up_search = performBottomUpSearch
    jobSub.formula_id_params.perform_denovo_below_mz = deNovoBelowMz
    jobSub.formula_id_params.enforced_formula_constraints = formulaConstraints
    jobSub.formula_id_params.detectable_elements = detectableElements
    jobSub.formula_id_params.formula_search_dbs = formulaSearchDBs
    jobSub.formula_id_params.ilp_timeout = timeOuts
    jobSub.fingerprint_prediction_params.enabled = doCSIFID if doSirius else False
    jobSub.structure_db_search_params.enabled = doCSIFID if doSirius else False
    jobSub.structure_db_search_params.structure_search_dbs = structureDBs
    jobSub.canopus_params.enabled = doClassyFire if jobSub.structure_db_search_params.enabled else False
    jobSub.ms_novelist_params.enabled = doMsNovelist if jobSub.structure_db_search_params.enabled else False
    jobSub.config_map = {'MS1MassDeviation.allowedMassDeviation': 
                         str(MS1accuracy_ppm)+' ppm'}
    return(jobSub)

# Start Sirius processing of configured job
def executeSirius(api, ps_info, jobSub):
    # Submit job
    job = api.jobs().start_job(project_id=ps_info.project_id, job_submission=jobSub)
    while True:
        if api.jobs().get_job(ps_info.project_id, job.id).progress.state != 'DONE':
            time.sleep(10)
        else:
            break

# Get Sirius results for import to CD
def retrieveSiriusResults(api, ps_info, jobSub):
    # Get feature IDs
    featureIds = [fid.aligned_feature_id for fid in api.features().get_aligned_features(ps_info.project_id)]
    
    # Initialize lists to populate results, depending on processes performed
    if jobSub.formula_id_params.enabled:
        formulaResults =[]
    if jobSub.structure_db_search_params.enabled: 
        structureResults =[]
    if jobSub.canopus_params.enabled:
        cmpdClassResults = []
    if jobSub.ms_novelist_params.enabled:
        deNovoStructureResults =[]
    
    # Retrieve results for all features    
    for featId in featureIds:
        # Get CD compound ID value
        CDcid = api.features().get_aligned_feature(ps_info.project_id, featId).external_feature_id
        # Get formula predictions
        if jobSub.formula_id_params.enabled:
            formulas = api.features().get_formula_candidates(ps_info.project_id, 
                                                             featId,
                                                              opt_fields=["statistics","compoundClasses"])
            formulas_dicts = [FormulaCandidate.to_dict(formResult) for 
                              formResult in formulas]
            for cmpdClass in formulas_dicts:
                if jobSub.canopus_params.enabled and cmpdClass['compoundClasses'] is not None:
                    cmpdClass_df = pd.DataFrame.from_dict(cmpdClass['compoundClasses']['classyFireLineage'])
                    cmpdClass_df = cmpdClass_df.drop(columns=['type','index'])
                    cmpdClass_df['formulaId'] = cmpdClass['formulaId']
                    cmpdClass_df['CDcid'] = CDcid
                    cmpdClassResults.append(cmpdClass_df)
                    
            
            formula_df = pd.DataFrame.from_dict(formulas_dicts)
            if not formula_df.empty:
                formula_df = formula_df.drop('compoundClasses', axis = 1)
                formula_df['CDcid'] = CDcid
                formulaResults.append(formula_df)
                
        # Get database structure predictions if available
        if jobSub.structure_db_search_params.enabled:
            structures = api.features().get_structure_candidates(ps_info.project_id, 
                                                                 featId,
                                                                 opt_fields=["dbLinks"])
            structures_dicts = [StructureCandidate.to_dict(structResult) for 
                                structResult in structures]
            structure_df = pd.DataFrame.from_dict(structures_dicts)
            if not structure_df.empty:
                structure_df['CDcid'] = CDcid
                structureResults.append(structure_df)
                
        # Get de novo structure predictions if available
        if jobSub.ms_novelist_params.enabled:
            deNovoStructures = api.features().get_de_novo_structure_candidates(ps_info.project_id, 
                                                                 featId,
                                                                 opt_fields=["dbLinks"])
            deNovoStructures_dicts = [StructureCandidate.to_dict(deNovoResult) for 
                                deNovoResult in deNovoStructures]
            deNovoStructure_df = pd.DataFrame.from_dict(deNovoStructures_dicts)
            if not deNovoStructure_df.empty:
                deNovoStructure_df['CDcid'] = CDcid
                deNovoStructureResults.append(deNovoStructure_df)
    
    results_dict = dict()
    if jobSub.formula_id_params.enabled:
        siriusFormulas = pd.concat(formulaResults, ignore_index=True)
        results_dict['Formulas'] = siriusFormulas
    if jobSub.structure_db_search_params.enabled: 
        siriusStructures = pd.concat(structureResults, ignore_index=True)
        results_dict['Structures'] = siriusStructures
    if jobSub.canopus_params.enabled:
        siriusCmpdClasses = pd.concat(cmpdClassResults, ignore_index = True)
        results_dict['Classes'] = siriusCmpdClasses
    if jobSub.ms_novelist_params.enabled:
        siriusDeNovoStructures = pd.concat(deNovoStructureResults, ignore_index=True)
        results_dict['deNovoStructures'] = siriusDeNovoStructures    
    return(results_dict)


# Shut down API  
def shutdownSirius():
    sdk = SiriusSDK()        
    sdk.shutdown_sirius()


        

   







