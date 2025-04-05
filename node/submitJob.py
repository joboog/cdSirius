#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 16:18:10 2024

Functions for working with the PySirius API to Sirius 6.0.  These functions
are called by the main functions within the cdSirius node.  Compound Discoverer
results are parsed by functions from formatSpectra for input to Sirius.

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
Function blocks below
"""
# Start a local Sirius instance and connect to it
def startSirius(siriusPath, siriusUser, siriusPW):
   sdk = SiriusSDK()
   api = sdk.start_sirius(sirius_path=os.path.abspath(siriusPath), 
                          port=8080,
                          headless=True)
   #api = sdk.connect("http://localhost:8080")

   time.sleep(10)

   # Set login information
   loginCreds = {'username': siriusUser,
                 'password': siriusPW,
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
           print(f"Sirius user {api_response.username} logged in")
       except Exception as e:
           print("Exception when calling LoginAndAccountApi->login: %s\n" % e) 
   return(api)


# Create project space
def makeProjectSpace(api, projectSpaceName, projectSpacePath):
    path = os.path.abspath(projectSpacePath)+"/"+projectSpaceName+".sirius"
    path = os.path.normpath(path)
    projectAPI = api.projects()
    ps_info = projectAPI.create_project(projectSpaceName, path_to_project=os.path.normpath(path))
    return(ps_info)
    

# Make features from a cdResult file and import into Sirius
# Creates a list of PySirius feature import objects from a 
# list of formatted feature dicts (siriusCompounds) created
# in formatSpectra.py
def importCDfeatures(cdResult, CheckedOnly, MinPeakRating, MaxMass, Limit, ps_info, api):
    siriusCompounds = makeFeatures(cdResult, 
                                   CheckedOnly, 
                                   MinPeakRating, 
                                   MaxMass, 
                                   Limit)
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
                 doCSIFID, structureDBs, PubChemFallback, doClassyFire, 
                 doMsNovelist, msNovelistCandidates, api):
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
    if PubChemFallback:
        jobSub.structure_db_search_params.expansive_search_confidence_mode=api.models().ConfidenceMode.APPROXIMATE
    else:
        jobSub.structure_db_search_params.expansive_search_confidence_mode=api.models().ConfidenceMode.OFF
    jobSub.canopus_params.enabled = doClassyFire if jobSub.structure_db_search_params.enabled else False
    jobSub.ms_novelist_params.enabled = doMsNovelist if jobSub.structure_db_search_params.enabled else False
    jobSub.config_map = {'MS1MassDeviation.allowedMassDeviation': 
                         str(MS1accuracy_ppm)+' ppm'}
    return(jobSub)

# Start Sirius processing of configured job
def executeSirius(api, ps_info, jobSub):
    # Submit job
    job = api.jobs().start_job(project_id=ps_info.project_id, job_submission=jobSub)
    command = api.jobs().get_job(project_id=ps_info.project_id, job_id=job.id, opt_fields=["command"])
    print("Sirius job started with configuration: "+command.command)
    while True:
        if api.jobs().get_job(ps_info.project_id, job.id).progress.state != 'DONE':
            time.sleep(10)
        else:
            print("Sirius job completed successfully")
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
            # Get compound class predictions if available
            for cmpdClass in formulas_dicts:
                if jobSub.canopus_params.enabled and cmpdClass['compoundClasses'] is not None:
                    cmpdClass_df = pd.DataFrame.from_dict(cmpdClass['compoundClasses']['classyFireLineage'])
                    cmpdClass_df.drop(columns=['type','index','levelIndex','parentId','parentName'], inplace = True)
                    cmpdClass_df['SiriusFormulas ID'] = cmpdClass['formulaId']
                    cmpdClass_df['Compounds ID'] = CDcid
                    cmpdClass_df.rename(columns={'id': 'Class ID',
                                                 'level': 'Classification Level',
                                                 'name': 'Classification Name',
                                                 'description': 'Description',
                                                 'probability': 'Probability'}, inplace = True)
                    cmpdClassResults.append(cmpdClass_df)
                    
            
            formula_df = pd.DataFrame.from_dict(formulas_dicts)
            if not formula_df.empty:
                formula_df['ppmError'] = pd.DataFrame(formula_df['medianMassDeviation'].tolist())['ppm']
                formula_df.drop(columns=['medianMassDeviation',
                                         'compoundClasses',
                                         'siriusScoreNormalized',
                                         'zodiacScore',
                                         'fragmentationTree',
                                         'annotatedSpectrum',
                                         'isotopePatternAnnotation',
                                         'lipidAnnotation',
                                         'predictedFingerprint',
                                         'canopusPrediction'], inplace = True)
                formula_df['Compounds ID'] = CDcid
                formula_df.rename(columns={'formulaId': 'SiriusFormulas ID',
                                           'molecularFormula': 'Formula',
                                           'adduct': 'Adduct',
                                           'rank': 'Rank',
                                           'siriusScore': 'Sirius Score',
                                           'isotopeScore': 'Isotope Score',
                                           'treeScore': 'Tree Score',
                                           'numOfExplainedPeaks': '# Explained Peaks',
                                           'numOfExplainablePeaks': '# Explainable Peaks',
                                           'totalExplainedIntensity': 'Explained Intensity',
                                           'ppmError': 'ΔMass [ppm]'}, inplace = True)
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
                dbLinks = structure_df['dbLinks']
                dbLinkList = []
                for dbLink in dbLinks:
                    PubChemID = next((item for item in dbLink if item["name"] == "PUBCHEM"), None)
                    if PubChemID:
                        PubChemID = PubChemID['id']
                    DSSToxID = next((item for item in dbLink if item["name"] == "DSSTox"), None)
                    if DSSToxID:
                        DSSToxID = DSSToxID['id']
                    dbLinkList.append({"PubChem ID": PubChemID, "DSSTox ID": DSSToxID})
                dbLink_df = pd.DataFrame(dbLinkList)
                structure_df = pd.concat([structure_df, dbLink_df], axis = 1)
                structure_df.drop(['dbLinks',
                                   'spectralLibraryMatches',
                                   'mcesDistToTopHit'], axis = 1, inplace = True)
                structure_df['Compounds ID'] = CDcid
                structure_df.rename(columns={'formulaId': 'SiriusFormulas ID',
                                             'inchiKey': 'InChIKey',
                                             'smiles': 'SMILES',
                                             'structureName': 'Name',
                                             'xlogP': 'Log Kow',
                                             'rank': 'Rank',
                                             'csiScore': 'CSI Score',
                                             'tanimotoSimilarity': 'Tanimoto Similarity',
                                             'molecularFormula': 'Formula',
                                             'adduct': 'Adduct'}, inplace = True)
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
                deNovodbLinks = deNovoStructure_df['dbLinks']
                deNovodbLinkList = []
                for deNovodbLink in deNovodbLinks:
                    deNovoPubChemID = next((item for item in deNovodbLink if item["name"] == "PUBCHEM"), None)
                    if deNovoPubChemID:
                        deNovoPubChemID = deNovoPubChemID['id']
                    deNovoDSSToxID = next((item for item in deNovodbLink if item["name"] == "DSSTox"), None)
                    if deNovoDSSToxID:
                        deNovoDSSToxID = deNovoDSSToxID['id']
                    deNovodbLinkList.append({"PubChem ID": deNovoPubChemID, "DSSTox ID": deNovoDSSToxID})
                deNovodbLink_df = pd.DataFrame(deNovodbLinkList)
                deNovoStructure_df = pd.concat([deNovoStructure_df, deNovodbLink_df], axis = 1)
                deNovoStructure_df.drop(['dbLinks',
                                         'spectralLibraryMatches',
                                         'mcesDistToTopHit'], axis = 1, inplace = True)
                deNovoStructure_df['Compounds ID'] = CDcid
                deNovoStructure_df.rename(columns={'formulaId': 'SiriusFormulas ID',
                                             'inchiKey': 'InChIKey',
                                             'smiles': 'SMILES',
                                             'structureName': 'Name',
                                             'xlogP': 'Log Kow',
                                             'rank': 'Rank',
                                             'csiScore': 'CSI Score',
                                             'tanimotoSimilarity': 'Tanimoto Similarity',
                                             'molecularFormula': 'Formula',
                                             'adduct': 'Adduct'}, inplace = True)
                deNovoStructureResults.append(deNovoStructure_df)
    
    results_dict = dict()
    if jobSub.formula_id_params.enabled:
        siriusFormulas = pd.concat(formulaResults, ignore_index=True)
        results_dict['SiriusFormulas'] = siriusFormulas
    if jobSub.structure_db_search_params.enabled: 
        siriusStructures = pd.concat(structureResults, ignore_index=True)
        siriusStructures['SiriusStructures ID'] = siriusStructures.index+1
        results_dict['SiriusStructures'] = siriusStructures
    if jobSub.canopus_params.enabled:
        siriusCmpdClasses = pd.concat(cmpdClassResults, ignore_index = True)
        siriusCmpdClasses['SiriusClasses ID'] = siriusCmpdClasses.index+1
        results_dict['SiriusClasses'] = siriusCmpdClasses
    if jobSub.ms_novelist_params.enabled:
        siriusDeNovoStructures = pd.concat(deNovoStructureResults, ignore_index=True)
        siriusDeNovoStructures['SiriusDeNovoStructures ID'] = siriusDeNovoStructures.index+1
        results_dict['SiriusDeNovoStructures'] = siriusDeNovoStructures    
    return(results_dict)


# Shut down API  
def shutdownSirius():
    sdk = SiriusSDK()        
    sdk.shutdown_sirius()


        

   







