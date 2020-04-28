#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 15:16:18 2020

@author: Warren

Script to find fragment yeilding best SuCOs score with
conformer of docked compounds
"""

# Copyright <2019> <University of Oxford>
# This code is licensed under MIT license (see LICENSE.txt for details)

import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit import RDConfig
import csv
import multiprocessing as mp

#################################################
#### Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
#    keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic')

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
#################################################

def get_FeatureMapScore(small_m, large_m, score_mode=FeatMaps.FeatMapScoreMode.All):
    featLists = []
    for m in [small_m, large_m]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fms[0].scoreMode = score_mode
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
    return fm_score

def get_frag_match(inspiration_frags,frag_scores):
    max_values = max(frag_scores, key=lambda x:x[1])
    
    max_score_frag = max_values[0]
    max_score = max_values[1]
    
    if max_score_frag in inspiration_frags:
        return 1, max_score_frag, max_score 
    else:
        return 0, max_score_frag, max_score
    
def get_avg_sucos(sucos_scores):
    return np.mean(sucos_scores)
    

def get_sucos(frag_sdf_folder, docked_sdf_file):
    
    path  = frag_sdf_folder + '/'
    
    frag_mol_list =  [Chem.MolFromMolFile((path + sdf_file), sanitize=True) for sdf_file in os.listdir(frag_sdf_folder)]
        
    docked_mol_list = Chem.SDMolSupplier(docked_sdf_file, sanitize=True)
     
    docked_mol_list = [mol for mol in docked_mol_list if mol is not None]
    
    all_frags_scores = []
    
    for docked_mol in docked_mol_list:
        docked_name = docked_mol.GetProp('_Name')
        print('Getting values for {}'.format(docked_name))
        
        sucos_scores = []
        frags_complete = []
        for frag_mol in frag_mol_list:
            ##############################################
            ####### Feature map
            ##############################################
            fm_score = get_FeatureMapScore(frag_mol, docked_mol)
            fm_score = np.clip(fm_score, 0, 1)
            ##############################################
               
            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(frag_mol, docked_mol,
                    allowReordering=False)
            protrude_dist = np.clip(protrude_dist, 0, 1)
            SuCOS_score = 0.5*fm_score + 0.5*(1 - protrude_dist)
            sucos_scores.append(SuCOS_score)
            frags_complete.append(frag_mol.GetProp('_Name'))
        
        frag_scores = list(zip(frags_complete,sucos_scores))
        insp_frags = docked_mol.GetProp('fragments')
        found, frag, score = get_frag_match(insp_frags, frag_scores)
                        
        # get frag mol using index of max frag

        frag_mol_index = frags_complete.index(frag)
        frag_mol = frag_mol_list[frag_mol_index]
        
        # Get frag and compound SMILES
        frag_SMILES = Chem.MolToSmiles(frag_mol)
        docked_SMILES = Chem.MolToSmiles(docked_mol)
               
        # Get avg sucos scores
        avg_score = get_avg_sucos(sucos_scores)
        
        # Get all scores
        all_frags_scores.append((docked_name, frag, 
                                 found,
                                 docked_SMILES, frag_SMILES,
                                 score, avg_score))
               
        with open('sucos_scores/JC_sucos_scores.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Compound_name','Fragment',
                             'Insp_frag_found',
                             'dockedSMILES','fragSMILES',
                             'SuCOS_score','Avg_SuCOS_score'])
            writer.writerows(all_frags_scores)
            

pool = mp.Pool(mp.cpu_count())

pool.starmap(get_sucos('sdf_files', 'covid_submissions_all_info-2020-04-23-docked.sdf'))

pool.close()    
    
            

