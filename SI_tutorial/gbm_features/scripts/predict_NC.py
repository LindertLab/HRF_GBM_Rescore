import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV,train_test_split,LeaveOneOut, KFold
import sys
import pickle
import lightgbm as lgb
import matplotlib.pyplot as plt
import shap
from sklearn.preprocessing import LabelEncoder,PolynomialFeatures
from scipy import stats
from sklearn.impute import SimpleImputer
import joblib

pdbid, chain = sys.argv[1], sys.argv[2]

res = ['CYS','MET','TRP','TYR','PHE','HIS','LEU',
       'ILE','ARG','LYS','VAL','THR','SER','PRO','GLU',
       'GLN','ASN','ASP','ALA','GLY','X']

score_term_res = ['TRP','TYR','PHE','HIS','LEU','ILE','ARG']

good_feat_names = ['RES_0_SCE', 'RES_0_PSSM_5', 'RES_n1_PSI_1', 'RES_0_PSSM_40', 'RES_0_PF', 'RES_0_PSSM_0', 'BL4', 'RES_0_PSSM_13', 'RES_n3_PSSM_26', 'RES_n3_PSSM_21', 'RES_n3_PSSM_2', 'RES_p1_PSI_2', 'RES_p3_PSSM_1', 'RES_0_PSSM_24', 'RES_p1_PSI_3', 'RES_n2_PSSM_20', 'RES_n2_PSSM_9', 'FP8', 'RES_0_MOD', 'RES_n2_PSI_3', 'BL10', 'RES_p2_PSSM_39', 'RES_p1_PSSM_35', 'RES_n3_PSSM_29', 'RES_p1_PSSM_40', 'RES_p3_PSI_2', 'RES_0_PSSM_28', 'RES_0_PSSM_15', 'RES_n1_PSSM_26', 'FP5', 'RES_p2_IDENT', 'RES_p1_PSSM_41', 'RES_n2_PSSM_40', 'RES_n1_IR', 'PP3', 'RES_n1_PSSM_29', 'RES_0_PSSM_1', 'RES_n1_PSSM_0','RES_0_PSSM_4']


model = pickle.load(open('../NC_regressor.pkl', 'rb'))
paper_feats = pd.read_pickle(f'{pdbid}_{chain}_input_features.pkl')
filtered_paper_feats = paper_feats[good_feat_names]

paper_preds = model.predict(filtered_paper_feats)
with open(f'{pdbid}_{chain}_predicted_NC.out','w') as outfile:
    outfile.write(f'#{pdbid}_{chain} predicted neighbor counts\n')
    with open(f'{pdbid}_{chain}_HRPF.csv','r') as infile:
        lines = infile.readlines()
        for line, pred in zip(lines, list(paper_preds)):
            spline = line.split()
            if spline[1] in score_term_res: 
                outfile.write(f'{spline[0]}\t{pred}\n')

