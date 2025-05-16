import glob
import os
import numpy as np
import sys
import pandas as pd
import pickle
from sklearn.preprocessing import LabelEncoder
import json
from collections import Counter

pdbid = sys.argv[1]
chain = sys.argv[2]

IRs =  {
        "CYS": 29.2, "MET": 20.5, "TRP": 17.4, "TYR": 12.0, "PHE": 11.2,
        "HIS": 9.3, "LEU": 4.4, "ILE": 4.4, "ARG": 2.9, "LYS": 2.2,
        "VAL": 1.9, "THR": 1.6, "SER": 1.4, "PRO": 1.0, "GLU": 0.69,
        "GLN": 0.66, "ASN": 0.44, "ASP": 0.42, "ALA": 0.14, "GLY": 0.04
        }

res_dict = {'C':'CYS','M':'MET','W':'TRP','Y':'TYR','F':'PHE','H':'HIS','L':'LEU',
       'I':'ILE','R':'ARG','K':'LYS','V':'VAL','T':'THR','S':'SER','P':'PRO','E':'GLU',
       'Q':'GLN','N':'ASN','D':'ASP','A':'ALA','G':'GLY','X':'X'}

res = ['CYS', 'MET', 'TRP', 'TYR', 'PHE', 'HIS', 'LEU', 'ILE', 'ARG', 'LYS', 'VAL', 'THR', 'SER', 'PRO', 'GLU', 'GLN', 'ASN', 'ASP', 'ALA', 'GLY', 'X']

le_res = LabelEncoder()
le_res.fit(res)

le_sse = LabelEncoder()

le_sse.fit(['C','E','H','X'])

le_sce = LabelEncoder()
le_sce.fit(['E','P2','P1','B3','B2','B1','X'])


class FeatureGenerator:
    def __init__(self, config_file):
        with open(config_file, 'r') as file:
            self.features = json.load(file)

    def generate_features(self, aa):
        generated_features = []
        for feature_name, feature_values in self.features.items():
            generated_features.append(feature_values.get(aa, None))
        return generated_features

def fetch_disopred(pdbid, chain, resnum, NA = False):
    if NA:
        return [np.nan, np.nan]
    with open(f'./{pdbid}_{chain}_disopred.out','r') as infile:
        for line in infile.readlines():
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            spline = line.split()
            if spline[0] == resnum:
                if spline[2] == '*':
                    return [1.0, float(spline[-1])]
                elif spline[2] == '.':
                    return [0.0, float(spline[-1])]

def fetch_psipred(pdbid, chain, resnum, NA = False):
    if NA:
        encoded = le_sse.transform(['X'])
        return [encoded[0], np.nan, np.nan, np.nan]
    with open(f'./{pdbid}_{chain}.ss2','r') as infile:
        for line in infile.readlines():
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            spline = line.split()
            if spline[0] == resnum:
                encoded = le_sse.transform([spline[2]])
                return [encoded[0], float(spline[3]), float(spline[4]), float(spline[5])]

def fetch_PSSM_Score(pdbid, chain, resnum, NA = False):
    pssm_score = []
    if NA:
        for i in range(42):
            pssm_score.append(np.nan)
        print(f'NAN len: {len(pssm_score)}')
        return pssm_score
    with open(f'./{pdbid}_{chain}.pssm_A','r') as infile:
        for line in infile.readlines()[3:-6]:
            spline = line.split()
            if spline[0] == resnum:
                for score in spline[2:22]:
                    score = int(score)
                    neg_score = score * -1
                    normalized_score = 1 / (1 + np.exp(neg_score))
                    pssm_score.append(normalized_score)
                for score in spline[22:42]:
                    pssm_score.append(float(score))
                pssm_score.append(float(spline[-2]))
                pssm_score.append(float(spline[-1]))
    return pssm_score 

def fetch_SCE(id_ch, res, NA = False):
    if NA:
        return np.nan
    df = pd.read_excel(f'./{id_ch}_mapped.xlsx')
    homolog_SCE = df[int(res)].tolist()

    counts = Counter(homolog_SCE)
    if counts.get("X", 0) == len(homolog_SCE):
        return np.nan

    counts.pop("X", None)
    
    most_common = counts.most_common()
    
    if not most_common:
        return np.nan  # No non-X values present

    max_count = most_common[0][1]
    tie = [val for val, count in most_common if count == max_count]

    if len(tie) > 1:
        return "X"

    return most_common[0][0]

def fetch_fasta(id_ch):
    seq = ''
    with open(f'{id_ch}.fasta','r') as infile:
        lines = infile.readlines()
        for line in lines[1:]:
            seq += line.strip()
    return seq


def main():
    feats = pd.DataFrame()
    id_ch = f'{pdbid}_{chain}'
    seq = fetch_fasta(id_ch)
    r_features = pd.read_csv(f'./{id_ch}_R_features.csv')
    with open(f'{pdbid}_{chain}_HRPF.csv','r') as infile:
        for line in infile.readlines():
            residue_feats = {} # Will hold all features for current residue
            spline = line.split()
            resnum = spline[0]
            condition = (r_features['pdbid'] == pdbid) & (r_features['chain'] == chain) & (r_features['labeling_pos'] == int(resnum))
            filtered_rows = r_features[condition].iloc[:,3:]
            row_dict = filtered_rows.iloc[0].to_dict()
            residue_feats.update(row_dict)

            neg_res_offset = -3
            for res in range(int(resnum) - 3, int(resnum)):
                res_prefix = f'RES_n{neg_res_offset * -1}'
                try:
                    neighbor_ident = res_dict[seq[res - 1]] # FASTA is zero-indexed
                except:
                    neighbor_ident = None
                if neighbor_ident is None:
                    encoded_res = le_res.transform(['X'])
                    residue_feats[f"{res_prefix}_IDENT"] = encoded_res[0]
                    residue_feats[f"{res_prefix}_IR"] = np.nan
                    sce = fetch_SCE(id_ch, str(res), NA = True)
                    residue_feats[f"{res_prefix}_SCE"] = sce

                    disopred = fetch_disopred(pdbid, chain, str(res), NA = True)
                    for i in range(len(disopred)):
                        residue_feats[f"{res_prefix}_DIS_{i}"] = disopred[i]
                    
                    psipred = fetch_psipred(pdbid, chain, str(res),NA = True)
                    for i in range(len(psipred)):
                        residue_feats[f"{res_prefix}_PSI_{i}"] = psipred[i]

                    pssm_score = fetch_PSSM_Score(pdbid, chain, str(res), NA = True)
                    for i in range(len(pssm_score)):
                        residue_feats[f"{res_prefix}_PSSM_{i}"] = pssm_score[i]
                    neg_res_offset += 1
                    
                    continue
                encoded = le_res.transform([neighbor_ident])
                residue_feats[f"{res_prefix}_IDENT"] = encoded[0]

                residue_feats[f"{res_prefix}_IR"] = float(IRs[neighbor_ident])

                sce = fetch_SCE(id_ch, str(res))
                if isinstance(sce, str):
                    encoded_sce = le_sce.transform([sce])
                    residue_feats[f"{res_prefix}_SCE"] = encoded_sce[0]
                else:
                    residue_feats[f"{res_prefix}_SCE"] = sce

                disopred = fetch_disopred(pdbid, chain, str(res))
                for i in range(len(disopred)):
                    residue_feats[f"{res_prefix}_DIS_{i}"] = disopred[i]
                psipred = fetch_psipred(pdbid, chain, str(res))
                for i in range(len(psipred)):
                    residue_feats[f"{res_prefix}_PSI_{i}"] = psipred[i]
                pssm_score = fetch_PSSM_Score(pdbid, chain, str(res))
                for i in range(len(pssm_score)):
                    residue_feats[f"{res_prefix}_PSSM_{i}"] = pssm_score[i]
                neg_res_offset += 1     

            encoded = le_res.transform([spline[1]])
            residue_feats["RES_0_IDENT"] = encoded[0]
            residue_feats["RES_0_MOD"] = float(spline[-1])
            residue_feats["RES_0_IR"] = float(IRs[spline[1]])
            residue_feats["RES_0_PF"] = float(IRs[spline[1]] / float(spline[-1]))
            residue_feats["RES_0_lnPF"] = np.log(float(IRs[spline[1]] / float(spline[-1])))
            homolog_SCE = fetch_SCE(id_ch, resnum)
            
            if isinstance(homolog_SCE, str):
                encoded_SCE = le_sce.transform([homolog_SCE])
                residue_feats["RES_0_SCE"] = encoded_SCE[0]
            else:
                residue_feats["RES_0_SCE"] = homolog_SCE

            disopred = fetch_disopred(pdbid, chain, resnum)
            for i in range(len(disopred)):
                    residue_feats[f"RES_0_DIS_{i}"] = disopred[i]

            psipred = fetch_psipred(pdbid, chain, resnum)
            for i in range(len(psipred)):
                    residue_feats[f"RES_0_PSI_{i}"] = psipred[i]

            pssm_score = fetch_PSSM_Score(pdbid, chain, resnum)
            for i in range(len(pssm_score)):
                residue_feats[f"RES_0_PSSM_{i}"] = pssm_score[i]
            
            pos_res_offset = 1
            for res in range(int(resnum)+1, int(resnum) + 4):
                res_prefix = f'RES_p{pos_res_offset}'
                try:
                    neighbor_ident = res_dict[seq[res - 1]] # FASTA is zero-indexed
                except:
                    neighbor_ident = None
                if neighbor_ident is None:
                    print('None')
                    encoded_res = le_res.transform(['X'])
                    residue_feats[f"{res_prefix}_IDENT"] = encoded_res[0]

                    sce = fetch_SCE(id_ch, str(res), NA = True)
                    residue_feats[f"{res_prefix}_SCE"] = sce

                    disopred = fetch_disopred(pdbid, chain, str(res), NA = True)
                    for i in range(len(disopred)):
                        residue_feats[f"{res_prefix}_DIS_{i}"] = disopred[i]

                    psipred = fetch_psipred(pdbid, chain, str(res),NA = True)
                    for i in range(len(psipred)):
                        residue_feats[f"{res_prefix}_PSI_{i}"] = psipred[i]

                    pssm_score = fetch_PSSM_Score(pdbid, chain, str(res), NA = True)
                    for i in range(len(pssm_score)):
                        residue_feats[f"{res_prefix}_PSSM_{i}"] = pssm_score[i]
                    pos_res_offset += 1
                    continue


                encoded = le_res.transform([neighbor_ident])
                residue_feats[f"{res_prefix}_IDENT"] = encoded[0]
                residue_feats[f"{res_prefix}_IR"] = float(IRs[neighbor_ident])

                sce = fetch_SCE(id_ch, str(res))
                if isinstance(sce, str):
                    encoded_sce = le_sce.transform([sce])
                    residue_feats[f"{res_prefix}_SCE"] = encoded_sce[0]
                else:
                    residue_feats[f"{res_prefix}_SCE"] = sce

                disopred = fetch_disopred(pdbid, chain, str(res))
                for i in range(len(disopred)):
                    residue_feats[f"{res_prefix}_DIS_{i}"] = disopred[i]
                psipred = fetch_psipred(pdbid, chain, str(res))
                for i in range(len(psipred)):
                    residue_feats[f"{res_prefix}_PSI_{i}"] = psipred[i]
                pssm_score = fetch_PSSM_Score(pdbid, chain, str(res))
                for i in range(len(pssm_score)):
                    residue_feats[f"{res_prefix}_PSSM_{i}"] = pssm_score[i]
                pos_res_offset += 1

            if feats.empty:
                feats = pd.DataFrame(columns = residue_feats.keys())
            feats.loc[len(feats)] = residue_feats
    feats.to_pickle(f'{id_ch}_input_features.pkl')

main()
