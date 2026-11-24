import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, average_precision_score, recall_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix

def construct_result(df, target):
    result = []
    for _, row in df.iterrows():
        if row['Possible_source'] == target: result.append(1)
        else: result.append(0)
    return np.array(result)


TARGET_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_groundtruth'
target_fn = list(filter(lambda f: f.endswith('.pkl') and not f.startswith('.'), os.listdir(TARGET_DIR)))
target_fn.sort()

RESULT_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PPR-Meta'
results_fn = list(filter(lambda f: f.endswith('.csv') and not f.startswith('.'), os.listdir(RESULT_DIR)))
results_fn.sort()

summary_df = pd.DataFrame(columns=['filename',
                                   'category',
                                   'accuracy',
                                   'accuracy_score',
                                   'balanced_accuracy_score',
                                   'precision',
                                   'precision_score',
                                   'average_precision_score',
                                   'recall',
                                   'recall_score',
                                   'f1_score',
                                   'matthews_corrcoef',
                                   'roc_auc_score'])
categories = ['chromosome', 'plasmid', 'phage', 'multiclass']
result_mapping = {'chromosome': 0, 'plasmid': 1, 'phage': 2}

mistakes_plasmid = []
mistakes_prokvirus = []
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), sep=',', header=0)
    result_2D = result.iloc[:, 2:5].to_numpy()
    classified_idx = [seqname.split()[0].strip('"') for seqname in list(result['Header'])]
    result = result.set_index(pd.Index(classified_idx))
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[:-4]}.fa.gz.pkl'))
    target = np.array([target_pkl.loc[seqname] for seqname in classified_idx]).flatten()
    missed_count = len(target_pkl) - len(target)
    result_2D = np.pad(result_2D,((0, max(0, missed_count)), (0, 0)), mode='constant', constant_values=0)
    
    missed_label = []
    for seqname in target_pkl.index:
        if seqname not in result.index:
            missed_label.append(target_pkl.loc[seqname])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    target_binary = np.zeros(len(target))
    target_binary[target==2] = 1
    result_binary = np.concatenate((construct_result(result, 'plasmid'), np.zeros(missed_count)))
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:]),
                                           categories[1],
                                           (result_binary==target_binary).sum()/len(target_binary),
                                           accuracy_score(target_binary, result_binary),
                                           balanced_accuracy_score(target_binary, result_binary),
                                           ((target_binary==1) & (result_binary==1)).sum()/(result_binary==1).sum(),
                                           precision_score(target_binary, result_binary, average='binary'),
                                           -1,
                                           ((target_binary==1) & (result_binary==1)).sum()/(target_binary==1).sum(),
                                           recall_score(target_binary, result_binary, average='binary'),
                                           f1_score(target_binary, result_binary, average='binary'),
                                           matthews_corrcoef(target_binary, result_binary),
                                           -1]],
                                           columns=['filename',
                                                    'category',
                                                    'accuracy',
                                                    'accuracy_score',
                                                    'balanced_accuracy_score',
                                                    'precision',
                                                    'precision_score',
                                                    'average_precision_score',
                                                    'recall',
                                                    'recall_score',
                                                    'f1_score',
                                                    'matthews_corrcoef',
                                                    'roc_auc_score'])])

    target_binary = np.zeros(len(target))
    target_binary[target==4] = 1
    result_binary = np.concatenate((construct_result(result, 'phage'), np.zeros(missed_count)))
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:]),
                                           categories[2],
                                           (result_binary==target_binary).sum()/len(target_binary),
                                           accuracy_score(target_binary, result_binary),
                                           balanced_accuracy_score(target_binary, result_binary),
                                           ((target_binary==1) & (result_binary==1)).sum()/(result_binary==1).sum(),
                                           precision_score(target_binary, result_binary, average='binary'),
                                           -1,
                                           ((target_binary==1) & (result_binary==1)).sum()/(target_binary==1).sum(),
                                           recall_score(target_binary, result_binary, average='binary'),
                                           f1_score(target_binary, result_binary, average='binary'),
                                           matthews_corrcoef(target_binary, result_binary),
                                           -1]],
                                           columns=['filename',
                                                    'category',
                                                    'accuracy',
                                                    'accuracy_score',
                                                    'balanced_accuracy_score',
                                                    'precision',
                                                    'precision_score',
                                                    'average_precision_score',
                                                    'recall',
                                                    'recall_score',
                                                    'f1_score',
                                                    'matthews_corrcoef',
                                                    'roc_auc_score'])])

    ppr_idx = np.logical_and(target!=0, target!=1)
    ppr_target = np.zeros(len(target))
    ppr_target[target==2] = 1
    ppr_target[target==4] = 2
    ppr_target = ppr_target[ppr_idx]
    result = np.array([result_mapping[i] for i in list(result['Possible_source'])])
    result = np.concatenate((result, np.zeros(missed_count)))
    ppr_result = result[ppr_idx]
    ppr_result_2D = result_2D[ppr_idx]
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:]),
                                           categories[3],
                                           (ppr_result==ppr_target).sum()/len(ppr_target),
                                           accuracy_score(ppr_target, ppr_result),
                                           balanced_accuracy_score(ppr_target, ppr_result),
                                           -1,
                                           precision_score(ppr_target, ppr_result, average='weighted'),
                                           average_precision_score(ppr_target, ppr_result_2D, average='weighted'),
                                           -1,
                                           recall_score(ppr_target, ppr_result, average='weighted'),
                                           f1_score(ppr_target, ppr_result, average='weighted'),
                                           matthews_corrcoef(ppr_target, ppr_result),
                                           roc_auc_score(ppr_target, ppr_result_2D, average='weighted', multi_class='ovr')]],
                                           columns=['filename',
                                                    'category',
                                                    'accuracy',
                                                    'accuracy_score',
                                                    'balanced_accuracy_score',
                                                    'precision',
                                                    'precision_score',
                                                    'average_precision_score',
                                                    'recall',
                                                    'recall_score',
                                                    'f1_score',
                                                    'matthews_corrcoef',
                                                    'roc_auc_score'])])

    # For plasmid
    mistake_plasmid = [(target[result==1]==3).sum(), (target[result==1]==4).sum(), (target[result==1]==0).sum(), (target[result==1]==1).sum(), (target[result!=1]==2).sum()]
    mistakes_plasmid.append(mistake_plasmid)
    # For prokaryotic virus
    mistake_prokvirus = [(target[result==2]==3).sum(), (target[result==2]==0).sum(), (target[result==2]==1).sum(), (target[result==2]==2).sum(), (target[result!=2]==4).sum()]
    mistakes_prokvirus.append(mistake_prokvirus)

# For plasmid
misclassified = pd.DataFrame(mistakes_plasmid, columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas'])
misclassified.to_csv('perf_summary/misclassified/PPR-Meta_plasmid.csv', index=False, header=True)
# For prokaryotic virus
misclassified = pd.DataFrame(mistakes_prokvirus, columns=['Prok->ProkVir', 'Euk->ProkVir', 'EukVir->ProkVir', 'Plas->ProkVir', 'ProkVir->NonProkVir'])
misclassified.to_csv('perf_summary/misclassified/PPR-Meta_prokvir.csv', index=False, header=True)

summary_df.to_csv('perf_summary/PPR-Meta.csv', index=False, header=True)
