import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, average_precision_score, recall_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix


TARGET_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_groundtruth'
target_fn = list(filter(lambda f: f.endswith('.pkl') and not f.startswith('.'), os.listdir(TARGET_DIR)))
target_fn.sort()

RESULT_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_genomad'
results_fn = list(filter(lambda f: f.endswith('_aggregated_classification') and not f.startswith('.'), os.listdir(RESULT_DIR)))
results_fn.sort()
results_fn = [f'{fn}/{fn}.tsv' for fn in results_fn]

# results_fn = list(filter(lambda f: f.endswith('_score_calibration') and not f.startswith('.'), os.listdir(RESULT_DIR)))
# results_fn.sort()
# results_fn = [f'{fn}/{fn.rsplit(sep='_', maxsplit=2)[0]}_calibrated_aggregated_classification.tsv' for fn in results_fn]

# results_fn = list(filter(lambda f: f.endswith('_nn_classification') and not f.startswith('.'), os.listdir(RESULT_DIR)))
# results_fn.sort()
# results_fn = [f'{fn}/{fn}.tsv' for fn in results_fn]

# results_fn = list(filter(lambda f: f.endswith('_marker_classification') and not f.startswith('.'), os.listdir(RESULT_DIR)))
# results_fn.sort()
# results_fn = [f'{fn}/{fn}.tsv' for fn in results_fn]

summary_df = pd.DataFrame(columns=['filename', 'category', 'accuracy', 'accuracy_score', 'balanced_accuracy_score',
                                   'precision', 'precision_score', 'average_precision_score', 'recall', 'recall_score',
                                   'f1_score', 'matthews_corrcoef', 'roc_auc_score'])
categories = ['chromosome', 'plasmid', 'virus', 'multiclass']

mistakes_plasmid =[]
mistakes_virus = []
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), sep='\t', header=0, index_col=0)
    result_2D = result
    classified_idx = list(result.index)
    result = result.to_numpy().argmax(axis=1)
    
    target_filename = [fn for fn in target_fn if fn.startswith(f[:2])][0]
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, target_filename))
    target = np.array([target_pkl.loc[seqname] for seqname in classified_idx]).flatten()

    result_binary = np.zeros(len(result))
    result_binary[result==0] = 1
    target_binary = np.zeros(len(target))
    target_binary[np.logical_or(target==3, target==0)] = 1
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:6]),
                                           categories[0],
                                           (target_binary==result_binary).sum()/len(target_binary),
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

    result_binary = np.zeros(len(result))
    result_binary[result==1] = 1
    target_binary = np.zeros(len(target))
    target_binary[target==2] = 1
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:6]),
                                           categories[1],
                                           (target_binary==result_binary).sum()/len(target_binary),
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

    result_binary = np.zeros(len(result))
    result_binary[result==2] = 1
    target_binary = np.zeros(len(target))
    target_binary[np.logical_or(target==4, target==1)] = 1
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:6]),
                                           categories[2],
                                           (target_binary==result_binary).sum()/len(target_binary),
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

    target_multiclass = np.zeros(len(target))
    target_multiclass[np.logical_or(target==0, target==3)] = 0
    target_multiclass[target==2] = 1
    target_multiclass[np.logical_or(target==1, target==4)] = 2
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:6]),
                                           categories[3],
                                           (target_multiclass==result).sum()/len(target_binary),
                                           accuracy_score(target_multiclass, result),
                                           balanced_accuracy_score(target_multiclass, result),
                                           -1,
                                           precision_score(target_multiclass, result, average='weighted'),
                                           average_precision_score(target_multiclass, result_2D, average='weighted'),
                                           -1,
                                           recall_score(target_multiclass, result_binary, average='weighted'),
                                           f1_score(target_multiclass, result_binary, average='weighted'),
                                           matthews_corrcoef(target_multiclass, result_binary),
                                           roc_auc_score(target_multiclass, result_2D, average='weighted', multi_class='ovr')]],
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
    # For virus
    mistake_virus = [(target[result==2]==3).sum(), (target[result==2]==0).sum(), (target[result==2]==2).sum(), (target[result!=2]==4).sum()+(target[result!=2]==1).sum()]
    mistakes_virus.append(mistake_virus)

# For plasmid
misclassified = pd.DataFrame(mistakes_plasmid, columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->NonPlas'])
misclassified.to_csv('perf_summary/misclassified/genomad_plasmid.csv', index=False, header=True)
# For virus
misclassified = pd.DataFrame(mistakes_virus, columns=['Prok->Vir', 'Euk->Vir', 'Plas->Vir', 'Vir->NonVir'])
misclassified.to_csv('perf_summary/misclassified/genomad_virus.csv', index=False, header=True)

summary_df.to_csv('perf_summary/genomad.csv', index=False, header=True)
