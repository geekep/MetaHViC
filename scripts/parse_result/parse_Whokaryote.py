import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, average_precision_score, recall_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix

def construct_euk(df):
    result = []
    euk = set(('eukaryote', 'organelle'))
    for _, row in df.iterrows():
        if row['predicted'] in euk: result.append(1)
        else: result.append(0)
    return np.array(result)

def construct_prok(df):
    result = []
    prok = set(('archaea', 'bacteria', 'prokaryote'))
    for _, row in df.iterrows():
        if row['predicted'] in prok: result.append(1)
        else: result.append(0)
    return np.array(result)


TARGET_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_groundtruth'
target_fn = list(filter(lambda f: f.endswith('.pkl') and not f.startswith('.'), os.listdir(TARGET_DIR)))
target_fn.sort()

RESULT_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_Whokaryote'
results_fn = list(filter(lambda f: not f.endswith('.fa') and not f.startswith('.'), os.listdir(RESULT_DIR)))
results_fn.sort()
results_fn = [f'{f}/whokaryote_predictions_S.tsv' for f in results_fn]
# results_fn = [f'{f}/whokaryote_predictions_T.tsv' for f in results_fn]

summary_df = pd.DataFrame(columns=['filename',
                                   'category',
                                   'accuracy',
                                   'accuracy_score',
                                   'balanced_accuracy_score',
                                   'precision',
                                   'precision_score',
                                   'recall',
                                   'recall_score',
                                   'f1_score',
                                   'matthews_corrcoef'])
categories = ['eukaryote']

mistakes = []
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), index_col=0)
    classified_idx = [seqname.split()[0] for seqname in list(result.index)]
    result = result.set_index(pd.Index(classified_idx))
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f[:-29]}.fa.gz.pkl'))
    target = np.array([target_pkl.loc[seqname] for seqname in classified_idx]).flatten()
    missed_count = len(target_pkl) - len(target)

    missed_label = []
    for seqname in target_pkl.index:
        if seqname not in result.index:
            missed_label.append(target_pkl.loc[seqname])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    target_binary = np.zeros(len(target))
    target_binary[target==0] = 1
    result_binary = np.concatenate((construct_euk(result), np.zeros(missed_count)))
    summary_df = pd.concat([summary_df, pd.DataFrame([['_'.join(re.findall(r'\d+', f)[1:]),
                                                       categories[0],
                                                       (result_binary==target_binary).sum()/len(target_binary),
                                                       accuracy_score(target_binary, result_binary),
                                                       balanced_accuracy_score(target_binary, result_binary),
                                                       ((target_binary==1) & (result_binary==1)).sum()/(result_binary==1).sum(),
                                                       precision_score(target_binary, result_binary, average='binary'),
                                                       ((target_binary==1) & (result_binary==1)).sum()/(target_binary==1).sum(),
                                                       recall_score(target_binary, result_binary, average='binary'),
                                                       f1_score(target_binary, result_binary, average='binary'),
                                                       matthews_corrcoef(target_binary, result_binary)]],
                                                     columns=['filename',
                                                              'category',
                                                              'accuracy',
                                                              'accuracy_score',
                                                              'balanced_accuracy_score',
                                                              'precision',
                                                              'precision_score',
                                                              'recall',
                                                              'recall_score',
                                                              'f1_score',
                                                              'matthews_corrcoef'])], axis=0)
    
    mistake = [(target[result_binary==1]==3).sum(), (target[result_binary==1]==4).sum(), (target[result_binary==1]==1).sum(), (target[result_binary==1]==2).sum(), (target_binary[result_binary==0]==0).sum()]
    mistakes.append(mistake)

misclassified = pd.DataFrame(mistakes, columns=['Prok->Euk', 'ProkVir->Euk', 'EukVir->Euk', 'Plas->Euk', 'Euk->Non-Euk'])
misclassified.to_csv('perf_summary/misclassified/Whokaryote.csv', index=False, header=True)
# misclassified.to_csv('perf_summary/misclassified/Whokaryot_plus_Tiara.csv', index=False, header=True)

summary_df.to_csv('perf_summary/Whokaryote.csv', index=False, header=True)
# summary_df.to_csv('perf_summary/Whokaryote_plus_Tiara.csv', index=False, header=True)
