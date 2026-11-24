import os
import re
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, average_precision_score, recall_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix

def construct_result(df):
    result = []
    for _, row in df.iterrows():
        if row['label'].startswith('plasmid'): result.append(1)
        else: result.append(0)
    return np.array(result)


TARGET_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_groundtruth'
target_fn = list(filter(lambda f: f.endswith('.pkl') and not f.startswith('.'), os.listdir(TARGET_DIR)))
target_fn.sort()

RESULT_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_PlasFlow'
results_fn = list(filter(lambda f: not f.endswith('.fasta') and not f.startswith('.'), os.listdir(RESULT_DIR)))
results_fn.sort()

summary_df = pd.DataFrame(columns=['filename', 'category', 'accuracy', 'accuracy_score', 'balanced_accuracy_score',
                                   'precision', 'precision_score', 'recall', 'recall_score',
                                   'f1_score', 'matthews_corrcoef'])
categories = ['plasmid']

mistakes = []
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), index_col=2)
    classified_idx = [seqname.split()[0] for seqname in list(result.index)]
    result = result.set_index(pd.Index(classified_idx))
    
    target_pkl = pd.read_pickle(os.path.join(TARGET_DIR, f'{f}.fa.gz.pkl'))
    target = np.array([target_pkl.loc[seqname] for seqname in classified_idx]).flatten()
    missed_count = len(target_pkl) - len(target)

    missed_label = []
    for seqname in target_pkl.index:
        if seqname not in result.index:
            missed_label.append(target_pkl.loc[seqname])
    missed_label = np.array(missed_label).flatten()
    target = np.concatenate((target, missed_label))

    target_binary = np.ones(len(target))
    target_binary[target!=2] = 0
    result_binary = np.concatenate((construct_result(result), np.zeros(missed_count)))
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
                                                                'matthews_corrcoef'])])

    mistake = [(target[result_binary==1]==3).sum(), (target[result_binary==1]==4).sum(), (target[result_binary==1]==0).sum(), (target[result_binary==1]==1).sum(), (target[result_binary==0]==2).sum()]
    mistakes.append(mistake)

misclass_df = pd.DataFrame(mistakes, columns=['Prok->Plas', 'ProkVir->Plas', 'Euk->Plas', 'EukVir->Plas', 'Plas->Non-Plas'])
misclass_df.to_csv('perf_summary/misclassified/PlasFlow.csv', index=False, header=True)

summary_df.to_csv('perf_summary/PlasFlow.csv', index=False, header=True)
