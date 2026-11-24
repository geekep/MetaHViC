import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, average_precision_score, recall_score, f1_score, matthews_corrcoef, roc_auc_score, confusion_matrix, precision_recall_curve, roc_curve, PrecisionRecallDisplay, RocCurveDisplay


def construct_target(nums):
    order = [[0], [1], [2], [3], [4]]
    return np.array(sum([order[i] * nums[i] for i in range(len(order))], []))

def construct_array(array, class_idx):
    target_array = np.zeros((array.shape[0], 2))
    target_array[array!=class_idx, 0] = 1
    target_array[array==class_idx, 1] = 1
    return target_array

def construct_vir_array(array):
    vir_array = np.zeros((array.shape[0], 2))
    vir_array[np.logical_or.reduce((array==0, array==2, array==3)), 0] = 1
    vir_array[np.logical_or(array==1, array==4), 1] = 1
    return vir_array

def construct_prok_array(array):
    prok_array = np.zeros((array.shape[0], 2))
    prok_array[np.logical_or(array==3, array==4), 0] = 1
    prok_array[np.logical_or.reduce((array==0, array==1, array==2)), 1] = 1
    return prok_array

def construct_euk_array(array):
    euk_array = np.zeros((array.shape[0], 2))
    euk_array[np.logical_or.reduce((array==2, array==3, array==4)), 0] = 1
    euk_array[np.logical_or(array==0, array==1), 1] = 1
    return euk_array

TARGET_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_groundtruth'
target_fn = list(filter(lambda f: f.endswith('.pkl') and not f.startswith('.'), os.listdir(TARGET_DIR)))
target_fn.sort()

RESULT_DIR = '/data/DeepMicroClass/data/sequence_subsampled_1000_contig_6000_class_5_IGLOO-Stack_undersampling_nb_stacks_1_fungi'
results_fn = list(filter(lambda f: f.endswith('_host_classification') and not f.startswith('.'), os.listdir(RESULT_DIR)))
results_fn.sort()
results_fn = [f'{fn}/{fn}.tsv' for fn in results_fn]

summary_df = pd.DataFrame(columns=['filename', 'category', 'accuracy', 'accuracy_score', 'balanced_accuracy_score',
                                   'precision', 'precision_score', 'average_precision_score', 'recall', 'recall_score',
                                   'f1_score', 'matthews_corrcoef', 'roc_auc_score'])
categories = ['prokaryote', 'prokaryote_virus', 'plasmid', 'eukaryote', 'eukaryote_virus',
              'virus', 'prok_group', 'euk_group', 'multiclass']

mistakes =[]
for f in results_fn:
    result = pd.read_table(os.path.join(RESULT_DIR, f), sep='\t', header=0, index_col=0)
    result_2D = result
    result = result.to_numpy().argmax(axis=1)

    nums_str = re.findall(r'\d+', f.rsplit('/', 1)[0])[1:]
    nums = [int(n) for n in nums_str]
    target = construct_target(nums)

    # fig, axes = plt.subplots(nrows=len(nums_str)+1, ncols=2, figsize=(6, 10), sharex=True, sharey=True)
    for class_idx in range(len(nums_str)):
        target_class = construct_array(target, class_idx).argmax(axis=1)
        result_class = construct_array(result, class_idx).argmax(axis=1)
        summary_df = pd.concat([summary_df,
                                pd.DataFrame([['_'.join(nums_str),
                                               categories[class_idx],
                                               (target_class==result_class).sum()/len(target_class),
                                               accuracy_score(target_class, result_class),
                                               balanced_accuracy_score(target_class, result_class),
                                               ((target_class==1) & (result_class==1)).sum()/(result_class==1).sum(),
                                               precision_score(target_class, result_class, average='binary'),
                                               -1,
                                               ((target_class==1) & (result_class==1)).sum()/(target_class==1).sum(),
                                               recall_score(target_class, result_class, average='binary'),
                                               f1_score(target_class, result_class, average='binary'),
                                               matthews_corrcoef(target_class, result_class),
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

        # precision, recall, _ = precision_recall_curve(y_true=target, probas_pred=result_2D, pos_label=class_idx)
        # PrecisionRecallDisplay(precision=precision, recall=recall, ).plot(ax=axes[class_idx, 0])
        # fpr, tpr, _ = roc_curve(y_true=target, y_score=result_2D, pos_label=class_idx)
        # RocCurveDisplay(fpr=fpr, tpr=tpr, ).plot(ax=axes[class_idx, 1])

    target_vir = construct_vir_array(target).argmax(axis=1)
    result_vir = construct_vir_array(result).argmax(axis=1)
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(nums_str),
                                           categories[-4],
                                           (target_vir==result_vir).sum()/len(target_vir),
                                           accuracy_score(target_vir, result_vir),
                                           balanced_accuracy_score(target_vir, result_vir),
                                           ((target_vir==1) & (result_vir==1)).sum()/(result_vir==1).sum(),
                                           precision_score(target_vir, result_vir, average='binary'),
                                           -1,
                                           ((target_vir==1) & (result_vir==1)).sum()/(target_vir==1).sum(),
                                           recall_score(target_vir, result_vir, average='binary'),
                                           f1_score(target_vir, result_vir, average='binary'),
                                           matthews_corrcoef(target_vir, result_vir),
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
    
    target_prok = construct_prok_array(target).argmax(axis=1)
    result_prok = construct_prok_array(result).argmax(axis=1)
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(nums_str),
                                           categories[-3],
                                           (target_prok==result_prok).sum()/len(target_prok),
                                           accuracy_score(target_prok, result_prok),
                                           balanced_accuracy_score(target_prok, result_prok),
                                           ((target_prok==1) & (result_prok==1)).sum()/(result_prok==1).sum(),
                                           precision_score(target_prok, result_prok, average='binary'),
                                           -1,
                                           ((target_prok==1) & (result_prok==1)).sum()/(target_prok==1).sum(),
                                           recall_score(target_prok, result_prok, average='binary'),
                                           f1_score(target_prok, result_prok, average='binary'),
                                           matthews_corrcoef(target_prok, result_prok),
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
    
    target_euk = construct_euk_array(target).argmax(axis=1)
    result_euk = construct_euk_array(result).argmax(axis=1)
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(nums_str),
                                           categories[-2],
                                           (target_euk==result_euk).sum()/len(target_euk),
                                           accuracy_score(target_euk, result_euk),
                                           balanced_accuracy_score(target_euk, result_euk),
                                           ((target_euk==1) & (result_euk==1)).sum()/(result_euk==1).sum(),
                                           precision_score(target_euk, result_euk, average='binary'),
                                           -1,
                                           ((target_euk==1) & (result_euk==1)).sum()/(target_euk==1).sum(),
                                           recall_score(target_euk, result_euk, average='binary'),
                                           f1_score(target_euk, result_prok, average='binary'),
                                           matthews_corrcoef(target_euk, result_euk),
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
    
    summary_df = pd.concat([summary_df,
                            pd.DataFrame([['_'.join(nums_str),
                                           categories[-1],
                                           (result==target).sum()/len(target),
                                           accuracy_score(target, result),
                                           balanced_accuracy_score(target, result),
                                           -1,
                                           precision_score(target, result, average='weighted'),
                                           average_precision_score(target, result_2D, average='weighted'),
                                           -1,
                                           recall_score(target, result, average='weighted'),
                                           f1_score(target, result, average='weighted'),
                                           matthews_corrcoef(target, result),
                                           roc_auc_score(target, result_2D, average='weighted', multi_class='ovr')]],
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

    # precision, recall, _ = precision_recall_curve(y_true=target, probas_pred=result_2D)
    # PrecisionRecallDisplay(precision=precision, recall=recall, ).plot(ax=axes[len(nums_str), 0])
    # fpr, tpr, _ = roc_curve(y_true=target, y_score=result_2D)
    # RocCurveDisplay(fpr=fpr, tpr=tpr, ).plot(ax=axes[len(nums_str), 1])
    # fig.supxlabel('False Positive Rate', fontsize=16)
    # fig.supylabel('True Positive Rate', fontsize=16)
    # fig.tight_layout(pad=0.8)
    # plt.yticks(np.arange(0, 1.01, 0.2))
    # plt.savefig('plot_figures/plot_benchmark/analysis_results/ROC.pdf', bbox_inches='tight')

    conf_mat = confusion_matrix(target, result)
    conf_mat = conf_mat[~np.eye(conf_mat.shape[0], dtype=bool)].reshape(conf_mat.shape[0], -1)
    mistakes.append(conf_mat.flatten())

label = ['Prok', 'ProkVir', 'Plasmid', 'Euk', 'EukVir']
label = [f'{label[i]}->{label[j]}' for i in range(len(label)) for j in range(len(label)) if i!=j]
misclass_df = pd.DataFrame(data=mistakes, columns=label)
misclass_df.to_csv(f'perf_summary/misclassified/IGLOO-Stack_undersampling_nb_stacks_1_fungi.csv', index=False, header=True)

summary_df.to_csv('perf_summary/IGLOO-Stack_undersampling_nb_stacks_1_fungi.csv', index=False, header=True)
