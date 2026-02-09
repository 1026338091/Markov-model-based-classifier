import math
import os
import matplotlib.pyplot as plt
import numpy as np

os.chdir('/Users/denghongyuan/PycharmProjects/python-Project-my')

modelCodons = ['TTT', 'TTC', 'TTA', 'TTG', 'CTT', 'CTC', 'CTA',
               'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTT', 'GTC',
               'GTA', 'GTG', 'TCT', 'TCC', 'TCA', 'TCG', 'AGT',
               'AGC', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC',
               'ACA', 'ACG', 'GCT', 'GCC', 'GCA', 'GCG', 'TAT',
               'TAC', 'CAT', 'CAC', 'CAA', 'CAG', 'AAT', 'AAC',
               'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 'GAG', 'TGT',
               'TGC', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',
               'GGT', 'GGC', 'GGA', 'GGG', 'TGG', 'TAA', 'TAG',
               'TGA']


def scoreModels_for_ROC(spacii_file="Spacii_2100.fa"):


    codingMatrix = getProbs("codingModel.tab")
    noncodingMatrix = getProbs("noncodingModel.tab")
    id2ancestorSeq = getSeq("Ancestor.fa")
    id2spaciiSeq = getSeq(spacii_file)

    allID = list(id2ancestorSeq.keys())
    results = []

    for ID in allID:
        cScore = 0
        nScore = 0

        ancestorSeq = id2ancestorSeq[ID]
        spaciiSeq = id2spaciiSeq[ID]

        for i in range(0, min(len(ancestorSeq), len(spaciiSeq)) - 2, 3):
            ancestorCodon = ancestorSeq[i:i + 3]
            spaciiCodon = spaciiSeq[i:i + 3]

            if ancestorCodon in modelCodons and spaciiCodon in modelCodons:
                ancestorIndex = modelCodons.index(ancestorCodon)
                spaciiIndex = modelCodons.index(spaciiCodon)

                codingProb = codingMatrix[ancestorIndex][spaciiIndex]
                noncodingProb = noncodingMatrix[ancestorIndex][spaciiIndex]

                if codingProb > 0:
                    cScore += math.log(codingProb)
                if noncodingProb > 0:
                    nScore += math.log(noncodingProb)

        # Log likelihood ratio
        likelihood_ratio = cScore - nScore

        # 真实标签
        true_label = 1 if "_c_" in ID else 0

        results.append({
            'ID': ID,
            'true_label': true_label,
            'LLR': likelihood_ratio
        })

    return results


def calculate_ROC(results):



    sorted_results = sorted(results, key=lambda x: x['LLR'], reverse=True)


    n_positive = sum(1 for r in results if r['true_label'] == 1)
    n_negative = sum(1 for r in results if r['true_label'] == 0)

    print(f"\nDataset Statistics:")
    print(f"  Total sequences: {len(results)}")
    print(f"  Coding sequences (positive class): {n_positive}")
    print(f"  Non-coding sequences (negative class): {n_negative}")

    tpr_list = [0.0]
    fpr_list = [0.0]
    thresholds = [float('inf')]

    tp = 0
    fp = 0


    for r in sorted_results:
        if r['true_label'] == 1:
            tp += 1
        else:
            fp += 1

        tpr = tp / n_positive if n_positive > 0 else 0
        fpr = fp / n_negative if n_negative > 0 else 0

        tpr_list.append(tpr)
        fpr_list.append(fpr)
        thresholds.append(r['LLR'])


    tpr_list.append(1.0)
    fpr_list.append(1.0)
    thresholds.append(float('-inf'))


    auc = 0
    for i in range(len(fpr_list) - 1):
        auc += (fpr_list[i + 1] - fpr_list[i]) * (tpr_list[i] + tpr_list[i + 1]) / 2


    j_scores = [tpr_list[i] - fpr_list[i] for i in range(len(tpr_list))]
    best_idx = j_scores.index(max(j_scores))

    return fpr_list, tpr_list, thresholds, auc, best_idx


def plot_ROC(fpr_list, tpr_list, thresholds, auc, best_idx):


    best_threshold = thresholds[best_idx]
    best_tpr = tpr_list[best_idx]
    best_fpr = fpr_list[best_idx]
    best_specificity = 1 - best_fpr


    fig, ax = plt.subplots(figsize=(10, 9))


    ax.plot(fpr_list, tpr_list, 'b-', linewidth=3, label=f'ROC Curve (AUC = {auc:.3f})')


    ax.plot([0, 1], [0, 1], 'r--', linewidth=2, label='Random Classifier (AUC = 0.500)')


    ax.plot(best_fpr, best_tpr, 'go', markersize=15,
            label=f'Optimal Threshold = {best_threshold:.2f}', zorder=5)


    ax.annotate(f'Best Point\nSensitivity={best_tpr:.3f}\nSpecificity={best_specificity:.3f}',
                xy=(best_fpr, best_tpr), xytext=(best_fpr + 0.15, best_tpr - 0.15),
                fontsize=10, ha='left',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=2))


    ax.fill_between(fpr_list, tpr_list, alpha=0.2, color='blue', label='AUC Area')


    ax.set_xlabel('False Positive Rate (1 - Specificity)', fontsize=14, fontweight='bold')
    ax.set_ylabel('True Positive Rate (Sensitivity)', fontsize=14, fontweight='bold')
    ax.set_title('ROC Curve: Coding vs Non-coding Classification\nUsing Spacii_2100.fa (Divergence Time ~2t)',
                 fontsize=15, fontweight='bold')


    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.05])


    ax.legend(loc='lower right', fontsize=11, framealpha=0.9)


    ax.grid(True, alpha=0.3, linestyle='--')


    ax.set_xticks(np.arange(0, 1.1, 0.1))
    ax.set_yticks(np.arange(0, 1.1, 0.1))

    plt.tight_layout()
    plt.savefig('HongyuanDeng_ROC_Curve.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ ROC curve saved to: HongyuanDeng_ROC_Curve.png")
    plt.show()


    print("\n" + "=" * 80)
    print("ROC Analysis Results:")
    print("=" * 80)
    print(f"AUC (Area Under Curve): {auc:.4f}")
    print(f"\nOptimal Threshold (Youden's J maximization):")
    print(f"  Likelihood Ratio Cutoff: {best_threshold:.4f}")
    print(f"  Sensitivity (True Positive Rate): {best_tpr:.4f}")
    print(f"  Specificity (True Negative Rate): {best_specificity:.4f}")
    print(f"  Youden's J Index: {best_tpr - best_fpr:.4f}")
    print("=" * 80)


    print("\nClassification Performance at Various Thresholds:")
    print(f"{'Threshold':<15} {'Sensitivity':<15} {'Specificity':<15} {'Balanced Acc':<15}")
    print("-" * 60)


    for idx in [0, len(thresholds) // 4, len(thresholds) // 2, 3 * len(thresholds) // 4, best_idx, len(thresholds) - 1]:
        if 0 <= idx < len(thresholds):
            threshold = thresholds[idx]
            sensitivity = tpr_list[idx]
            specificity = 1 - fpr_list[idx]
            balanced_acc = (sensitivity + specificity) / 2

            marker = " ← OPTIMAL" if idx == best_idx else ""
            print(f"{threshold:<15.2f} {sensitivity:<15.3f} {specificity:<15.3f} {balanced_acc:<15.3f}{marker}")

    print("=" * 80)


def getProbs(f1):
    with open(f1) as f:
        pMatrix = []
        for line in f:
            tmp = line.rstrip().split("\t")
            tmp = [float(i) for i in tmp]
            pMatrix.append(tmp)
    return pMatrix


def getSeq(filename):
    with open(filename) as f:
        id2seq = {}
        currkey = ""
        for line in f:
            if line.startswith(">"):
                currkey = line[1:].split("|")[0].strip()
                id2seq[currkey] = ""
            else:
                id2seq[currkey] += line.rstrip()
    return id2seq


if __name__ == "__main__":
    print("=" * 80)
    print("ROC Curve Analysis")
    print("Dataset: Spacii_2100.fa (more divergent sequences)")
    print("=" * 80)


    results = scoreModels_for_ROC("Spacii_2100.fa")


    fpr_list, tpr_list, thresholds, auc, best_idx = calculate_ROC(results)


    plot_ROC(fpr_list, tpr_list, thresholds, auc, best_idx)
