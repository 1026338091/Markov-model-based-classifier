import matplotlib.pyplot as plt
import math
import os

# set
DATA_DIR = "/Users/denghongyuan/PycharmProjects/python-Project-my"

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


def scoreModels(spacii_file="Spacii.fa"):


    # path
    codingMatrix = getProbs(os.path.join(DATA_DIR, "codingModel.tab"))
    noncodingMatrix = getProbs(os.path.join(DATA_DIR, "noncodingModel.tab"))
    id2ancestorSeq = getSeq(os.path.join(DATA_DIR, "Ancestor.fa"))
    id2spaciiSeq = getSeq(os.path.join(DATA_DIR, spacii_file))

    allID = list(id2ancestorSeq.keys())
    results = []

    print("\nSequence Classification Results")
    print("=" * 110)
    print(
        f"{'Sequence ID':<20} {'True Label':<15} {'Prediction':<15} {'Coding Score':<15} {'Non-coding Score':<15} {'LLR':<10} {'Result'}")
    print("-" * 110)

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

        likelihood_ratio = cScore - nScore
        true_label = "Coding" if "_c_" in ID else "Non-coding"
        prediction = "Coding" if cScore > nScore else "Non-coding"
        correct = "✓" if prediction == true_label else "✗"

        results.append({
            'ID': ID,
            'true_label': true_label,
            'prediction': prediction,
            'cScore': cScore,
            'nScore': nScore,
            'LLR': likelihood_ratio,
            'correct': correct
        })

        print(
            f"{ID:<20} {true_label:<15} {prediction:<15} {cScore:<15.2f} {nScore:<15.2f} {likelihood_ratio:<10.2f} {correct}")


    correct_count = sum(1 for r in results if r['correct'] == "✓")
    total = len(results)

    print("\n" + "=" * 110)
    print(f"Accuracy: {correct_count}/{total} = {correct_count / total:.2%}")
    print("=" * 110)

    return results


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

    print("=" * 110)
    print("PART 1: Classification with Spacii.fa (divergence time t)")
    print("=" * 110)
    results_t = scoreModels("Spacii.fa")


    with open(os.path.join(DATA_DIR, "HongyuanDeng_Model.txt"), "w") as f:
        f.write("Coding Model Classification Results\n")
        f.write("=" * 110 + "\n\n")
        f.write("PART 1: Using Spacii.fa (divergence time t)\n")
        f.write("=" * 110 + "\n")
        f.write(
            f"{'Sequence ID':<20} {'True Label':<15} {'Prediction':<15} {'Coding Score':<15} {'Non-coding Score':<15} {'LLR':<10} {'Result'}\n")
        f.write("-" * 110 + "\n")

        for r in results_t:
            f.write(
                f"{r['ID']:<20} {r['true_label']:<15} {r['prediction']:<15} {r['cScore']:<15.2f} {r['nScore']:<15.2f} {r['LLR']:<10.2f} {r['correct']}\n")

        correct_count = sum(1 for r in results_t if r['correct'] == "✓")
        f.write("\n" + "=" * 110 + "\n")
        f.write(f"Accuracy: {correct_count}/{len(results_t)} = {correct_count / len(results_t):.2%}\n")


    print("\n\n" + "=" * 110)
    print("PART 2: Classification with Spacii_2100.fa (divergence time ~2t)")
    print("=" * 110)
    results_2t = scoreModels("Spacii_2100.fa")


    with open(os.path.join(DATA_DIR, "HongyuanDeng_Model.txt"), "a") as f:
        f.write("\n\n" + "=" * 110 + "\n")
        f.write("PART 2: Using Spacii_2100.fa (divergence time ~2t)\n")
        f.write("=" * 110 + "\n")
        f.write(
            f"{'Sequence ID':<20} {'True Label':<15} {'Prediction':<15} {'Coding Score':<15} {'Non-coding Score':<15} {'LLR':<10} {'Result'}\n")
        f.write("-" * 110 + "\n")

        for r in results_2t:
            f.write(
                f"{r['ID']:<20} {r['true_label']:<15} {r['prediction']:<15} {r['cScore']:<15.2f} {r['nScore']:<15.2f} {r['LLR']:<10.2f} {r['correct']}\n")

        correct_count = sum(1 for r in results_2t if r['correct'] == "✓")
        f.write("\n" + "=" * 110 + "\n")
        f.write(f"Accuracy: {correct_count}/{len(results_2t)} = {correct_count / len(results_2t):.2%}\n")

    print("\n✓ Results saved to HongyuanDeng_Model.txt")


    print("\n" + "=" * 110)
    print("Misclassified Sequences Analysis:")
    print("=" * 110)

    print("\nSpacii.fa misclassifications:")
    for r in results_t:
        if r['correct'] == "✗":
            print(f"  {r['ID']}: True={r['true_label']}, Predicted={r['prediction']}, LLR={r['LLR']:.2f}")

    print("\nSpacii_2100.fa misclassifications:")
    for r in results_2t:
        if r['correct'] == "✗":
            print(f"  {r['ID']}: True={r['true_label']}, Predicted={r['prediction']}, LLR={r['LLR']:.2f}")