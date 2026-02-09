# Coding Sequence Classification Using Markov Models

## Project Overview
This project implements a **Markov model-based classifier** to distinguish coding (protein-coding) sequences from non-coding sequences by analyzing codon substitution patterns during evolution. The classifier uses log-likelihood ratios to determine whether evolutionary constraints characteristic of coding regions are present in evolved sequences.

## Author
Hongyuan Deng  
Course: Bioinformatics  
Institution: Northeastern University

## Biological Background

### The Problem
When genes evolve over time, some maintain their protein-coding function while others lose it and become pseudogenes. This project addresses the question: **Can we identify which sequences are still functional by analyzing their mutation patterns?**

### Key Concepts

**Coding Sequences (Genes)**
- Encode functional proteins
- Under purifying selection pressure
- Prefer synonymous mutations (silent changes that don't alter amino acids)
- Example: TTT → TTC (both code for Phenylalanine)

**Non-coding Sequences (Pseudogenes)**
- Lost protein-coding function
- No selection pressure
- Accumulate random mutations
- Both synonymous and non-synonymous mutations occur freely

**Codon Substitution Models**
- **Coding Model**: Transition probabilities for codons under selection
- **Non-coding Model**: Transition probabilities for random codon changes
- Models trained from real genomic data

## Project Structure

```
.
├── HongyuanDeng_Model.py          # Main classification script
├── HongyuanDeng_ROC.py            # ROC curve analysis
├── HongyuanDeng.py                # UPGMA phylogenetic tree
├── codingModelClassifier.py       # Template code
├── upgma.py                       # UPGMA template
├── Ancestor.fa                    # Ancestral sequences (40 sequences)
├── Spacii.fa                      # Evolved sequences (time t)
├── Spacii_2100.fa                 # More diverged sequences (time ~2t)
├── codingModel.tab                # 64×64 coding transition matrix
├── noncodingModel.tab             # 64×64 non-coding transition matrix
├── HongyuanDeng_Model.txt         # Classification results
├── HongyuanDeng_UPGMA.txt         # Phylogenetic tree output
└── HongyuanDeng_ROC_Curve.png     # ROC curve visualization
```

## Methodology

### Part 1: Markov Model Classification

**Algorithm:**
1. Extract codon pairs from ancestral and evolved sequences
2. For each codon transition (Ancestor → Spacii):
   - Calculate probability under coding model
   - Calculate probability under non-coding model
3. Compute log-likelihood scores:
   - `cScore = Σ log P(transition | coding model)`
   - `nScore = Σ log P(transition | non-coding model)`
4. Calculate Log-Likelihood Ratio (LLR):
   - `LLR = cScore - nScore`
5. Classification:
   - `LLR > 0` → Coding sequence (coding model fits better)
   - `LLR < 0` → Non-coding sequence (non-coding model fits better)

**Mathematical Foundation:**
```
LLR = log[P(D|Coding) / P(D|Non-coding)]
    = log P(D|Coding) - log P(D|Non-coding)
    = Σ log P(codon_i → codon_j | Coding) - Σ log P(codon_i → codon_j | Non-coding)
```

### Part 2: ROC Curve Analysis

**Purpose:** Evaluate classifier performance across different thresholds

**Metrics:**
- **TPR (Sensitivity)**: Proportion of coding sequences correctly identified
- **FPR (1 - Specificity)**: Proportion of non-coding sequences incorrectly classified
- **AUC**: Area Under the ROC Curve (overall classifier quality)
- **Optimal Threshold**: Determined by maximizing Youden's J index

### Part 3: UPGMA Phylogenetic Tree

**Purpose:** Construct evolutionary relationships between species

**Algorithm:**
- Unweighted Pair Group Method with Arithmetic Mean (UPGMA)
- Hierarchical clustering based on genetic distances
- Output: Newick format phylogenetic tree

## Results

### Classification Performance

#### Dataset 1: Spacii.fa (Divergence time t)
```
Accuracy: 40/40 = 100%
```
- **Perfect classification!**
- All 20 coding sequences correctly identified (LLR > 0)
- All 20 non-coding sequences correctly identified (LLR < 0)

**Sample Results:**
| Sequence ID | True Label | Prediction | cScore | nScore | LLR | Result |
|-------------|------------|------------|--------|---------|-----|---------|
| SimSeq_c_1 | Coding | Coding | -76.91 | -82.05 | 5.14 | ✓ |
| SimSeq_n_2 | Non-coding | Non-coding | -248.07 | -176.80 | -71.27 | ✓ |
| SimSeq_c_14 | Coding | Coding | -54.93 | -77.06 | 22.13 | ✓ |

#### Dataset 2: Spacii_2100.fa (Divergence time ~2t)
```
Accuracy: 21/40 = 52.5%
```
- **Performance degradation with increased divergence time**
- Longer evolutionary time obscures selection signals
- Many coding sequences misclassified as non-coding

**Key Insight:** Classifier works well for moderately diverged sequences but loses power as evolutionary distance increases.

### ROC Curve Analysis (Spacii_2100.fa)

**Performance Metrics:**
- **AUC = 0.947** (Excellent discrimination ability)
- **Optimal Threshold = -108.89**
- **Sensitivity = 1.000** (100% of coding sequences detected)
- **Specificity = 0.810** (81% of non-coding correctly excluded)

**Interpretation:**
- AUC close to 1.0 indicates strong classification power
- Even with increased divergence, model retains good discriminative ability
- Trade-off between sensitivity and specificity can be adjusted via threshold

### UPGMA Phylogenetic Tree

**Input Species:**
- M_Spacii, T_Pain, G_Unit, Q_Doba, R_Mani, A_Finch

**Output (Newick format):**
```
(M_Spacii:6.9,(((T_Pain:1.0,G_Unit:1.0):2.0,Q_Doba:3.0):1.1,(R_Mani:2.0,A_Finch:2.0):2.1):2.8);
```

**Biological Interpretation:**
- T_Pain and G_Unit are most closely related (distance = 2)
- R_Mani and A_Finch form a sister group (distance = 4)
- M_Spacii is the most divergent species

## Implementation Details

### Key Functions

#### `scoreModels(spacii_file)`
Classifies sequences using Markov model comparison.

**Parameters:**
- `spacii_file`: Path to evolved sequences (Spacii.fa or Spacii_2100.fa)

**Returns:**
- List of classification results with scores and predictions

**Algorithm:**
1. Load probability matrices and sequences
2. For each sequence pair (ancestor-spacii):
   - Extract codons (triplets)
   - Look up transition probabilities in matrices
   - Sum log probabilities
3. Compare scores and classify

#### `calculate_ROC(results)`
Computes ROC curve metrics from classification results.

**Returns:**
- FPR list, TPR list, thresholds, AUC, optimal threshold index

#### `plot_ROC()`
Generates publication-quality ROC curve visualization.

#### `UPGMA(distanceMatrix, speciesList)`
Implements UPGMA algorithm for phylogenetic tree construction.

**Key Steps:**
1. Find closest pair in distance matrix
2. Merge pair and update matrix
3. Calculate new distances (average of merged clusters)
4. Repeat until one cluster remains

## Usage

### Running the Classification
```python
python HongyuanDeng_Model.py
```

**Output:**
- Console: Classification results table
- File: `HongyuanDeng_Model.txt` (detailed results)

### Running ROC Analysis
```python
python HongyuanDeng_ROC.py
```

**Output:**
- Console: ROC metrics and optimal threshold
- File: `HongyuanDeng_ROC_Curve.png` (ROC plot)

### Running UPGMA
```python
python HongyuanDeng.py
```

**Output:**
- Console: Step-by-step tree construction
- File: `HongyuanDeng_UPGMA.txt` (Newick tree)

## Dependencies
```python
import math
import matplotlib.pyplot as plt
import numpy as np
```

## File Formats

### FASTA Files
```
>SimSeq_c_1
ATGTTAATATGGGACCAAGC...
>SimSeq_n_2
ATGTACAGATTGCTTCGGTC...
```
- Sequence IDs contain labels: `_c_` = coding, `_n_` = non-coding

### Probability Matrices (.tab)
- Tab-delimited 64×64 matrices
- Rows: Ancestral codons
- Columns: Evolved codons
- Values: Transition probabilities

### Distance Matrix
- Symmetric matrix of pairwise evolutionary distances
- Used for UPGMA clustering

## Key Findings

### 1. Selection Pressure Detection
**Coding sequences show distinctive patterns:**
- Higher probability of synonymous substitutions
- Lower probability of non-synonymous changes
- Detectable signal even after moderate divergence

### 2. Time-Dependent Signal Decay
**Classification accuracy decreases with divergence time:**
- Time t: 100% accuracy
- Time ~2t: 52.5% accuracy
- **Biological insight:** Selection signals weaken over evolutionary time

### 3. Model Performance
**ROC analysis reveals:**
- Strong discriminative power (AUC = 0.947)
- Can achieve 100% sensitivity with 81% specificity
- Optimal for comparative genomics within moderate evolutionary distances

## Biological Significance

### Why Do Coding Sequences Show Different Patterns?

**Genetic Code Redundancy:**
- 61 codons encode 20 amino acids
- Multiple codons per amino acid (synonymous codons)
- Example: Leucine encoded by 6 codons (TTA, TTG, CTT, CTC, CTA, CTG)

**Selection Pressure:**
- Coding sequences must maintain protein function
- Non-synonymous mutations often deleterious → removed by selection
- Synonymous mutations neutral → accumulate freely
- Result: Biased codon substitution patterns

**Pseudogenes:**
- Genes that lost function (e.g., through frameshift or stop codon)
- No longer under selection
- Accumulate mutations randomly
- Codon usage becomes random

## Limitations and Considerations

### 1. Divergence Time Sensitivity
- Classifier performs best on moderately diverged sequences
- Very ancient divergences overwhelm selection signals
- Optimal for comparisons within ~100-200 million years

### 2. Model Assumptions
- Assumes constant mutation rates
- Assumes similar evolutionary processes across lineages
- May not account for recent adaptive evolution

### 3. Data Requirements
- Requires known ancestral sequences
- Needs training data to build probability matrices
- Performance depends on model quality

## Extensions and Future Work

### Potential Improvements
1. **Higher-order Markov models**: Consider longer codon contexts
2. **Position-specific models**: Account for codon position biases
3. **Phylogenetic context**: Incorporate tree structure in classification
4. **Machine learning integration**: Combine with modern ML approaches

### Additional Analyses
- Gene ontology enrichment of retained vs. lost genes
- Chromosomal distribution of pseudogenes
- Correlation with expression levels
- Tissue-specific selection patterns

## References

### Key Concepts
- **Markov Models in Genomics**: Hidden Markov Models for gene finding
- **Molecular Evolution**: Neutral theory and selection
- **Phylogenetics**: UPGMA and distance-based methods

### Related Tools
- **PAML**: Phylogenetic Analysis by Maximum Likelihood
- **HyPhy**: Hypothesis testing using phylogenies
- **MCDP**: Markov Chain for genomic interval analysis

## Interpretation Guide

### Understanding LLR Values

**Strong Coding Signal (LLR >> 0):**
- Example: SimSeq_c_14 (LLR = 22.13)
- High confidence the sequence is under selection
- Likely encodes a functional protein

**Strong Non-coding Signal (LLR << 0):**
- Example: SimSeq_n_2 (LLR = -71.27)
- High confidence the sequence is not under selection
- Likely a pseudogene or non-functional region

**Ambiguous Cases (LLR ≈ 0):**
- Difficult to classify
- May be weakly selected or recently lost function
- Require additional evidence

### ROC Curve Interpretation

**AUC = 0.947:**
- 94.7% probability that a randomly chosen coding sequence ranks higher than a non-coding one
- Excellent classifier performance
- Much better than random (AUC = 0.5)

**Optimal Threshold (-108.89):**
- Maximizes balance between sensitivity and specificity
- At this threshold: 100% sensitivity, 81% specificity
- Can be adjusted based on application needs

## Conclusions

### Main Findings
1. **Markov models successfully detect coding constraints** in moderately diverged sequences (100% accuracy at time t)
2. **Signal degrades with evolutionary time** (52.5% accuracy at time ~2t)
3. **ROC analysis shows robust performance** (AUC = 0.947) even for highly diverged sequences
4. **UPGMA successfully reconstructs** species relationships from distance data

### Biological Insights
- Selection pressure on coding sequences creates detectable patterns in codon evolution
- These patterns persist for moderate evolutionary timescales
- Beyond ~2x divergence time, neutral drift obscures functional signals
- Combination of multiple evidence types improves classification

## Acknowledgments
- Course instructors for project design
- Provided probability matrices based on empirical genomic data
- Template code structure

## License
Educational use only - Course Project

---

*Last Updated: January 2026*
