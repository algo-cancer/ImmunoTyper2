# Gene Prefix Consistency Metrics

This repository contains JSON data that quantifies the accuracy of immunoglobulin and T-cell receptor gene calls based on their prefix consistency values. The analysis was performed on 40 samples from the 1000 Genomes Project, with ground truth provided by HPRC assemblies.

## What is Prefix Consistency?

Prefix consistency is a confidence metric developed to assess the accuracy of allele calls derived from ImmunoTyper2. It works by measuring the stability of allele presence across alternative solutions of the integer linear programming (ILP) model.

The metric is calculated by:
1. Obtaining the optimal ILP solution and its objective value
2. Generating alternative solutions by constraining the objective value to be progressively worse (2%, 4%, 6%, and 8%)
3. Sequentially evaluating allele presence starting from the optimal solution
4. Counting consecutive solution bands in which the allele is present

An allele's prefix consistency value is the count of consecutive solution bands in which it appears (starting from optimal), terminated at the first absence. Alleles only present in the optimal solution receive a value of 0, while those present across all solutions receive a value of 5.

Higher prefix consistency values indicate greater confidence in the allele call. Thresholds can be empirically determined to classify calls as high or low confidence.

## Exact Values vs. Thresholds: Understanding the Difference

### Exact Prefix Consistency Values

Exact prefix consistency values (in `prefix_consistency_metrics`) measure performance for gene calls that have *exactly* that specific value. For example, `"2"` refers to gene calls that are present in exactly the optimal solution plus the first two alternative solutions, but absent in the third alternative solution.

**Example:** 
If IGHV1-2 has the following exact metrics:
```
"2": {
  "value": 2,
  "num_samples": 5,
  "tp": 4,
  "fp": 1,
  "precision": 0.8
}
```
This means that of all IGHV1-2 calls with a prefix consistency value of *exactly* 2, 80% were correct (4 out of 5).

### Prefix Consistency Thresholds

Threshold metrics (in `prefix_consistency_threshold_metrics`) measure performance for gene calls that have a prefix consistency value *greater than or equal to* that number. For example, `"2"` refers to gene calls with prefix consistency values of 2, 3, 4, or 5.

**Example:**
If IGHV1-2 has the following threshold metrics:
```
"2": {
  "value": 2,
  "num_samples": 12,
  "tp": 10,
  "fp": 2,
  "precision": 0.8333
}
```
This means that of all IGHV1-2 calls with a prefix consistency value of 2 or higher, 83.33% were correct (10 out of 12).

### Practical Application

These two metrics serve different purposes:

1. **Exact metrics** tell you how gene calls at each specific consistency level perform, which is useful for understanding the distribution of accuracy across different confidence levels.

2. **Threshold metrics** tell you how gene calls would perform if you applied a particular cutoff value in your analysis, which is useful for determining optimal filtering thresholds.

For example, if you set a threshold of prefix_consistency ≥ 3, you would only use gene calls with values of 3, 4, or 5. The threshold metrics show you what precision (and other metrics) you could expect with such a filter.

## JSON Data Structure

The JSON file contains comprehensive metrics on how prefix consistency relates to accuracy for different gene types and individual genes:

### Top Level Structure

```json
{
  "gene_types": {
    "IGHV": { ... },
    "IGLV": { ... },
    "IGKV": { ... },
    "TRAV": { ... },
    "TRGV": { ... },
    "TRDV": { ... }
  }
}
```

### Gene Type Level

Each gene type contains two main components:

```json
"IGHV": {
  "genes": { ... },  // Individual gene metrics
  "summary_metrics": { ... }  // Aggregated statistics
}
```

#### Summary Metrics

Summary metrics provide aggregated statistics across all genes of this type:

```json
"summary_metrics": {
  "mean_ppv": 0.8845,  // Average precision across genes
  "median_ppv": 0.9231,  // Median precision
  "mean_sensitivity": 0.7654,  // Average recall/sensitivity
  "median_sensitivity": 0.8124,  // Median sensitivity
  "mean_f_beta": 0.8901,  // Average F-beta score (beta=0.5)
  "median_f_beta": 0.9111,  // Median F-beta score
  "mean_passing_proportion": 0.6432,  // Average proportion of calls passing optimal threshold
  "median_passing_proportion": 0.6789,  // Median proportion passing
  "threshold_distribution": {  // Count of genes with each optimal threshold
    "1": 5,
    "2": 12,
    "3": 8,
    "4": 3,
    "5": 1
  }
}
```

### Gene Level

Each gene contains metrics across different prefix consistency values:

```json
"IGHV1-2": {
  "prefix_consistency_metrics": { ... },  // Metrics for EXACT values
  "prefix_consistency_threshold_metrics": { ... },  // Metrics for THRESHOLD values
  "total_samples": 25,  // Total samples with this gene call
  "optimal_threshold": 3  // Best threshold value for accuracy
}
```

#### Prefix Consistency Metrics

These metrics show performance for calls with EXACTLY this prefix consistency value:

```json
"prefix_consistency_metrics": {
  "0": {
    "value": 0,  // The exact prefix consistency value
    "num_samples": 3,  // Number of samples with exactly this value
    "tp": 2,  // True positive count
    "fp": 1,  // False positive count
    "precision": 0.6667  // Proportion of calls that are correct (tp/(tp+fp))
  },
  "1": { ... },
  // values 2-5 follow same pattern
}
```

#### Prefix Consistency Threshold Metrics

These metrics show performance for calls with prefix consistency values GREATER THAN OR EQUAL TO this threshold:

```json
"prefix_consistency_threshold_metrics": {
  "0": {
    "value": 0,  // The threshold value
    "num_samples": 15,  // Number of samples ≥ this value
    "tp": 10,  // True positive count
    "fp": 5,  // False positive count
    "precision": 0.6667  // Proportion of calls that are correct (tp/(tp+fp))
  },
  "1": { ... },
  // thresholds 2-5 follow same pattern
}
```
