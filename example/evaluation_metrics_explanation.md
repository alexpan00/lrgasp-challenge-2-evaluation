# How the evaluation metrics works
## Spearman Correlation Coefficient (SCC)
SCC evaluates the monotonic relationship between the estimation and gold standard, which is based on the rank for gene isoform abundance rather than the raw data. It is computed as<br>
![Eeq](https://latex-staging.easygenerator.com/gif.latex?SC%7B%7BC%7D_%7B%5CTheta%20%2C%5Chat%7B%5CTheta%20%7D%7D%7D%3D%5Cfrac%7B%5Coperatorname%7Bcov%7D%5Cleft%28%20r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%2Cr%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%20%5Cright%29%7D%7B%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%7D%7D%5Ccdot%20%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%7D%7D%7D)<br>
where ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D) and ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D) are the ranks of ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%5CTheta%20%7D) and ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%5Chat%7B%5CTheta%20%7D%7D), respectively, and ![formula](https://latex-staging.easygenerator.com/gif.latex?%7B%5Coperatorname%7Bcov%7D%5Cleft%28%20r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%2Cr%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%20%5Cright%29%7D) is the covariance of the rank variables,![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%7D%7D)  and ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7Bs%7D_%7Br%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%7D) are the sample coefficient of variations of  ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D) and ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D), respectively.<br>
![Spearman' rho](figures/spearman_correlation.png)<br>
Spearman Correlation Coefficient (SCC) between the estimation and gold standard. The SCC reveals gene CDY1 can be accurately quantified but not gene METTL9.
## Abundance Recovery Rate
Abundance Recovery Rate (ARR) represents the percentage of the estimated value to the real abundance, which is calculated by <br>
![Eq](https://latex.codecogs.com/gif.latex?AR%7B%7BR%7D_%7Bi%7D%7D%3D%5Cfrac%7B%7B%7B%7B%5Chat%7B%5Ctheta%20%7D%7D%7D_%7Bi%7D%7D%7D%7B%7B%7B%5Ctheta%20%7D_%7Bi%7D%7D%7D%5Ctimes%20100%25%2C%5Ctext%7B%20%7D%5Cleft%28%20i%3D1%2C2%2C%5Ccdots%20%2CI%20%5Cright%29)<br>
An accurate estimation should have the ARR value close to 100%.
## Median Relative Difference
Median Relative Difference (MRD) represents the median of relative difference for abundance estimates in all gene isoforms, which can be calculated by<br>
![Eq](https://latex.codecogs.com/gif.latex?MRD%3Dmedian%5Cleft%5C%7B%20%5Cfrac%7B%5Cleft%7C%20%7B%7B%5Ctheta%20%7D_%7Bi%7D%7D-%7B%7B%7B%5Chat%7B%5Ctheta%20%7D%7D%7D_%7Bi%7D%7D%20%5Cright%7C%7D%7B%7B%7B%5Ctheta%20%7D_%7Bi%7D%7D%7D%2C%5Ctext%7B%20%7D%5Cleft%28%20i%3D1%2C2%2C%5Ccdots%20%2CI%20%5Cright%29%20%5Cright%5C%7D)<br>
A small MRD value indicates good performance of abundance estimation.
## Normalized Root Mean Square Error
Normalized Root Mean Square Error (NRMSE) gives a measure of the extent to which the relationship deviates from a one-to-one linear relationship. It can be computed by<br>
![Eq](https://latex-staging.easygenerator.com/gif.latex?NRMSE%3D%5Cfrac%7B%5Csqrt%7B%5Cfrac%7B1%7D%7BI%7D%5Csum%5Climits_%7Bi%3D1%7D%5E%7BI%7D%7B%7B%7B%5Cleft%28%20%7B%7B%5Ctheta%20%7D_%7Bi%7D%7D-%7B%7B%7B%5Chat%7B%5Ctheta%20%7D%7D%7D_%7Bi%7D%7D%20%5Cright%29%7D%5E%7B2%7D%7D%7D%7D%7D%7B%7B%7Bs%7D_%7B%5CTheta%20%7D%7D%7D)<br>
where ![Eq](https://latex.codecogs.com/gif.latex?s_%7B%5CTheta%7D) is the sample coefficient of variation of ![Eq](https://latex.codecogs.com/gif.latex?%5CTheta). A good performance of abundance estimation should have a small value of NRMSE.
Since quantification performance could be influenced by the exon-isoform structure and the abundance, we evaluate the quantification of different sets of genes/transcripts with different isoform features, including isoform numbers, exon numbers, gold standard abundance values and a customized  statistic  K-value representing the complexity of exon-isoform structures.
## K-value
![K-value](figures/K-value.png)<br>
K-value is the condition number of the exon-isoform binary matrix, which can be used to measure the complexity of exon-isoform structures for each gene.
## Irreproducibility
![Irreproducibility](figures/Irreproducibility.png)<br>
By fitting the coefficient of variation versus average isoform abundance into a smooth curve, it can be shown that Method B has a lower coefficient of variation and higher Irreproducibility.
## Consistency
![Consistency](figures/consistency.png)<br>
By setting an expression threshold (e.g., 1 in this toy example), we can define which set of genes express (in blue) or not (in yellow). This statistic is to measure the consistency of the expressed gene sets between replicates. 
## Resolution entropy
![Resolution entropy](figures/resolution_entropy.png)<br>
(A) The software output only a few certain discrete values has lower resolution entropy as it cannot capture the continuous and subtle difference of gene expressions. (B) The software with continuous output values has higher resolution entropy.
## Fold-change-based evaluation
![Fold-change-based evaluation](figures/fold-change-based-evaluation.png)<br>
By counting the # of transcripts that are differentially expressed, statistics such as precision, recall, and accuracy can be defined and calculated.
## Evaluation with different gene features
![Evaluation with different gene features](figures/split.jpg)<br>
Explore the relationship between the above evaluation metrics and important features like K-value, expression level, isoform length, and number of exons by splitting the transcripts into groups on these features.
