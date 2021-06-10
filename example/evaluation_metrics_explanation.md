# How the evaluation metrics works
## Spearman Correlation Coefficient (SCC)
SCC evaluates the monotonic relationship between the estimation and gold standard, which is based on the rank for gene isoform abundance rather than the raw data. It is computed as<br>
![Eeq](https://latex-staging.easygenerator.com/gif.latex?SC%7B%7BC%7D_%7B%5CTheta%20%2C%5Chat%7B%5CTheta%20%7D%7D%7D%3D%5Cfrac%7B%5Coperatorname%7Bcov%7D%5Cleft%28%20r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%2Cr%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%20%5Cright%29%7D%7B%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%7D%7D%5Ccdot%20%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%7D%7D%7D)<br>
where ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D) and ![Eq](https://latex-staging.easygenerator.com/gif.latex?r%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D) are the ranks of ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%5CTheta%20%7D) and ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%5Chat%7B%5CTheta%20%7D%7D), respectively, and ![formula](https://latex-staging.easygenerator.com/gif.latex?%7B%5Coperatorname%7Bcov%7D%5Cleft%28%20r%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%2Cr%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%20%5Cright%29%7D) is the covariance of the rank variables,![Eq](https://latex-staging.easygenerator.com/gif.latex?%7B%7Bs%7D_%7Br%7B%7Bg%7D_%7B%5CTheta%20%7D%7D%7D%7D)  and ![Eq](https://latex-staging.easygenerator.com/gif.latex?%7Bs%7D_%7Br%7B%7Bg%7D_%7B%7B%5Chat%7B%5CTheta%20%7D%7D%7D%7D%7D) are the sample standard deviations of  and , respectively.<br>
![Spearman' rho](figures/spearman_correlation.png)<br>
Spearman Correlation Coefficient (SCC) between the estimation and gold standard. The SCC reveals gene CDY1 can be accurately quantified but not gene METTL9.
## Reproducibility
![Reproducibility](figures/reproducibility.png)<br>
By fitting the standard deviation versus average isoform abundance into a smooth curve, it can be shown that Method B has a lower standard deviation and higher reproducibility.
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
