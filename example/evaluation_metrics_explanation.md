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
where ![Eq](https://latex.codecogs.com/svg.image?s_%7B%5CTheta%7D) is the sample coefficient of variation of ![Eq](https://latex.codecogs.com/gif.latex?%5CTheta). A good performance of abundance estimation should have a small value of NRMSE.
Since quantification performance could be influenced by the exon-isoform structure and the abundance, we evaluate the quantification of different sets of genes/transcripts with different isoform features, including isoform numbers, exon numbers, gold standard abundance values and a customized  statistic  K-value representing the complexity of exon-isoform structures.
## Irreproducibility
The irreproducibility statistic characterizes the average coefficient of variation of abundance estimates among different replicates, which is calculated by
![Eq](https://latex.codecogs.com/svg.image?IM%3D%5Csqrt%7B%5Cfrac%7B1%7D%7BIG%7D%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5Csum_%7Bg%3D1%7D%5E%7BG%7D%20CV_%7Big%7D%5E%7B2%7D%20%20%7D%20#0)<br>
Here, ![Eq](https://latex.codecogs.com/svg.image?CV_%7Big%7D#0) is the coefficient of variation of ![Eq](https://latex.codecogs.com/svg.image?log%5Cleft%20(%20%5Chat%7B%5Ctheta%7D_%7Bigr%7D%20%2B1%20%20%5Cright%20)%20%5Cleft%20(%20r%3D1%2C2%2C%5Ccdots%20%2CR%20%5Cright%20)%20#0), which is calculated by ![Eq](https://latex.codecogs.com/svg.image?CV_%7Big%7D%3D%5Cfrac%7Bs_%7Big%7D%7D%7Bu_%7Big%7D%7D%20#0), where ![](https://latex.codecogs.com/svg.image?s_%7Big%7D#0) and ![](https://latex.codecogs.com/svg.image?u_%7Big%7D#0) are the sample standard deviation and mean of  abundance estimates, 
![Irreproducibility](figures/reproducibility.png)<br>

By fitting the coefficient of variation versus average isoform abundance into a smooth curve, it can be shown that Method B has a lower coefficient of variation and higher reproducibility.
## Consistency
Consistency measure examines the similarity of abundance profiles between mutual pairs of replicates, which is defined as:
![](https://latex.codecogs.com/svg.image?C%5Cleft%20(%20%5Calpha%20%20%5Cright%20)%20%3D%5Cfrac%7B1%7D%7BIG%5Ccdot%20C_%7BR%7D%5E%7B2%7D%20%7D%20%5Csum_%7Bi%3D1%7D%5E%7BI%7D%20%5Csum_%7Bg%3D1%7D%5E%7BG%7D%20%5Csum_%7B1%5Cle%20r_1%3Cr_2%5Cle%20R%7D%5E%7B%7D%20P%5Cleft%20(%20%5Cleft%20%5C%7B%20log%5Cleft%20(%20%5Chat%7B%5Ctheta%7D_%7Bigr_1%7D%20%2B1%20%5Cright%20)%20%3C%20%5Calpha%2C%20log%5Cleft%20(%20%5Chat%7B%5Ctheta%7D_%7Bigr_2%7D%20%2B1%20%5Cright%20)%20%3C%20%5Calpha%5Cright%20%5C%7D%20or%20%5Cleft%20%5C%7B%20log%5Cleft%20(%20%5Chat%7B%5Ctheta%7D_%7Bigr_1%7D%20%2B1%20%5Cright%20)%20%5Cge%20%20%5Calpha%2C%20log%5Cleft%20(%20%5Chat%7B%5Ctheta%7D_%7Bigr_2%7D%20%2B1%20%5Cright%20)%20%5Cge%20%20%5Calpha%5Cright%20%5C%7D%20%5Cright%20)%20%20#0), where ![](https://latex.codecogs.com/svg.image?%5Calpha#0) is a customized threshold defining whether a transcript is expressed or not. 
![Consistency](figures/consistency.png)<br>
By setting an expression threshold (e.g., 1 in this toy example), we can define which set of genes express (in blue) or not (in yellow). This statistic is to measure the consistency of the expressed gene sets between replicates. 
## Resolution entropy
For a given sample, a Resolution Entropy (RE) statistic characterizes the resolution of abundance estimation:![](https://latex.codecogs.com/svg.image?RE%3D-%5Csum_%7Bm%3D1%7D%5E%7BM%7DP_m%20In%5Cleft%20(%20P_m%20%5Cright%20)%2C%20%5Ctext%7Bwhere%20%7D%20P_m%20%3D%20%5Cfrac%7Bn_m%7D%7B%5Csum_%7Bj%3D1%7D%5E%7BM%7D%20n_j%20%7D.#0) <br>
Here, the abundance estimates are binned into M groups, where ![](https://latex.codecogs.com/svg.image?n_m#0) represents the number of transcript isoforms with the abundance estimate ![](https://latex.codecogs.com/svg.image?%5Cwidehat%7B%5CTheta%7D%20%5Cin%20%5Bm%5Ccdot%20%5Calpha%20%2C%20%5Cleft%20(%20m%2B1%20%20%5Cright%20)%5Ccdot%20%5Calpha%29%20#0), and ![](https://latex.codecogs.com/svg.image?%5Calpha%3Dmax%5Cleft%20(%20%5Cwidehat%7B%5CTheta%7D%20%20%5Cright%20)%20%2FM#0).![](https://latex.codecogs.com/svg.image?RE%3D0#0) if all transcript isoforms have the same estimated abundance values, while it obtains a large value when the estimates are uniformly distributed among M groups. <br>
![Resolution entropy](figures/resolution_entropy.png)<br>
(A) The software output only a few certain discrete values has lower resolution entropy as it cannot capture the continuous and subtle difference of gene expressions. (B) The software with continuous output values has higher resolution entropy.
## Fold-change-based evaluation
![Fold-change-based evaluation](figures/fold-change-based-evaluation.png)<br>
By counting the # of transcripts that are differentially expressed, statistics such as precision, recall, and accuracy can be defined and calculated.
## K-value
![K-value](figures/K-value.png)<br>
K-value is the condition number of the exon-isoform binary matrix, which can be used to measure the complexity of exon-isoform structures for each gene.
## Evaluation with different gene features
![Evaluation with different gene features](figures/split.jpg)<br>
Explore the relationship between the above evaluation metrics and important features like K-value, expression level, isoform length, and number of exons by splitting the transcripts into groups on these features.
