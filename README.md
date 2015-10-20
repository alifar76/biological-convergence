# biological-convergence
This is a work in progress.

Background
------

Using a three-model approach (based on Poisson, Negative Binomial and Zero-inflated negative binomial regression models) to find significant OTUs between two treatment groups is one way of analysing the dataset. In this approach, the treatment groups are explanatory variables and OTU counts are response variables. It is implemented in [NegBinSig-Test](https://github.com/alifar76/NegBinSig-Test).

Another approach is to use logistic regression modelling based on penalized regression (as implemented via Lasso and Elastic Nets). In these approaches, OTU counts become explanatory variables and treatment groups become response variable as implemented in [MicrobeNets](https://github.com/alifar76/MicrobeNets).

If one, from a biological standpoint, considers statistical models as nothing more than screening tools, one might be interested in knowing whether there is any biological convergance between various models. In such a case, it will be useful to see if we find a handful of predictors as significant, regardless of the test used. This is what this pipeline attempts to do from the output produced by the two analysis scripts described earlier.

Required R packages
------

- [pscl](http://cran.r-project.org/web/packages/pscl/index.html)
- [MASS](http://cran.r-project.org/web/packages/MASS/index.html)
- [foreach](http://cran.r-project.org/web/packages/foreach/index.html)
- [doMC](http://cran.r-project.org/web/packages/doMC/index.html)
- [glmnet](http://cran.r-project.org/web/packages/glmnet/index.html)

Running the script
------

For help regarding the input parameters to run the pipline, which is present in ```src/``` folder, please type the following:

```python bio_convergence.py -h```

Input
------

The pipeline accepts an OTU table file, generated via [QIIME 1.8.0 (stable public release)](http://qiime.org/), (in tab-delimited format) as input and a standard mapping/metadata file compatible with QIIME. 

The publicly available dataset of the [lean-obese study](http://www.ncbi.nlm.nih.gov/pubmed/19043404), as obtained via [qiita](http://qiita.ucsd.edu), is provided with the pipeline as an example.


Output
------

The main output of the script is called **convergence_trt1_trt2.txt**, where **trt1** and **trt2** are the levels of the metadata variable used for prediction.
