# biological-convergence
This is a work in progress.

Background
------

Using a three-model approach (based on Poisson, Negative Binomial and Zero-inflated negative binomial regression models) to find significant OTUs between two treatment groups is one way of analysing the dataset. In this approach, the treatment groups are explanatory variables and OTU counts are response variables. It is implemented in [NegBinSig-Test](https://github.com/alifar76/NegBinSig-Test).

Another approach is to use logistic regression modelling based on penalized regression (as implemented via Lasso and Elastic Nets). In these approaches, OTU counts become explanatory variables and treatment groups become response variable as implemented in [MicrobeNets](https://github.com/alifar76/MicrobeNets).

If one, from a biological standpoint, considers statistical models as nothing more than screening tools, one might be interested in knowing whether there is any biological convergance between various models. In such a case, it will be useful to see if we find a handful of predictors as significant, regardless of the test used. This is what this pipeline attempts to do from the output produced by the two analysis scripts.

Running the script
------

For help regarding the input parameters to run the pipline, which is present in src/ folder, please type the following:

```python bio_convergence.py -h```

Input
------

The pipeline requires an QIIME compatible OTU table and mapping file as its input. The publicly available dataset of the lean/obese study, as obtained via qiita, is provided with the pipeline.


Output
------

The main output of the script is called **convergence_trt1_trt2.txt**, where **trt1** and **trt2** levels of the metadata variable used for prediction.
