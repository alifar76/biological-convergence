# biological-convergence
This is a work in progress.

Background
------

Using a three-model approach (based on Poisson, Negative Binomial and Zero-inflated negative binomial regression models) to find significant OTUs between two treatment groups is one way of analysing the dataset. In this approach, the treatment groups are explanatory variables and OTU counts are response variables. It is implemented in [NegBinSig-Test](https://github.com/alifar76/NegBinSig-Test).

Another approach is to use logistic regression modelling based on penalized regression (as implemented via Lasso and Elastic Nets). In these approaches, OTU counts become explanatory variables and treatment groups become response variable as implemented in [MicrobeNets](https://github.com/alifar76/MicrobeNets).

If one, from a biological standpoint, considers statistical models as nothing more than screening tools, one might be interested in knowing whether there is any biological convergance between various models. In such a case, it will be useful to see if we find a handful of predictors as significant, regardless of the test used. This is what this pipeline attempts to do from the output produced by the two analysis scripts.

Running the script
------

There is 1 script in the folder src. It is called ```comparer.py```. 

To run the script, simply type:

```python bio_convergence.py```

For help regarding the input parameters, please type:

```python bio_convergence.py -h``

Output Explained
------

The main output of the script is called:


5) **Taxonomy**: Indicates the taxonomy of the OTU.
