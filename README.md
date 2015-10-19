# EN-NegBin-Compare
A simple Python script to compare the output generated via [NegBinSig-Test](https://github.com/alifar76/NegBinSig-Test) and [MicrobeNets](https://github.com/alifar76/MicrobeNets), to see if there are any common OTUs/predictors between the two analysis approaches.

Background
------

Using a three-model approach (based on Poisson, Negative Binomial and Zero-inflated negative binomial regression models) to find significant OTUs between two treatment groups is one way of analysing the dataset. In this approach, the treatment groups are explanatory variables and OTU counts are response variables.

Another approach is to use logistic regression modelling based on penalized regression (as implemented via Lasso and Elastic Nets). In these approaches, OTU counts become explanatory variables and treatment groups become response variable. 

If one, from a biological standpoint, considers statistical models as nothing more than screening tools, one might be interested in knowing whether there is any biological convergance between various models. In such a case, it will be useful to see if we find a handful of predictors as significant, regardless of the test used. This is what this simple script does from the output produced by the two analysis scripts.

Running the script
------

There is 1 script in the folder src. It is called ```comparer.py```. 

To run the script, simply type:

```python comparer.py```

Please ensure that the two files are in the same folder as the script. Otherwise, the script will crash.

Output Explained
------

The output of the script contains information about OTUs that are common between the two analysis approaches. Currently, there are 5 columns in the output file prodcued. The columns and their descriptions of the output file are as follows:

1) **OTU_ID**: Indicates the OTU ID.

2) **Mean Difference**: Indicates the mean difference between treatment groups for the OTU. 

3) **NegBin Group**: Indicates which treatment group is the OTU significant in using the three model based regression modelling.

4) **EN Group**: Indicates which treatment group is the OTU significant in using logistic regression model based on elastic nets.

5) **Taxonomy**: Indicates the taxonomy of the OTU.
