import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import boruta_py

X = pd.read_csv('otutab.csv', index_col=0).values
y = pd.read_csv('mapfile_test.csv', index_col=0).values

# Converter
metahead = pd.read_csv('mapfile_test.csv', index_col=0)
colname = 'temp'
metahead['temp'].T.to_dict().values()



y = y.reshape((len(y),))

# define random forest classifier, with utilising all cores and
# sampling in proportion to y labels
rf = RandomForestClassifier(n_jobs=-1, class_weight='auto', max_depth=5)

# define Boruta feature selection method
feat_selector = boruta_py.BorutaPy(rf, n_estimators='auto', verbose=2)

# find all relevant features
feat_selector.fit(X, y)

# Get index of selected feature
select_feat = np.where( feat_selector.support_ == True )
header = pd.read_csv('otutab.csv', index_col=0).columns.values
signif_feature = [list(header)[y] for y in select_feat]



# check selected features
feat_selector.support_

# check ranking of features
feat_selector.ranking_

# call transform() on X to filter it down to selected features
X_filtered = feat_selector.transform(X)


