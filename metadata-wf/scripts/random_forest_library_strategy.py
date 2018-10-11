"""Use RandomForest to classify library strategy."""
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split, learning_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

features = pd.read_parquet(snakemake.input['features'])
Y = pd.read_parquet(snakemake.input['labels']).library_strategy.reindex(features.index)


# Set classes with fewer than 50 samples to OTHER
Ycnts = Y.value_counts()
drop_label = Ycnts.index[Ycnts < 50].tolist()

for label in drop_label:
    Y[Y == label] = 'OTHER'
    
# Split out OTHER
Y_OTHER = Y[Y == 'OTHER'].copy()
features_OTHER = features.reindex(Y_OTHER.index).dropna()

Y = Y[Y != 'OTHER'].copy()
features = features.reindex(Y.index)

# Encode labels and get starting matrices
encoder = LabelEncoder()
Y_enc = encoder.fit_transform(Y)
X = features.values

# Split into trianing and test data keeping class proportions similar (stratify)
X_train, X_test, Y_train, Y_test = train_test_split(X, Y_enc, stratify=Y_enc, random_state=42)

# Initialize classifier
clf =  RandomForestClassifier(n_estimators=200, oob_score=True)

# Now I am trying to figure out how training many samples are needed to 
# have good performance. To do this I use Learning Curve.
train_sizes, train_scores, test_scores = learning_curve(
    clf, 
    X_train, 
    Y_train, 
    cv=5, 
    train_sizes=np.linspace(0.001, 1, 20), 
    n_jobs=-1, 
)

train_mean = np.mean(train_scores, axis=1)
train_std = np.std(train_scores, axis=1)

test_mean = np.mean(test_scores, axis=1)
test_std = np.std(test_scores, axis=1)

plt.plot(train_sizes, train_mean, color='blue', marker='o', markersize=5, label='training accuracy')
plt.fill_between(train_sizes, train_mean + train_std, train_mean - train_std, alpha=0.15, color='blue')

plt.plot(train_sizes, test_mean, color='green', linestyle='--', marker='s', markersize=5, label='validation accuracy')
plt.fill_between(train_sizes, test_mean + test_std, test_mean - test_std, alpha=0.15, color='green')

plt.xlabel('Training Sample Size')
plt.ylabel('Score')
plt.legend(loc='lower right')
fig = plt.gcf()
fig.savefig(snakemake.output['learning_curve'])

# Get a global view of model performance
clf.fit(X_train, Y_train)
oob = clf.oob_score_
Y_pred = clf.predict(X_test)
with open(snakemake.output['metrics'], 'w') as fh:
    fh.write(f'OOB: {oob}\n')
    fh.write(classification_report(Y_test, Y_pred, target_names=encoder.classes_))

importances = pd.Series(clf.feature_importances_, index=features.columns).sort_values(ascending=False)
importances.to_csv(snakemake.output['feature_importances'], sep='\t')

# Make predictions on all data and OTHER using a 20-fold cross validation approach. 
# Here I am using only 20% of the data for training given the learning curve results above.
n = 20
Y_res = np.empty((Y_enc.shape[0], n), dtype=object)
Y_res_other = np.empty((Y_OTHER.shape[0], n), dtype=object)
for i in range(n):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y_enc, stratify=Y_enc, test_size=.8)
    clf.fit(X_train, Y_train)
    Y_pred = clf.predict(X)
    Y_res[:, i] = encoder.inverse_transform(Y_pred)
    
    Y_pred2 = clf.predict(features_OTHER.values)
    Y_res_other[:, i] = encoder.inverse_transform(Y_pred2)

YPred = pd.DataFrame(Y_res, index=features.index)
YPred.columns = [f'RF{i}' for i in YPred.columns]
YPred.to_parquet(snakemake.output['pred'])

YPredOther = pd.DataFrame(Y_res_other, index=features_OTHER.index)
YPredOther.columns = [f'RF{i}' for i in YPredOther.columns]
YPredOther.to_parquet(snakemake.output['pred_other'])