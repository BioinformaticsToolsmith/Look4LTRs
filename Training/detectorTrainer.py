'''
Trains Look4LTRs detector module
'''

import argparse
import os


parser = argparse.ArgumentParser(description='Train Look4LTRs detector')
parser.add_argument('--feature_file', required=True, nargs='+', help='Feature file; should contain label and features for each stretch pair')
parser.add_argument('--output', required=True, help='Output config file for trained model')


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Get arguments and validate 
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
args = parser.parse_args()
feature_file_list = args.feature_file
output_file = args.output

for feature_file in feature_file_list:
    if not os.path.exists(feature_file):
        print("Feature file {} does not exist".format(feature_file))
        exit(1)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Imports (put after arg parse to avoid unnecessary imports if args are wrong)
#@#@#@#@#@#@#@#@#@#@#@#@#@#@


import warnings
warnings.filterwarnings("ignore")

import gc

import statistics

from sklearn.model_selection import train_test_split
import numpy as np

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import SGDClassifier

from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import cross_validate
from scipy.stats import uniform
from sklearn.metrics import make_scorer
from sklearn.metrics import recall_score, f1_score, precision_score, accuracy_score

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Functions
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

def specificity_score(y, y_pred):
    return recall_score(y, y_pred, pos_label=0, zero_division=0)


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Loading Data
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Loading feature files in...")
# loading in all feature tables and concatenating them
feature_table = []
for feature_file in feature_file_list:
    feature_table.append(np.loadtxt(feature_file, delimiter=','))

feature_table = np.concatenate(feature_table, axis=0)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Splitting into training, validation and testing (70-20-10)
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Splitting Data...")
# Splitting into training and testing with stratification
X = feature_table[:, 1:]
y = feature_table[:, 0].astype(np.int8)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.1, stratify=y)

# Splitting training into training and validation
train_val_ratio = 0.2/0.9
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=train_val_ratio, stratify=y_train)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Removing old data for memory
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

del feature_table

gc.collect()


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Creating Pipelines
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

# list of column indexes to standardize (last column is binary so dont)
standardize_columns = list(range(X_train.shape[1]-1))

ct = ColumnTransformer( [('standardize', StandardScaler(), standardize_columns)], remainder='passthrough' )

pipeline = Pipeline([
    ('preprocessing', ct),
    ('sgd', SGDClassifier())
])

# Parameters to randomize in RandomizedSearchCV
sgd_params = [{'sgd__loss': ['hinge', 'log_loss', 'modified_huber', 'squared_hinge', 'perceptron'],
          'sgd__alpha': uniform(), 
          'sgd__penalty': ['l2', 'l1']
             }, 
             {'sgd__l1_ratio': uniform(), 'sgd__penalty': ['elasticnet']
             },
             {'sgd__loss': ['hinge', 'log_loss', 'modified_huber', 'squared_hinge', 'perceptron'],
          'sgd__penalty': ['l2', 'l1'], 
          'sgd__learning_rate': ['constant'],
          'sgd__eta0': uniform()
             },
             {'sgd__loss': ['hinge', 'log_loss', 'modified_huber', 'squared_hinge', 'perceptron'],
          'sgd__penalty': ['l2', 'l1'], 
          'sgd__learning_rate': ['invscaling'],
          'sgd__eta0': uniform(),
          'sgd__power_t': uniform()
             },
             {'sgd__loss': ['hinge', 'log_loss', 'modified_huber', 'squared_hinge', 'perceptron'],
          'sgd__penalty': ['l2', 'l1'], 
          'sgd__learning_rate': ['adaptive'],
          'sgd__early_stopping': [True],
          'sgd__eta0': uniform()
             },
            ]

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Performing RandomizedSearchCV
#@#@#@#@#@#@#@#@#@#@#@#@#@#@


print("Performing random search cross-validation...")
scorers = {
    'recall': make_scorer(recall_score, zero_division=0),
    'precision': make_scorer(precision_score, zero_division=0),
    'f1': make_scorer(f1_score, zero_division=0),
    'specificity': make_scorer(specificity_score)
}

min_class_size = np.min(np.bincount(y_train))

n_splits = min(10, min_class_size)
sgd_search = RandomizedSearchCV(pipeline, sgd_params, n_iter=10, scoring=scorers, refit='f1', cv=n_splits)

sgd_search.fit(X_train, y_train)

# Best parameters for the model
best_params = sgd_search.best_params_
best_params = {k.replace('sgd__', ''): v for k, v in best_params.items()}


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Using Best Parameters to Train Model 1000 times
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

print("Retraining model to get best weights...")
best_f1 = -1
best_model = None

for i in range(10):
    pipeline = Pipeline([
    ('preprocessing', ct),
    ('sgd', SGDClassifier(**best_params))
    ])

    pipeline.fit(X_train, y_train)
    y_pred = pipeline.predict(X_val)
    f1 = f1_score(y_val, y_pred)
    if f1 > best_f1:
        best_f1 = f1
        best_model = pipeline

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Evaluate Model on Test Data
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

y_pred = best_model.predict(X_test)
recall = recall_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)
f1 = f1_score(y_test, y_pred)
specificity = specificity_score(y_test, y_pred)
accuracy = accuracy_score(y_test, y_pred)

print("Best model used the following parameters: {}".format(best_params))
print("Recall: {}".format(recall))
print("Precision: {}".format(precision))
print("F1: {}".format(f1))
print("Specificity: {}".format(specificity))
print("Accuracy: {}".format(accuracy))


#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Standardizing Features according to entire dataset
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

ct.fit(X)

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Getting mean and standard deviation of features
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

standard_scaler = ct.named_transformers_['standardize']


mean = standard_scaler.mean_
std = standard_scaler.scale_

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Getting weights of model and intercept
#@#@#@#@#@#@#@#@#@#@#@#@#@#@

weights = best_model.named_steps['sgd'].coef_
intercept = best_model.named_steps['sgd'].intercept_

#@#@#@#@#@#@#@#@#@#@#@#@#@#@
# Saving to file
#@#@#@#@#@#@#@#@#@#@#@#@#@#@
with open(output_file, 'w') as f:
    # Writing mean values
    f.write("mean: " + " ".join(map(str, mean)) + "\n")
    # Writing standard deviation values
    f.write("std: " + " ".join(map(str, std)) + "\n")
    # Writing the intercept
    f.write("intercept: " + " ".join(map(str, intercept)) + "\n")
    # Writing weight vectors (assuming a single row of weights for simplicity)
    f.write("weights: " + " ".join(map(str, weights[0])) + "\n")
    
    
    
