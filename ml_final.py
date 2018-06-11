import pandas as pd
import scipy
import numpy as np

from sklearn.model_selection import train_test_split
from skmultilearn.problem_transform import BinaryRelevance, LabelPowerset, ClassifierChain
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score

from sklearn import metrics
from sklearn.externals import joblib

data = pd.read_csv("new_data.csv")
y = data[['antiviral', 'antibacterial', 'antifungal']]
to_drop = ['ID','Sequence','antiviral', 'antibacterial', 'antifungal']
X = data.drop(to_drop,axis=1)

X_train, X_test, y_train, y_test = train_test_split(X, y,test_size=0.1,random_state=55)

clf= BinaryRelevance(GradientBoostingClassifier(random_state=55))
clf.fit(X_train, y_train)
y_pred =clf.predict(X_test).tocoo()

print y_pred

print accuracy_score(y_pred, y_test)

# Start pickling
joblib.dump(clf, 'ml_model.pkl')
