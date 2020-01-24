import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.naive_bayes import GaussianNB, ComplementNB, MultinomialNB, BernoulliNB, CategoricalNB

from python_code.get_acceptors_and_donors import get_acceptors_and_donors, DEFAULT_DATA_FILE, DEFAULT_A, DEFAULT_B, dna_data_split
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample

import itertools


def get_data_to_classifier(true_seq, false_seq, a: int , b: int) -> pd.DataFrame:
    columns = ["class"]
    for x in range(a+b):
        columns.append(f"{x}_A")
        columns.append(f"{x}_C")
        columns.append(f"{x}_G")
        columns.append(f"{x}_T")
    df = pd.DataFrame(columns=columns)

    seq_row = []

    for id, seq in enumerate(true_seq):
        seq_row.append(1)
        for nucleotide in seq:
            seq_row += nucleotide_to_table(nucleotide)
        df = df.append(pd.DataFrame([seq_row], columns=columns))
        seq_row = []

    for id, seq in enumerate(false_seq):
        seq_row.append(0)
        for nucleotide in seq:
            seq_row += nucleotide_to_table(nucleotide)
        df = df.append(pd.DataFrame([seq_row], columns=columns))
        seq_row = []

    df = df.reset_index(drop=True)

    df.to_csv(f"{a}_{b}.csv")

    return df


def nucleotide_to_table(nucleotide: str):
    if nucleotide == "A":
        return [1, 0, 0, 0]
    elif nucleotide == "C":
        return [0, 1, 0, 0]
    elif nucleotide == "G":
        return [0, 0, 1, 0]
    elif nucleotide == "T":
        return [0, 0, 0, 1]
    else:
        return [1, 1, 1, 1]


def get_only_specific_values(positions, df):
    names = ["class"]
    for x in positions:
        names.append(f"{x}_A")
        names.append(f"{x}_C")
        names.append(f"{x}_T")
        names.append(f"{x}_G")
    return df[names]


def simple_classify():
    df = pd.read_csv("test.csv")
    df = df.drop("Unnamed: 0", axis=1)

    model = RandomForestClassifier()
    pred = cross_val_predict(model, df.drop("class", axis=1), df.loc[:, "class"], cv=20)

    mat = confusion_matrix(pred, df.loc[:, "class"])

    names = np.unique(pred)
    sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
                xticklabels=names, yticklabels=names)
    plt.xlabel('Truth')
    plt.ylabel('Predicted')
    plt.show()

def classify():
    '''
    with open(DEFAULT_DATA_FILE) as f:
        sequences = dna_data_split(f)
    true_d, false_d, true_a, false_a = get_acceptors_and_donors(DEFAULT_A, DEFAULT_B, sequences)

    df = get_data_to_classifier(true_d, false_d, 10, 10)
    '''

    df = pd.read_csv("10_10.csv")
    df = df.drop("Unnamed: 0", axis=1)

    #model = RandomForestClassifier()
    #rfe = RFE(model, 20)

    #rfe = rfe.fit(df.drop("class", axis=1), df.loc[:, "class"])

    #print(rfe.support_)
    #print(rfe.ranking_)

    df_majority = df[df.loc[: , "class"] == 0]
    df_minority = df[df.loc[:, "class"] == 1]

    #print(df_majority)
    #print(df_minority)

    df_majority_down_sampled = resample(df_majority,
                                       replace=False,
                                       n_samples=len(df_minority),
                                       random_state=123)
    #print(df_majority_downsampled)

    features = pd.concat([df_minority, df_majority_down_sampled])
    features = features.reset_index(drop=True)

    features_class = features.loc[:, "class"]

    features = features.drop("class", axis=1)
    features = features.iloc[:, [3, 4, 8, 11, 15, 19, 20, 24, 27, 29, 32, 38, 48, 52, 58, 63, 67, 71, 75, 79]]
    print(features)

    model = RandomForestClassifier()

    model.fit(features, features_class)

    #pred = cross_val_predict(model, features, features_class, cv=20)

    pred = model.predict(df.drop("class", axis=1).iloc[:, [3, 4, 8, 11, 15, 19, 20, 24, 27, 29, 32, 38, 48, 52, 58, 63, 67, 71, 75, 79]])

    mat = confusion_matrix(pred, df.loc[:,"class"])

    names = np.unique(pred)
    sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
                xticklabels=names, yticklabels=names)
    plt.xlabel('Truth')
    plt.ylabel('Predicted')
    plt.show()



if __name__ == "__main__":
    simple_classify()


'''
plt.figure(figsize=(30, 20))

cor = df.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.savefig("result.png", dpi=500)

'''
'''
model = GaussianNB()

#X_train, X_test, y_train, y_test = train_test_split(tmp_df.drop("class", axis=1), tmp_df.loc[:,"class"], train_size=0.7, random_state = 42)

# Train the model
#model.fit(X_train, y_train)

#pred = model.predict(X_test)

# Predict Output
pred = cross_val_predict(model, tmp_df.drop("class", axis=1), tmp_df.loc[:,"class"], cv=20)

# Plot Confusion Matrix
#mat = confusion_matrix(pred, y_test)
mat = confusion_matrix(pred, tmp_df.loc[:,"class"])
'''
''' 
names = np.unique(pred)
sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
            xticklabels=names, yticklabels=names)
plt.xlabel('Truth')
plt.ylabel('Predicted')
plt.show()
'''

'''
pred = model.predict(X_test)

# Plot Confusion Matrix
mat = confusion_matrix(pred, y_test)
names = np.unique(pred)
sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
            xticklabels=names, yticklabels=names)
plt.xlabel('Truth')
plt.ylabel('Predicted')
plt.show()
'''

'''
x = df.drop("class", 1)
y = df["class"]


plt.figure()
cor = df.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.show()
 '''
