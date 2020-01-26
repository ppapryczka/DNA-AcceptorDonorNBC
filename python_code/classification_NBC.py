import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split, cross_val_predict
from sklearn.naive_bayes import GaussianNB
import os

from python_code.get_acceptors_and_donors import (
    get_acceptors_and_donors,
    DEFAULT_DATA_FILE,
    dna_data_read,
)
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import RFE
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.utils import resample


def get_data_to_classifier(true_seq, false_seq, a: int, b: int) -> pd.DataFrame:
    file_path = f"{a}_{b}.csv"

    if os.path.isfile(file_path):
        columns = ["class"]
        for x in range(a + b):
            columns.append(f"{x}_A")
            columns.append(f"{x}_C")
            columns.append(f"{x}_G")
            columns.append(f"{x}_T")

        data_df = pd.DataFrame(columns=columns)

        for seq in true_seq:
            row = [1] + sum(list(map(nucleotide_to_table, list(seq))), [])
            data_df = data_df.append(pd.DataFrame([row], columns=columns))

        for seq in false_seq:
            row = [0] + sum(list(map(nucleotide_to_table, list(seq))), [])
            data_df = data_df.append(pd.DataFrame([row], columns=columns))

        data_df = data_df.reset_index(drop=True)

        data_df.to_csv(file_path)
        return data_df
    else:
        data_df = pd.read_csv(file_path)
        data_df = data_df.drop("Unnamed: 0", axis=1)

        return data_df


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


def simple_classify():
    df_dna = pd.read_csv("30_30.csv")
    df_dna = df_dna.drop("Unnamed: 0", axis=1)
    model = RandomForestClassifier()

    x_data = df_dna.drop("class", axis=1)
    """
    selected_columns = [14, 26, 29, 41, 77, 85, 89, 108, 109, 112, 118, 128, 132, 138, 143, 147, 151, 155, 157, 159, 163, 167, 171, 177, 179, 183, 187, 191, 195, 199, 203, 207, 211, 215, 219, 223, 227, 231, 235, 239]
    col = []
    for x in selected_columns:
        col.append(x_data.columns[x])
    x_data = x_data.loc[:, col]
    """
    model = ExtraTreesClassifier()
    model.fit(x_data, df_dna.loc[:, "class"])
    print(model.feature_importances_)

    """
    prediction = cross_val_predict(model, x_data , df_dna.loc[:, "class"], cv=10)

    mat = confusion_matrix(prediction, df_dna.loc[:, "class"])

    names = np.unique(prediction)
    sns.heatmap(
        mat,
        square=True,
        annot=True,
        fmt="d",
        cbar=False,
        xticklabels=names,
        yticklabels=names,
    )
    plt.xlabel("Truth")
    plt.ylabel("Predicted")
    plt.show()
    """


def classify():
    """
    with open(DEFAULT_DATA_FILE) as f:
        sequences = dna_data_split(f)
    true_d, false_d, true_a, false_a = get_acceptors_and_donors(DEFAULT_A, DEFAULT_B, sequences)

    df = get_data_to_classifier(true_d, false_d, 10, 10)
    """

    df = pd.read_csv("30_30.csv")
    df = df.drop("Unnamed: 0", axis=1)

    model = RandomForestClassifier()
    rfe = RFE(model, 40)

    rfe = rfe.fit(df.drop("class", axis=1), df.loc[:, "class"])

    print(rfe.support_)
    print(rfe.ranking_)


if __name__ == "__main__":
    """
    a = 20
    b = 20
    with open(DEFAULT_DATA_FILE) as f:
        sequences = dna_data_read(f)
    true_d, false_d, true_a, false_a = get_acceptors_and_donors(a, b, sequences)

    df_dna = get_data_to_classifier(true_d, false_d, a, b)
    print(df_dna)
    """
    features = [
        128,
        22,
        15,
        91,
        77,
        46,
        36,
        136,
        101,
        16,
        113,
        31,
        97,
        32,
        1,
        73,
        70,
        14,
        158,
        25,
        106,
        45,
        110,
        35,
        86,
        10,
        1,
        51,
        63,
        1,
        126,
        103,
        82,
        21,
        176,
        23,
        124,
        26,
        79,
        85,
        66,
        1,
        90,
        109,
        122,
        44,
        69,
        30,
        61,
        40,
        60,
        59,
        55,
        39,
        65,
        100,
        76,
        42,
        133,
        54,
        107,
        11,
        12,
        56,
        121,
        3,
        74,
        72,
        132,
        8,
        123,
        71,
        27,
        24,
        83,
        104,
        102,
        1,
        156,
        68,
        84,
        6,
        152,
        4,
        108,
        1,
        17,
        89,
        99,
        1,
        134,
        48,
        117,
        9,
        185,
        29,
        49,
        38,
        94,
        43,
        19,
        7,
        159,
        50,
        2,
        13,
        184,
        53,
        1,
        1,
        165,
        34,
        1,
        138,
        191,
        28,
        41,
        188,
        1,
        37,
        200,
        199,
        198,
        197,
        201,
        194,
        196,
        195,
        1,
        157,
        146,
        93,
        1,
        189,
        130,
        129,
        162,
        193,
        1,
        139,
        160,
        187,
        175,
        1,
        20,
        140,
        151,
        1,
        88,
        58,
        171,
        1,
        87,
        18,
        135,
        1,
        105,
        1,
        141,
        1,
        75,
        78,
        179,
        1,
        169,
        81,
        161,
        1,
        116,
        96,
        181,
        1,
        80,
        33,
        167,
        5,
        142,
        1,
        143,
        1,
        144,
        67,
        168,
        1,
        112,
        57,
        173,
        1,
        120,
        64,
        164,
        1,
        118,
        47,
        178,
        1,
        114,
        62,
        186,
        1,
        115,
        92,
        163,
        1,
        52,
        145,
        180,
        1,
        155,
        111,
        174,
        1,
        95,
        154,
        177,
        1,
        131,
        149,
        192,
        1,
        137,
        172,
        183,
        1,
        98,
        119,
        182,
        1,
        127,
        150,
        190,
        1,
        147,
        170,
        153,
        1,
        125,
        148,
        166,
        1,
    ]
    res = []
    for id, f in enumerate(features):
        if f == 1:
            res.append(id)
    print(res)

    simple_classify()

"""
plt.figure(figsize=(30, 20))

cor = df.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.savefig("result.png", dpi=500)

"""
"""
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
"""
""" 
names = np.unique(pred)
sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
            xticklabels=names, yticklabels=names)
plt.xlabel('Truth')
plt.ylabel('Predicted')
plt.show()
"""

"""
pred = model.predict(X_test)

# Plot Confusion Matrix
mat = confusion_matrix(pred, y_test)
names = np.unique(pred)
sns.heatmap(mat, square=True, annot=True, fmt='d', cbar=False,
            xticklabels=names, yticklabels=names)
plt.xlabel('Truth')
plt.ylabel('Predicted')
plt.show()
"""

"""
x = df.drop("class", 1)
y = df["class"]


plt.figure()
cor = df.corr()
sns.heatmap(cor, annot=True, cmap=plt.cm.Reds)
plt.show()
 """

"""
        '''
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
        '''
"""
"""

    df_majority = df[df.loc[:, "class"] == 0]
    df_minority = df[df.loc[:, "class"] == 1]

    # print(df_majority)
    # print(df_minority)

    df_majority_down_sampled = resample(
        df_majority, replace=False, n_samples=len(df_minority), random_state=123
    )
    # print(df_majority_downsampled)

    features = pd.concat([df_minority, df_majority_down_sampled])
    features = features.reset_index(drop=True)

    features_class = features.loc[:, "class"]

    features = features.drop("class", axis=1)
    features = features.iloc[
        :, [3, 4, 8, 11, 15, 19, 20, 24, 27, 29, 32, 38, 48, 52, 58, 63, 67, 71, 75, 79]
    ]
    print(features)

    model = RandomForestClassifier()

    model.fit(features, features_class)

    # pred = cross_val_predict(model, features, features_class, cv=20)

    pred = model.predict(
        df.drop("class", axis=1).iloc[
            :,
            [
                3,
                4,
                8,
                11,
                15,
                19,
                20,
                24,
                27,
                29,
                32,
                38,
                48,
                52,
                58,
                63,
                67,
                71,
                75,
                79,
            ],
        ]
    )

    mat = confusion_matrix(pred, df.loc[:, "class"])

    names = np.unique(pred)
    sns.heatmap(
        mat,
        square=True,
        annot=True,
        fmt="d",
        cbar=False,
        xticklabels=names,
        yticklabels=names,
    )
    plt.xlabel("Truth")
    plt.ylabel("Predicted")
    plt.show()

"""
