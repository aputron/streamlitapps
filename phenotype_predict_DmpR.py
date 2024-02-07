#%%
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd
import tensorflow.keras.backend as K
from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten
from tensorflow.keras.models import Sequential
from tensorflow import keras as k
from keras_preprocessing.sequence import pad_sequences
from sklearn.preprocessing import LabelEncoder, OneHotEncoder


#%%
def mcc(label_real, label_predicted):
    """Custom Metric for Accuracy and Precision

    Parameters:
    ----------
    label_real: tensorflow.Tensor
        Correct labels
    label_predicted: tensorflow.Tensor
        Predicted labels

    Returns
    -------
    float
        Custom metric
    """
    true_pos = K.sum(K.round(K.clip(label_real * label_predicted, 0, 1)))
    true_neg = K.sum(K.round(K.clip((1 - label_real) * (1 - label_predicted), 0, 1)))
    false_pos = K.sum(K.round(K.clip((1 - label_real) * label_predicted, 0, 1)))
    false_neg = K.sum(K.round(K.clip(label_real * (1 - label_predicted), 0, 1)))
    number = true_pos * true_neg - false_pos * false_neg
    denominator = K.sqrt((true_pos + false_pos) * 
                        (true_pos + false_neg) * 
                        (true_neg + false_pos) * 
                        (true_neg + false_neg))
    return number / (denominator + K.epsilon())

def one_hot_encoding_aa(dataset):
    """One hot encoding of Amino Acid sequences

    Parameters:
    ----------
    dataset: List[str]
        List of Amino Acid sequences

    Returns
    -------
    numpy.ndarray
        One Hot Encoding of the amino acids (has dimensions dataset_size x seq_len x 21)
    """
    integer_encoder = LabelEncoder()
    one_hot_encoder = OneHotEncoder(categories='auto')
    amino_acids = "ARNDCQEGHILKMFPSTWYV*"
    input_features = []

    # fix the encoded categories
    ie = integer_encoder.fit_transform(list(amino_acids)) #.toarray().reshape(-1, 1)
    ie = np.array(ie).reshape(-1, 1)
    oe = one_hot_encoder.fit_transform(ie)

    for sequence in dataset:
        if type(sequence) == str:
            integer_encoded = integer_encoder.transform(list(sequence)) #.toarray().reshape(-1, 1)
            integer_encoded = np.array(integer_encoded).reshape(-1,1)
            one_hot_encoded = one_hot_encoder.transform(integer_encoded)
            input_features.append(one_hot_encoded.toarray())
            
    input_features = pad_sequences(input_features, padding="post")
    input_features = np.stack(input_features)
    
    return input_features

def one_hot_encoding(dataset, maxlen):
    """One hot encoding of DNA sequences

    Parameters:
    ----------
    dataset: List[str]
        List of DNA sequences
    maxlen: int
        Max length of padded sequences

    Returns
    -------
    numpy.ndarray
        One Hot Encoding of the DNA sequences (has dimensions dataset_size x maxlen x 4)
    """
    integer_encoder = LabelEncoder()
    one_hot_encoder = OneHotEncoder(categories='auto')
    input_features = []
    for sequence in dataset:
        if type(sequence) == str:
            integer_encoded = integer_encoder.fit_transform(list(sequence)) #.toarray().reshape(-1, 1)
            integer_encoded = np.array(integer_encoded).reshape(-1,1)
            one_hot_encoded = one_hot_encoder.fit_transform(integer_encoded)
            input_features.append(one_hot_encoded.toarray())

    input_features = pad_sequences(input_features, padding="post", maxlen=maxlen)
    input_features = np.stack(input_features)
    
    return input_features

candidates = []
description = []
for record in Bio.SeqIO.parse("Insilico/function_test_4-5_aa.fasta", "fasta"):
    candidat = str(record.seq)
    candidates.append(candidat)
    description.append(record.description)
candidate = one_hot_encoding(candidates, 1784)

def cnn_1d(train_features):
    model = Sequential()
    model.add(Conv1D(filters=49, kernel_size=3, input_shape=(train_features, 4)))
    model.add(MaxPooling1D(pool_size=4))
    model.add(Conv1D(filters=49, kernel_size=3, input_shape=(train_features, 4)))
    model.add(MaxPooling1D(pool_size=4))
    # model.add(Flatten())
    model.add(Conv1D(filters=49, kernel_size=3, input_shape=(train_features, 4)))
    model.add(MaxPooling1D(pool_size=4))
    #model.add(Conv1D(filters=400, kernel_size=3, input_shape=(train_features, 4)))
    #model.add(MaxPooling1D(pool_size=4))
    model.add(Flatten())
    model.add(Dense(16, activation='elu'))
    model.add(Dense(1, activation='sigmoid'))

    model.compile(loss='binary_crossentropy', optimizer='sgd', metrics=['binary_accuracy', mcc])

    # if you need to convert the model to a multi-class classification model
    #model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['categorical_accuracy'])

    model.summary()
    return model

dependencies = {
    'mcc': mcc
}

#signal_model = keras.models.load_model('CNN_aa/trained_models/mcc_BC_opt_aaprop_signal_v0.1.h5', compile=False)
#bkg_model = keras.models.load_model('CNN_aa/trained_models/mcc_BC_opt_aaprop_bkg_v0.2.h5', compile=False)

#signal_model = keras.models.load_model('CNN_aa/trained_models/mcc_BC_opt_aaprop_signal_v0.1.h5', compile=False)
background_model = keras.models.load_model('CNN_multilabel/trained_models/binary_classification_CNN1D_bkg_v2.3.h5', compile=False)
#signal_model.load_weights("CNN_multilabel/val_signal_best_weights_4.h5")
signal_model = cnn_1d(1784)
signal_model.load_weights("CNN_multilabel/opt_sig_best.h5")

# with h5py.File('CNN_multilabel/opt_sig_best.h5', 'r') as f:
#     # Iterate through the layers in the model
#     for layer_name in f.keys():
#         layer = f[layer_name]
#         print(layer.values)
#         for i in layer:
#             print(i)

signal_prediction = signal_model.predict(candidate)#.flatten().tolist()
bkg_prediction = background_model.predict(candidate)#.flatten().tolist()

THRESHOLD_SIGNAL = 0.52
THRESHOLD_BACKGROUND = 0.52

# pred_labels = []
# for i in signal_prediction:
#     if i[0] > i[1]:
#         pred = "low"
#     else:
#         pred = "high"
#     pred_labels.append(pred)
# pred_labels = []
# for i in np.round(signal_prediction):
#      pred_labels.append("high" if i > THRESHOLD else "low")

# def compare_predict(softmax_sig=False, softmax_bkg=False):
#     if softmax_sig:
#         sig_pred = [("low" if i[0] > i[1] else "high") for i in signal_prediction]
#     else:
#         sig_pred = [("low" if i < THRESHOLD_SIGNAL else "high") for i in np.round(signal_prediction)]

#     if softmax_bkg:
#         bkg_pred = [("low" if i[0] > i[1] else "high") for i in bkg_prediction]
#     else:
#         bkg_pred = [("low" if i < THRESHOLD_BACKGROUND else "high") for i in np.round(bkg_prediction)]

def multi_pheno_predict(signal_prediction, bkg_prediction, signal_first=False):
    bkg_predi = [i[1] for i in bkg_prediction]
    bkg_pred = pd.DataFrame(bkg_predi, columns=["Background"])#, index=description)
    sig_pred = pd.DataFrame(signal_prediction, columns=["Signal"])#, index=description)

    prediction = pd.concat([sig_pred, bkg_pred], axis=1)

    if signal_first:
        sort = prediction.sort_values(['Signal', 'Background'], ascending=[False, True])
    else:
        sort = prediction.sort_values(['Background', 'Signal'], ascending=[True, False])

    hslb, hshb, lshb, lslb = [], [], [], []

    for index, row in sort.iterrows():
        if round(row['Background']) == 0:
            if not "*" in description[index]:
                hslb.append(index)
        elif round(row['Background']) == 1:
            hshb.append(index)

    for index, row in sort[::-1].iterrows():
        if round(row['Background']) == 0:
            lslb.append(index)
        elif round(row['Background']) == 1:
            lshb.append(index)
                
    hslb_f = hslb[:101]
    hshb_f = hshb[:100]
    lshb_f = lshb[:100]
    lslb_f = lslb[:100]

    # for idx, row in prediction.iterrows():
    #     sig_val, bkg_val = row['Signal'], row['Background']
    #     if round(sig_val) == 0 and round(bkg_val) == 0:
    #         lslb.append(idx)
    #     elif round(sig_val) == 0 and round(bkg_val) == 1:
    #         lshb.append(idx)
    #     elif round(sig_val) == 1 and round(bkg_val) == 0:
    #         if not "*" in description[index]:
    #             hslb.append(idx)
    #     elif round(sig_val) == 1 and round(bkg_val) == 1:
    #         hshb.append(idx)

    # lslb.sort(key=lambda idx: (prediction.iloc[idx]['Background'], prediction.iloc[idx]['Signal']))
    # lshb.sort(key=lambda idx: (prediction.iloc[idx]['Background'], -prediction.iloc[idx]['Signal']))
    # hslb.sort(key=lambda idx: (-prediction.iloc[idx]['Background'], prediction.iloc[idx]['Signal']))
    # hshb.sort(key=lambda idx: (-prediction.iloc[idx]['Background'], -prediction.iloc[idx]['Signal']))
    # Rename lslb to lslb_f (also others)

prediction.iloc[hslb_f]

prediction.iloc[lshb]
prediction.iloc[lslb]


fasta_dna = []
fasta_aa = []
for i in hslb_f:
    signal = signal_prediction[i][0]
    bkg = bkg_predi[i]
    dna_read = SeqRecord(Seq(candidates[i]), id=description[i], description=f"{signal}  {bkg}")
    fasta_dna.append(dna_read)
Bio.SeqIO.write(fasta_dna, "high_candidates_dna.fasta", "fasta")


fasta_dna = []
fasta_aa = []
for i in range(len(hslb_f)):
    diff = dif(ref_aa, candidates_aa[i])
    description=[]
    for x in diff:
        des = f"{ref_aa[x]}{x + 1}{candidates_aa[i][x]}"
        description.append(des)
    description = "|".join(description)
    dna_read = SeqRecord(Seq(candidates_dna[i]), id = "DmpR_insilico_mutant", description=description)
    fasta_dna.append(dna_read)
    aa_read = SeqRecord(Seq(candidates_aa[i]), id = "DmpR_insilico_mutant", description=description)
    fasta_aa.append(aa_read)
Bio.SeqIO.write(fasta_dna, f"{filename}_dna.fasta", "fasta")
Bio.SeqIO.write(fasta_aa, f"{filename}_aa.fasta", "fasta")
    