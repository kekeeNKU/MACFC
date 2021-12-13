from pandas import *
from numpy import *


def auto_norm(data):
    # data:（sample,feature）
    mins = data.min(0)
    maxs = data.max(0)
    ranges = maxs - mins
    row = data.shape[0]
    norm_data = data - tile(mins, (row, 1))
    norm_data = norm_data / tile(ranges, (row, 1))
    return norm_data


def load_arcene():
    data = read_csv('./UCI/arcene_train.csv', header=None)
    label = read_csv('./UCI/arcene_train_labels.csv', header=None)
    f, c = data.values, label.values
    values = auto_norm(unique(hstack((f, c)), axis=0))
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_credit():
    data = read_csv('./UCI/credit.csv', header=None)
    values = unique(data.values, axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_madelon():
    data = read_csv('./UCI/madelon.csv', header=None)
    values = unique(transpose(data.values), axis=0)
    f = auto_norm(values[:, 1:])  # data normalization for constant value
    c = values[:, 0]
    return transpose(f), c


def load_musk():
    data = loadtxt('./UCI/clean1.txt', delimiter=',', dtype=str)
    data = data[:, 2:].astype(float)
    values = unique(data, axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_ult():
    data = read_csv('./UCI/Meter A.csv', header=None)
    values = unique(data.values, axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_parkinson():
    data = read_csv('./UCI/ParkinsonDisease.csv')
    values = unique(data.values[1:, :], axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_sonar():
    data = read_csv('./UCI/sonar.csv', header=None)
    values = unique(data.values, axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c


def load_TCGA(file):
    # data = read_csv('./TCGA/genomicMatrix-KIRC.xls', dtype=str, delimiter='\t', header=None)
    data = read_csv(file, dtype=str, delimiter='\t', header=None)
    label = [data.iloc[0, :][0]] + [i.strip().split('-')[-1] for i in data.iloc[0, :][1:]]
    dl = [label.index(i) for i in label if i not in ['sample', '01', '11']]
    data.iloc[0, :] = array(label)
    data = data.drop(columns=dl)
    tfdata = transpose(data.values)
    tfvalues = tfdata[1:, :].astype(float)
    FName = tfdata[0, 1:]
    c = tfvalues[:, 0]
    f = tfvalues[:, 1:]
    return FName, transpose(f), c


def load_mc_data(file):
    # data = read_csv('./multi_class/DLBCL.csv', header=None)
    data = read_csv(file, header=None)
    values = unique(data.values, axis=0)
    f = auto_norm(values[:, :-1])  # data normalization for constant value
    c = values[:, -1]
    # return feature*sample matrix and class labels
    return transpose(f), c

