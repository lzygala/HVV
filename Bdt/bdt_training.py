import numpy as np
import glob
import pandas as pd
import xgboost as xgb
import tqdm
import math
import os

#Train BDT on selected signal
def train_bdt(df, features, year = "", sig_name = "WWH_ucsd_C2V_4",max_depth=5,eta=0.3,mcw=1,verb=False):
    if year != "":
        df_training = df[(df.split_train) & (df.year==year) & ((df.sample_isSignal == 0) | (df.Sample_name == sig_name))].sample(frac=1.)
        df_testing = df[(~df.split_train) & (df.year==year) & ((df.sample_isSignal == 0) | (df.Sample_name == sig_name))].sample(frac=1.)
    else:
        df_training = df[(df.split_train) & ((df.sample_isSignal == 0) | (df.Sample_name == sig_name))].sample(frac=1.)
        df_testing = df[(~df.split_train) & ((df.sample_isSignal == 0) | (df.Sample_name == sig_name))].sample(frac=1.)

    split_test_integral = df_testing.Xsec_genWeight.sum()
    split_train_integral = df_training.Xsec_genWeight.sum()
    
    xgbd_testing = xgb.DMatrix(
        df_testing[features.keys()], 
        label=df_testing.sample_isSignal,
        weight=np.abs(df_testing.split_weight),
        nthread=-1
    )
    xgbd_training = xgb.DMatrix(
        df_training[features.keys()], 
        label=df_training.sample_isSignal, 
        weight=np.abs(df_training.split_weight),
        nthread=-1
    )
    
    X_train, y_train = df_training[features].copy(), df_training.sample_isSignal.values
    w_train = (np.abs(df_training.split_weight)).values
    X_test, y_test = df_testing[features].copy(), df_testing.sample_isSignal.values
    w_test = (np.abs(df_testing.split_weight)).values
    
    X_train.rename(columns=features, inplace=True)
    X_test.rename(columns=features, inplace=True)

    #display(df_training_l)

    #evallist = [(xgbd_training, "train"), (xgbd_testing, "eval")]

    nRounds = 500
    params = {}
    params["objective"] = "binary:logistic"
    params["max_depth"] = max_depth
    params["eval_metric"] = "auc"
    params["verbosity"] = 1
    params["nthread"] = -1
    params["learning_rate"] = eta
    #params["gamma"] = gamma
    #params["max_leaves"] = mln
    #params["alpha"] = alpha
    #params["lambda"] = lda
    params["min_child_weight"] = mcw


    sumWeights_training_signal = np.abs(xgbd_training.get_weight()[xgbd_training.get_label() == 1]).sum()
    sumWeights_training_background = np.abs(xgbd_training.get_weight()[xgbd_training.get_label() == 0]).sum()
    params["scale_pos_weight"] = sumWeights_training_background/sumWeights_training_signal

    xgb_reg = xgb.XGBClassifier(**params)
    xgb_reg.fit(X_train, y_train, eval_set=[(X_test, y_test)], sample_weight = w_train, sample_weight_eval_set = [w_test], verbose=False, early_stopping_rounds=10)
    
    bestScore = xgb_reg.best_score
    print("Best Iteration: {}".format(bestScore))
    

    #print("Training Complete")
    return xgb_reg
    
#Train BDT on WWH+WZH
def train_bdt_mixed(df, features, year = "", signal_type = "C2V_4", max_depth=5,eta =0.3,mcw=1,verb=False):
    if year != "":
        df_training = df[(df.split_train) & (df.year==year) & ((df.sample_isSignal == 0) | (df.sig_type == signal_type))].sample(frac=1.)
        df_testing = df[(~df.split_train) & (df.year==year) & ((df.sample_isSignal == 0) | (df.sig_type == signal_type))].sample(frac=1.)
    else:
        df_training = df[(df.split_train) & ((df.sample_isSignal == 0) | (df.sig_type == signal_type))].sample(frac=1.)
        df_testing = df[(~df.split_train) & ((df.sample_isSignal == 0) | (df.sig_type == signal_type))].sample(frac=1.)

    if(len(df_training.index) == 0 or len(df_testing.index) == 0):
        print("0 occurred: ",year)
        return
    
    X_train, y_train = df_training[features].copy(), df_training.sample_isSignal.values
    w_train = (np.abs(df_training.split_weight)).values
    X_test, y_test = df_testing[features].copy(), df_testing.sample_isSignal.values
    w_test = (np.abs(df_testing.split_weight)).values
    
    X_train.rename(columns=features, inplace=True)
    X_test.rename(columns=features, inplace=True)

    xgbd_testing = xgb.DMatrix(
        df_testing[features.keys()], 
        label=df_testing.sample_isSignal,
        weight=np.abs(df_testing.split_weight),
        nthread=-1
    )
    xgbd_training = xgb.DMatrix(
        df_training[features.keys()], 
        label=df_training.sample_isSignal, 
        weight=np.abs(df_training.split_weight),
        nthread=-1
    )
    #evallist = [(xgbd_training, "train"), (xgbd_testing, "eval")]

    nRounds = 500
    params = {}
    params["objective"] = "binary:logistic"
    params["max_depth"] = max_depth
    params["eval_metric"] = "auc"
    params["verbosity"] = 1
    params["nthread"] = -1
    params["learning_rate"] = eta
    params["min_child_weight"] = mcw
    params["early_stopping_rounds"] =  10
    params["n_estimators"] =  500

    sumWeights_training_signal = np.abs(xgbd_training.get_weight()[xgbd_training.get_label() == 1]).sum()
    sumWeights_training_background = np.abs(xgbd_training.get_weight()[xgbd_training.get_label() == 0]).sum()
    params["scale_pos_weight"] = sumWeights_training_background/sumWeights_training_signal

    xgb_reg = xgb.XGBClassifier(**params)
    xgb_reg.fit(X_train, y_train, eval_set=[(X_test, y_test)], sample_weight = w_train, sample_weight_eval_set = [w_test], verbose=False, early_stopping_rounds=10)
    #xgb_reg.fit(X_train, y_train, eval_set=[(X_test, y_test)], sample_weight = w_train,  verbose=True)
    
    bestScore = xgb_reg.best_score
    print("Best Iteration: {}".format(bestScore))
    #print(xgb_reg.get_params())

    #print("Training Complete")
    return xgb_reg