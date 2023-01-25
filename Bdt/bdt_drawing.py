
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def draw_bdt_score(df):
    fig, axes = plt.subplots(1,2,figsize=(20, 10))

    bkg_weights = df.loc[df["sample_isSignal"] == 0]["Xsec_genWeight"] 
    sig_weights = df.loc[df["sample_isSignal"] == 1]["Xsec_genWeight"]

    plotData_background = df.loc[df["sample_isSignal"] == False]["bdt"]
    plotData_signal = df.loc[df["sample_isSignal"] == True]["bdt"]

    axes[0].hist(
        plotData_background,
        bins=np.linspace(0, 1, 101),
        weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events",
        log=False
    )

    axes[0].hist(
        plotData_signal,
        bins=np.linspace(0, 1, 101),
        weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights):0.1f} Events",
        log=False
    )

    axes[0].set_xlabel("BDT Score", size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=16)
    
    axes[1].hist(
        plotData_background,
        bins=np.linspace(0, 1, 101),
        weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events",
        log=True
    )

    axes[1].hist(
        plotData_signal,
        bins=np.linspace(0, 1, 101),
        weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights):0.1f} Events",
        log=True
    )

    axes[1].set_xlabel("BDT Score", size=18)
    axes[1].set_ylabel("Log(Events)", size=18)
    axes[1].legend(fontsize=16)
    
def draw_bdt_score_overlay(dfs={}):
    fig, axes = plt.subplots(1,2,figsize=(20, 10))

    for year, df in dfs.items():
        bkg_weights = df.loc[df["sample_isSignal"] == 0]["Xsec_genWeight"] 
        sig_weights = df.loc[df["sample_isSignal"] == 1]["Xsec_genWeight"]

        plotData_background = df.loc[df["sample_isSignal"] == False]["bdt"]
        plotData_signal = df.loc[df["sample_isSignal"] == True]["bdt"]

        axes[0].hist(
            plotData_background,
            bins=np.linspace(0, 1, 101),
            weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
            histtype="step",
            label=year + f" Background : {sum(bkg_weights):0.1f} Events",
            log=False
        )

        axes[0].hist(
            plotData_signal,
            bins=np.linspace(0, 1, 101),
            weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
            histtype="step",
            #color="r",
            label=year + f"Signal : {sum(sig_weights):0.1f} Events",
            log=False
        )

        axes[0].set_xlabel("BDT Score", size=18)
        axes[0].set_ylabel("Events", size=18)
        axes[0].legend(fontsize=16)

        axes[1].hist(
            plotData_background,
            bins=np.linspace(0, 1, 101),
            weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
            histtype="step",
            label=year + f"Background : {sum(bkg_weights):0.1f} Events",
            log=True
        )

        axes[1].hist(
            plotData_signal,
            bins=np.linspace(0, 1, 101),
            weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
            histtype="step",
            #color="r",
            label=year + f"Signal : {sum(sig_weights):0.1f} Events",
            log=True
        )

        axes[1].set_xlabel("BDT Score", size=18)
        axes[1].set_ylabel("Log(Events)", size=18)
        axes[1].legend(fontsize=16)

def get_bdt_fom(df, cut):
    s = df.loc[(df["sample_isSignal"] == 1) & (df["bdt"] > cut)]["abs_Xsec_genWeight"].sum()
    b = df.loc[(df["sample_isSignal"] == 0) & (df["bdt"] > cut)]["abs_Xsec_genWeight"].sum()
    #print(s,b)
    #return (s/math.sqrt(s+b))
    #return (s/math.sqrt(b))
    return ( math.sqrt(  2*(  ((s+b)*math.log(1+(s/b))-s)  )  )  )
        
def draw_bscore_fom_overlay(dfs={},xmin=0.0,xmax=1.0,xbins=100):        
    fig, axes = plt.subplots(1,1,figsize=(10, 10))

    for year, df in dfs.items():
        
        fom = {}
        for i in range(1, xbins+1):
            cut = xmin + (((xmax-xmin) / xbins) * i)
            fom[cut]=get_bdt_fom(df,cut)
        axes.hist(
            fom.keys(), weights=fom.values(),bins=np.linspace(xmin, xmax, xbins+1),
            histtype="step",
            label=year,
            log=False
        )

        axes.set_xlabel("BDT Score Cut", size=18)
        axes.set_ylabel("S / Sqrt(S+B)", size=18)
        axes.legend(fontsize=16, loc="lower right")
    plt.show()
        
def draw_plots(df, var, bdt_cut, name, bins,fs=16):
    fig, axes = plt.subplots(1,2,figsize=(16, 8))

    plotData_background = df.loc[df["sample_isSignal"] == False][var]
    plotData_signal = df.loc[df["sample_isSignal"] == True][var]

    plotData_background_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)][var]
    plotData_signal_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)][var]

    bkg_weights = df.loc[(df["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights = df.loc[(df["sample_isSignal"] == True)]["Xsec_genWeight"]

    bkg_weights_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"] 
    sig_weights_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"]

    axes[0].hist(
        plotData_background,
        bins=bins,
        weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events : {len(bkg_weights)} Raw Events",
        log=False
    )

    axes[0].hist(
        plotData_signal,
        bins=bins,
        weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights):0.1f} Events : {len(sig_weights)} Raw Events",
        log=False
    )

    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=fs)

    axes[1].hist(
        plotData_background_cut,
        bins=bins,
        weights=df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights_cut):0.1f} Events : {len(bkg_weights_cut)} Raw Events",
        log=False
    )

    axes[1].hist(
        plotData_signal_cut,
        bins=bins,
        weights=df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights_cut):0.1f} Events : {len(sig_weights_cut)} Raw Events",
        log=False
    )

    axes[1].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=fs)
    
def draw_plots_3(df, dfs, var, bdt_cut, name, bins):
    fig, axes = plt.subplots(1,3,figsize=(24, 8))

    plotData_background = df.loc[df["sample_isSignal"] == False][var]
    plotData_signal = df.loc[df["sample_isSignal"] == True][var]

    plotData_background_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)][var]
    plotData_signal_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)][var]

    bkg_weights = df.loc[(df["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights = df.loc[(df["sample_isSignal"] == True)]["Xsec_genWeight"]

    bkg_weights_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"] 
    sig_weights_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"]

    axes[0].hist(
        plotData_background,
        bins=bins,
        weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events : {len(plotData_signal)} Raw Events",
        log=False
    )

    axes[0].hist(
        plotData_signal,
        bins=bins,
        weights=df.loc[df["sample_isSignal"] == True]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights):0.0001f} Events : {len(plotData_background)} Raw Events",
        log=False
    )

    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=16)

    axes[1].hist(
        plotData_background_cut,
        bins=bins,
        weights=df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights_cut):0.0001f} Events : {len(plotData_signal_cut)} Raw Events",
        log=False
    )

    axes[1].hist(
        plotData_signal_cut,
        bins=bins,
        weights=df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"],
        histtype="step",
        color="r",
        label=f"Signal : {sum(sig_weights_cut):0.1f} Events : {len(plotData_background_cut)} Raw Events",
        log=False
    )

    axes[1].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=16)

        
    for df1 in dfs:
        Sample_name = df1.iloc[-1].at["Sample_name"]
        df_cut = df1[df1.eval("bdt > "+str(bdt_cut))].copy()
        if df_cut.empty:
            continue
        axes[2].hist(
            df_cut[var],
            bins=bins,
            weights=df_cut["Xsec_genWeight"],
            histtype="step",
            label=f"{Sample_name:s} : {df_cut.Xsec_genWeight.sum():0.0001f} Events : {str(len(df_cut)):s} Raw Events",
            log=False
        ) 

    axes[2].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[2].set_ylabel("Events", size=18)
    axes[2].legend(fontsize=8)
    
def draw_plots_separate_background(dfs, var, bdt_cut, name, bins):
    fig, axes = plt.subplots(1,2,figsize=(16, 8))
    
    for df in dfs:
        Sample_name = df.iloc[-1].at["Sample_name"]
        axes[0].hist(
            df[var],
            bins=bins,
            weights=df["Xsec_genWeight"],
            histtype="step",
            label=f"{Sample_name:s} : {df.Xsec_genWeight.sum():0.1f} Events",
            log=False
        ) 
        
    for df in dfs:
        Sample_name = df.iloc[-1].at["Sample_name"]
        df_cut = df[df.eval("bdt > "+str(bdt_cut))].copy()
        if df_cut.empty:
            continue
        axes[1].hist(
            df_cut[var],
            bins=bins,
            weights=df_cut["Xsec_genWeight"],
            histtype="step",
            label=f"{Sample_name:s} : {df_cut.Xsec_genWeight.sum():0.1f} Events : {str(len(df_cut)):s} Raw Events",
            log=False
        ) 
    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=8)
    
    axes[1].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=8)
    
def draw_bdt_score_quant(df):
    fig, axes = plt.subplots(1,1,figsize=(18, 10))
    
    df_bkg = df.loc[df["sample_isSignal"] == 0][['Xsec_genWeight', 'bdt', 'Sample_name']].copy()
    df_sig = df.loc[df["sample_isSignal"] == 1][['Xsec_genWeight', 'bdt', 'Sample_name']].copy()
    
    df_sig_WWH = df.loc[(df["sample_isSignal"] == 1) & (df["Sample_name"] == "WWH_ucsd_C2V_Reweight")][['Xsec_genWeight', 'bdt', 'Sample_name']].copy()
    df_sig_WZH = df.loc[(df["sample_isSignal"] == 1) & (df["Sample_name"] == "WZH_ucsd_C2V_Reweight")][['Xsec_genWeight', 'bdt', 'Sample_name']].copy()
    
    sig_weight_total = df_sig["Xsec_genWeight"].sum()
    bkg_weight_total = df_bkg["Xsec_genWeight"].sum()
    
    df_bkg["Xsec_genWeight_norm"] = df_bkg["Xsec_genWeight"] / bkg_weight_total
    df_sig["Xsec_genWeight_norm"] = df_sig["Xsec_genWeight"] / sig_weight_total
    
    UniqueNames_bkg = df_bkg.Sample_name.unique()
    UniqueNames_sig = df_sig.Sample_name.unique()
    print(UniqueNames_sig)
    
    quant_bins = [0.0]
    sig_quants = []
    quant_names = ["0"]
    
    df_quant = df.loc[df["sample_isSignal"] == 1][['Xsec_genWeight', 'bdt']].copy()
    #df_quant = df_quant.reset_index(drop=True)
    df_quant = df_quant.sort_values('bdt')
    df_quant = df_quant.reset_index(drop=True)
    #display(df_quant)
    
    sig_sum = 0.0
    for ind in df_quant.index:
        sig_sum += df_quant['Xsec_genWeight'][ind]
        if ind == len(df_quant)-1:
            quant_bins.append(1.0)
            sig_quants.append(sig_sum/sig_weight_total)
            break
        if sig_sum/sig_weight_total >= 0.1:
                quant_bins.append((df_quant['bdt'][ind]))# + df_quant['bdt'][ind])/2)
                quant_names.append(f"{((df_quant['bdt'][ind+1] + df_quant['bdt'][ind])/2):0.4f}")
                sig_quants.append(sig_sum/sig_weight_total)
                #print(sig_sum,sig_weight_total,sig_sum/sig_weight_total)
                sig_sum = 0.0
    #print(quant_bins)
    #print(sig_quants)
    #print(sum(sig_quants))
    
    h,e = np.histogram(df_bkg["bdt"], 
                       bins=quant_bins, 
                       weights=df_bkg["Xsec_genWeight"])
    
    h1,e1 = np.histogram(df_sig["bdt"], 
                       bins=quant_bins, 
                       weights=df_sig["Xsec_genWeight"])
    
    h2,e2 = np.histogram(df_bkg["bdt"], 
                       bins=quant_bins, 
                       weights=df_bkg["Xsec_genWeight_norm"])
    
    h3,e3 = np.histogram(df_sig["bdt"], 
                       bins=quant_bins, 
                       weights=df_sig["Xsec_genWeight_norm"])
    
    hr,er = np.histogram(df_bkg["bdt"], 
                       bins=quant_bins)
    
    hwwh,ewwh = np.histogram(df_sig_WWH["bdt"], 
                       bins=quant_bins, 
                       weights=df_sig_WWH["Xsec_genWeight"])
    hwzh,ewzh = np.histogram(df_sig_WZH["bdt"], 
                       bins=quant_bins, 
                       weights=df_sig_WZH["Xsec_genWeight"])
    
    
    rectswwh = axes.bar(range(len(quant_bins)-1),hwwh, width=1, edgecolor='none',color='none',
                align="edge",)
    rectswzh = axes.bar(range(len(quant_bins)-1),hwzh, width=1, edgecolor='none',color='none',
                align="edge",)
    
    labels_bkgind_raw = {}
    labels_bkgind_scaled = {}
    for name in UniqueNames_bkg:
        df_b_tmp = df_bkg[df_bkg.Sample_name == name].copy()
        htmp,etmp = np.histogram(df_b_tmp["bdt"], bins=quant_bins)
        htmp2,etmp2 = np.histogram(df_b_tmp["bdt"], bins=quant_bins, 
                       weights=df_b_tmp["Xsec_genWeight"])
        rects_tmp = axes.bar(range(len(quant_bins)-1),htmp, width=1, edgecolor='none',color='none',
                    align="edge",)
        rects_tmp2 = axes.bar(range(len(quant_bins)-1),htmp2, width=1, edgecolor='none',color='none',
                    align="edge",)
        labels_tmp = [i.get_height() for i in rects_tmp]
        labels_tmp2 = [f"          {i.get_height():.2e}" for i in rects_tmp2]
        labels_bkgind_raw[name] = labels_tmp
        labels_bkgind_scaled[name] = labels_tmp2
    
    rects1 = axes.bar(range(len(quant_bins)-1),h3, width=1, edgecolor='none',color='none',
                align="edge",)
    
    #rects = axes[0].patches

    # Make some labels.

   # for rect,label in zip(rects1,labels_sig):
        #height = rect.get_height()
        #print(height)
        #axes[0].text(
            #rect.get_x() + rect.get_width() / 2, height, label, ha="center", va="bottom",color="r"
        #)
    
    rects2 = axes.bar(range(len(quant_bins)-1),h2, width=1, edgecolor='none',color='none',
                align="edge",)
    
    labels_bkg = [f"{i.get_height()*100:0.2f}%" for i in rects2]
    print(labels_bkg)
        
    rects3 = axes.bar(range(len(quant_bins)-1),h1, width=1, edgecolor='r', color='none',log=True,
                align="edge",
                tick_label=quant_names,
                label=f"Signal : {sig_weight_total:0.1f} Events")
    
    labels_sig = [f"{i.get_height()*100:0.2f}%" for i in rects1]
    
    #print("c")
    #for rect,label in zip(rects3,labels_sig):
        #height = rect.get_height()
        #print(height)
        #axes.text(
            #rect.get_x() + rect.get_width() / 2, height, label, ha="center", va="bottom",color="r"
        #)
    

    rects4 = axes.bar(range(len(quant_bins)-1),h, width=1, edgecolor='b', color='none',log=True,
                align="edge",
                tick_label=quant_names,
                label=f"Background : {bkg_weight_total:0.1f} Events")
    
    rects5 = axes.bar(range(len(quant_bins)-1),hr, width=1, edgecolor='none', color='none',
                align="edge",)
    
    labels_bkg_total = []
    for i in range(10):
        label = f"{rects2[i].get_height()*100:0.2f}%"
        label += f"\n{rects4[i].get_height():.2e}"
        label += f"\nRaw: {rects5[i].get_height()}\n"
        for name in UniqueNames_bkg:
            if labels_bkgind_raw[name][i] == 0:
                continue
            label += "\n" + name + ": " + str(labels_bkgind_raw[name][i])
            label += "\n" + str(labels_bkgind_scaled[name][i])
        labels_bkg_total.append(label)
        
    
    labels_sig_total = []
    for i in range(10):
        label = f"{rects1[i].get_height()*100:0.2f}%"
        label += f"\nosWWH: {rectswwh[i].get_height():.2e}"
        label += f"\nWZH: {rectswzh[i].get_height():.2e}\n"
        labels_sig_total.append(label)
    
    for rect,label in zip(rects4,labels_bkg_total):
        height = 0.06
        #print(height)
        axes.text(
            rect.get_x(), abs(height), label, ha="left", va="top",color="b"
        )
        
    for rect,label in zip(rects3,labels_sig_total):
        height = rect.get_height()
        #print(height)
        axes.text(
            rect.get_x() + rect.get_width() / 2, abs(height), label, ha="center", va="bottom",color="r"
        )
    

    axes.set_xlabel("BDT Score", size=18)
    axes.set_ylabel("Log(Events)", size=18)
    axes.legend(fontsize=16)
    
    

    fig.tight_layout()
    
    
       
def draw_plots_tmp(df, var, bdt_cut, name, bins,fs=16):
    fig, axes = plt.subplots(1,2,figsize=(16, 8))

    plotData_background = df.loc[df["sample_isSignal"] == False][var]
    plotData_signal = df.loc[df["sample_isSignal"] == True][var]

    plotData_background_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)][var]
    plotData_signal_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)][var]

    bkg_weights = df.loc[(df["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights = df.loc[(df["sample_isSignal"] == True)]["Xsec_genWeight"]

    bkg_weights_cut = df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"] 
    sig_weights_cut = df.loc[(df["sample_isSignal"] == True) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"]

    axes[0].hist(
        plotData_background,
        bins=bins,
        weights=df.loc[df["sample_isSignal"] == False]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events : {len(bkg_weights)} Raw Events",
        log=False
    )
    for sig_name in ["WWH_ucsd_C2V_Reweight"]:
        sig_samp = df[((df.Sample_name == sig_name))].copy()
        plotData_sig_samp = sig_samp[var]
        sig_weights = sig_samp["Xsec_genWeight"]
                    

        axes[0].hist(
            plotData_sig_samp,
            bins=bins,
            weights=sig_samp["Xsec_genWeight"],
            histtype="step",
            #color="r",
            label=f"{sig_name} : {sum(sig_weights):0.1f} Events : {len(sig_weights)} Raw Events",
            log=False
        )

    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=fs)

    axes[1].hist(
        plotData_background_cut,
        bins=bins,
        weights=df.loc[(df["sample_isSignal"] == False) & (df["bdt"] >= bdt_cut)]["Xsec_genWeight"],
        histtype="step",
        label=f"Background : {sum(bkg_weights_cut):0.1f} Events : {len(bkg_weights_cut)} Raw Events",
        log=False
    )
    for sig_name in ["WWH_ucsd_C2V_Reweight"]:
        sig_samp = df[((df.Sample_name == sig_name) & (df["bdt"] >= bdt_cut))].copy()
        plotData_sig_samp = sig_samp[var]
        sig_weights = sig_samp.loc[(sig_samp["bdt"] >= bdt_cut)]["Xsec_genWeight"]
                    

        axes[1].hist(
            plotData_sig_samp,
            bins=bins,
            weights=sig_samp.loc[(sig_samp["bdt"] >= bdt_cut)]["Xsec_genWeight"],
            histtype="step",
            #color="r",
            label=f"{sig_name} : {sum(sig_weights):0.1f} Events : {len(sig_weights)} Raw Events",
            log=False
        )


    axes[1].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=fs)
    
    
      
def draw_plots_tmp_2(df, var, bdt_cut, name, bins,fs=16):
    fig, axes = plt.subplots(1,2,figsize=(16, 8))
    
    df_cut = df[df["bdt"] >= bdt_cut].copy()
    df_wwh = df[(df["sample_isSignal"] == True) & (df['Sample_name'].str.contains('WWH'))].copy()
    df_wzh = df[(df["sample_isSignal"] == True) & (df['Sample_name'].str.contains('WZH'))].copy()

    df_wwh_cut = df_cut[(df_cut["sample_isSignal"] == True) & (df_cut['Sample_name'].str.contains('WWH'))].copy()
    df_wzh_cut = df_cut[(df_cut["sample_isSignal"] == True) & (df_cut['Sample_name'].str.contains('WZH'))].copy()

    
    plotData_background = df.loc[df["sample_isSignal"] == False][var]
    plotData_signal = df.loc[df["sample_isSignal"] == True][var]
    plotData_signal_WWH = df_wwh[var]
    plotData_signal_WZH = df_wzh[var]

    plotData_background_cut = df_cut.loc[(df_cut["sample_isSignal"] == False)][var]
    plotData_signal_cut = df_cut.loc[(df_cut["sample_isSignal"] == True)][var]
    plotData_signal_cut_WWH = df_wwh_cut[var]
    plotData_signal_cut_WZH = df_wzh_cut[var]

    bkg_weights = df.loc[(df["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights = df.loc[(df["sample_isSignal"] == True)]["Xsec_genWeight"]
    sig_weights_WWH = df_wwh["Xsec_genWeight"]
    sig_weights_WZH = df_wzh["Xsec_genWeight"]

    bkg_weights_cut = df_cut.loc[(df_cut["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights_cut = df_cut.loc[(df_cut["sample_isSignal"] == True)]["Xsec_genWeight"]
    sig_weights_cut_WWH = df_wwh_cut["Xsec_genWeight"]
    sig_weights_cut_WZH = df_wzh_cut["Xsec_genWeight"]

    #print(len(plotData_signal_WWH), len(sig_weights_WWH,))
    #print(sum(bkg_weights), sum(sig_weights), sum(sig_weights_WWH),sum(sig_weights_WZH))
    axes[0].hist(
        plotData_background,
        bins=bins,
        weights=bkg_weights,
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events : {len(bkg_weights)} Raw Events",
        #log=True
    )
    axes[0].hist(
        plotData_signal,
        bins=bins,
        weights=sig_weights,
        histtype="step",
        #color="r",
        label=f"Signal : {sum(sig_weights):0.1f} Events : {len(sig_weights)} Raw Events",
        #log=True
    )
    axes[0].hist(
            plotData_signal_WWH,
            bins=bins,
            weights=sig_weights_WWH,
            histtype="step",
            #color="r",
            label=f"WWH : {sum(sig_weights_WWH):0.1f} Events : {len(sig_weights_WWH)} Raw Events",
            #log=True
        )
    axes[0].hist(
            plotData_signal_WZH,
            bins=bins,
            weights=sig_weights_WZH,
            histtype="step",
            #color="r",
            label=f"WZH : {sum(sig_weights_WZH):0.1f} Events : {len(sig_weights_WZH)} Raw Events",
            #log=True
        )

    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=fs)

    axes[1].hist(
        plotData_background_cut,
        bins=bins,
        weights=bkg_weights_cut,
        histtype="step",
        label=f"Background : {sum(bkg_weights_cut):0.1f} Events : {len(bkg_weights_cut)} Raw Events",
        #log=True
    )
    axes[1].hist(
            plotData_signal_cut,
            bins=bins,
            weights=sig_weights_cut,
            histtype="step",
            #color="r",
            label=f"Signal : {sum(sig_weights_cut):0.1f} Events : {len(sig_weights_cut)} Raw Events",
            #log=True
        )
    axes[1].hist(
            plotData_signal_cut_WWH,
            bins=bins,
            weights=sig_weights_cut_WWH,
            histtype="step",
            #color="r",
            label=f"WWH : {sum(sig_weights_cut_WWH):0.1f} Events : {len(sig_weights_cut_WWH)} Raw Events",
            #log=True
        )
    axes[1].hist(
            plotData_signal_cut_WZH,
            bins=bins,
            weights=sig_weights_cut_WZH,
            histtype="step",
            #color="r",
            label=f"WZH : {sum(sig_weights_cut_WZH):0.1f} Events : {len(sig_weights_cut_WZH)} Raw Events",
            #log=True
        )



    axes[1].set_xlabel(name + " (BDT Score > "+str(bdt_cut)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=fs) 

    
 
def draw_plots_tmp_3(df, df_2,df_3, var, bdt_cut, bdt_cut_2,bdt_cut_3, name, bins,fs=16):
    fig, axes = plt.subplots(1,2,figsize=(16, 8))
    
    df_full = pd.concat([df,df_2,df_3])
    df_cut_1 = df[df["bdt"] >= bdt_cut].copy()
    df_cut_2 = df_2[df_2["bdt"] >= bdt_cut_2].copy()
    df_cut_3 = df_3[df_3["bdt"] >= bdt_cut_3].copy()
    df_full_cut = pd.concat([df_cut_1,df_cut_2,df_cut_3])
    df_wwh = df_full[(df_full["sample_isSignal"] == True) & (df_full['Sample_name'].str.contains('WWH'))].copy()
    df_wzh = df_full[(df_full["sample_isSignal"] == True) & (df_full['Sample_name'].str.contains('WZH'))].copy()

    df_wwh_cut = df_full_cut[(df_full_cut["sample_isSignal"] == True) & (df_full_cut['Sample_name'].str.contains('WWH'))].copy()
    df_wzh_cut = df_full_cut[(df_full_cut["sample_isSignal"] == True) & (df_full_cut['Sample_name'].str.contains('WZH'))].copy()

    
    plotData_background = df_full.loc[df_full["sample_isSignal"] == False][var]
    plotData_signal = df_full.loc[df_full["sample_isSignal"] == True][var]
    plotData_signal_WWH = df_wwh[var]
    plotData_signal_WZH = df_wzh[var]
    
    unique = (df_full_cut[(df_full_cut["sample_isSignal"] == False)]["Sample_name"].unique())
    disp = {}
    for i in unique:
        disp[i] = df_full_cut[df_full_cut["Sample_name"]==i]["Xsec_genWeight"].sum()
    print(disp)

    plotData_background_cut = df_full_cut.loc[(df_full_cut["sample_isSignal"] == False)][var]
    plotData_signal_cut = df_full_cut.loc[(df_full_cut["sample_isSignal"] == True)][var]
    plotData_signal_cut_WWH = df_wwh_cut[var]
    plotData_signal_cut_WZH = df_wzh_cut[var]

    bkg_weights = df_full.loc[(df_full["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights = df_full.loc[(df_full["sample_isSignal"] == True)]["Xsec_genWeight"]
    sig_weights_WWH = df_wwh["Xsec_genWeight"]
    sig_weights_WZH = df_wzh["Xsec_genWeight"]

    bkg_weights_cut = df_full_cut.loc[(df_full_cut["sample_isSignal"] == False)]["Xsec_genWeight"] 
    sig_weights_cut = df_full_cut.loc[(df_full_cut["sample_isSignal"] == True)]["Xsec_genWeight"]
    sig_weights_cut_WWH = df_wwh_cut["Xsec_genWeight"]
    sig_weights_cut_WZH = df_wzh_cut["Xsec_genWeight"]

    #print(len(plotData_signal_WWH), len(sig_weights_WWH,))
    #print(sum(bkg_weights), sum(sig_weights), sum(sig_weights_WWH),sum(sig_weights_WZH))
    axes[0].hist(
        plotData_background,
        bins=bins,
        weights=bkg_weights,
        histtype="step",
        label=f"Background : {sum(bkg_weights):0.1f} Events : {len(bkg_weights)} Raw Events",
        #log=True
    )
    axes[0].hist(
        plotData_signal,
        bins=bins,
        weights=sig_weights,
        histtype="step",
        #color="r",
        label=f"Signal : {sum(sig_weights):0.1f} Events : {len(sig_weights)} Raw Events",
        #log=True
    )
    axes[0].hist(
            plotData_signal_WWH,
            bins=bins,
            weights=sig_weights_WWH,
            histtype="step",
            #color="r",
            label=f"WWH : {sum(sig_weights_WWH):0.1f} Events : {len(sig_weights_WWH)} Raw Events",
            #log=True
        )
    axes[0].hist(
            plotData_signal_WZH,
            bins=bins,
            weights=sig_weights_WZH,
            histtype="step",
            #color="r",
            label=f"WZH : {sum(sig_weights_WZH):0.1f} Events : {len(sig_weights_WZH)} Raw Events",
            #log=True
        )

    axes[0].set_xlabel(name, size=18)
    axes[0].set_ylabel("Events", size=18)
    axes[0].legend(fontsize=fs)

    axes[1].hist(
        plotData_background_cut,
        bins=bins,
        weights=bkg_weights_cut,
        histtype="step",
        label=f"Background : {sum(bkg_weights_cut):0.1f} Events : {len(bkg_weights_cut)} Raw Events",
        #log=True
    )
    axes[1].hist(
            plotData_signal_cut,
            bins=bins,
            weights=sig_weights_cut,
            histtype="step",
            #color="r",
            label=f"Signal : {sum(sig_weights_cut):0.1f} Events : {len(sig_weights_cut)} Raw Events",
            #log=True
        )
    axes[1].hist(
            plotData_signal_cut_WWH,
            bins=bins,
            weights=sig_weights_cut_WWH,
            histtype="step",
            #color="r",
            label=f"WWH : {sum(sig_weights_cut_WWH):0.1f} Events : {len(sig_weights_cut_WWH)} Raw Events",
            #log=True
        )
    axes[1].hist(
            plotData_signal_cut_WZH,
            bins=bins,
            weights=sig_weights_cut_WZH,
            histtype="step",
            #color="r",
            label=f"WZH : {sum(sig_weights_cut_WZH):0.1f} Events : {len(sig_weights_cut_WZH)} Raw Events",
            #log=True
        )



    axes[1].set_xlabel(name + " (EE BDT Score > "+str(bdt_cut)+", EMu BDT Score > "+str(bdt_cut_2)+", MuMu BDT Score > "+str(bdt_cut_3)+")", size=18)
    axes[1].set_ylabel("Events", size=18)
    axes[1].legend(fontsize=fs) 
