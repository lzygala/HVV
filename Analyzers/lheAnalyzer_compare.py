#!/usr/bin/env python
import sys, math, ROOT, copy

def KILL(log):
    print '\n@@@ FATAL -- '+log+'\n'
    raise SystemExit

def lhep_pdgID  (line): return int  (line.split()[ 0])
def lhep_status (line): return int  (line.split()[ 1])
def lhep_mother1(line): return int  (line.split()[ 2])
def lhep_mother2(line): return int  (line.split()[ 3])
def lhep_px     (line): return float(line.split()[ 6])
def lhep_py     (line): return float(line.split()[ 7])
def lhep_pz     (line): return float(line.split()[ 8])
def lhep_E      (line): return float(line.split()[ 9])
def lhep_M      (line): return float(line.split()[10])

def print_lhep(l):
    print lhep_pdgID  (l),
    print lhep_status (l),
    print lhep_mother1(l),
    print lhep_mother2(l),
    print lhep_px     (l),
    print lhep_py     (l),
    print lhep_pz     (l),
    print lhep_E      (l),
    print lhep_M      (l)

    return

def draw_plots(hist_list, plot_name):
    canvas = ROOT.TCanvas("c","c",3600,2400)
    i = 1
    maximum = 0

    legend = ROOT.TLegend(0.7,0.7,0.9,0.9,"","NDC")
    ROOT.gStyle.SetOptStat(0)

    for label in hist_list:
        hist_list[label].Scale(1./hist_list[label].Integral())
        if hist_list[label].GetMaximum() > maximum:
            maximum = hist_list[label].GetMaximum()


    for label in hist_list:
        hist_list[label].SetLineColor(i)
        #hist_list[label].SetFillColorAlpha(i,0.1)

        legend.AddEntry(hist_list[label], label, "l")

        if i == 1:
            hist_list[label].SetTitle(plot_name)
            hist_list[label].SetMaximum(maximum*1.01)
            hist_list[label].Draw("hist")
        else:
            hist_list[label].Draw("hist, sames")

        i += 1

    legend.Draw("same")

    canvas.SaveAs(plot_name+'.png')
    canvas.SaveAs(plot_name+'.pdf')


### main
if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)


    ## define relevant histograms
    h_higgs_m = {}
    h_higgs_pt = {}
    h_wp_m = {}
    h_wp_pt = {}
    h_wm_m = {}
    h_wm_pt = {}
    h_j1_pt = {}
    h_j2_pt = {}
    h_jet_invmass = {}
    h_jet_invmass_extended = {}
    h_j1_eta = {}
    h_j2_eta = {}
    h_jet_etasep = {}

    files = { ("chw_0.01000000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_0.01000000/unweighted_events_10k_merged.lhe", 'r')),
              ("chw_0.5000000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_0.5000000/unweighted_events_10k_merged.lhe", 'r')),
              ("chw_0.9999000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_0.9999000/unweighted_events_10k_merged.lhe", 'r')),
              ("chw_1.000000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_1.000000/unweighted_events_10k_merged.lhe", 'r')),
              ("chw_5.000000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_5.000000/unweighted_events_10k_merged.lhe", 'r')),
              ("chw_10.000000", open("/eos/cms/store/user/lzygala/HVV/pp_hw+w-jj_5k_SMNP0_EFT_VBFCut/chw_10.000000/unweighted_events_10k_merged.lhe", 'r'))}

    for k, ifile in files: # one histo per weight class
        label = k

        h_higgs_m[k]  = ROOT.TH1F('h_higgs_m_'+label , 'h_higgs_m_'+label ,100, 0, 2000)
        h_higgs_pt[k]  = ROOT.TH1F('h_higgs_pt_'+label , 'h_higgs_pt_'+label ,100, 0, 2000)
        h_wp_m[k]  = ROOT.TH1F('h_wp_m_'+label , 'h_wp_m_'+label ,100, 0, 2000)
        h_wp_pt[k]  = ROOT.TH1F('h_wp_pt_'+label , 'h_wp_pt_'+label ,100, 0, 2000)
        h_wm_m[k]  = ROOT.TH1F('h_wm_m_'+label , 'h_wm_m_'+label ,100, 0, 2000)
        h_wm_pt[k]  = ROOT.TH1F('h_wm_pt_'+label , 'h_wm_pt_'+label ,100, 0, 2000)
        h_j1_pt[k]  = ROOT.TH1F('h_j1_pt_'+label , 'h_j1_pt_'+label ,100, 0, 2000)
        h_j2_pt[k]  = ROOT.TH1F('h_j2_pt_'+label , 'h_j2_pt_'+label ,100, 0, 2000)
        h_j1_eta[k]  = ROOT.TH1F('h_j1_eta_'+label , 'h_j1_eta_'+label ,100, -5, 5)
        h_j2_eta[k]  = ROOT.TH1F('h_j2_eta_'+label , 'h_j2_eta_'+label ,100, -5, 5)
        h_jet_invmass[k]  = ROOT.TH1F('h_jet_invmass_'+label , 'h_jet_invmass_'+label ,120, 0, 120)
        h_jet_invmass_extended[k]  = ROOT.TH1F('h_jet_invmass_extended_'+label , 'h_jet_invmass_extended_'+label ,100, 0, 2000)
        h_jet_etasep[k]  = ROOT.TH1F('h_jet_etasep_'+label , 'h_jet_etasep_'+label ,100, 0, 10)



    ###
    for k, ifile in files:

        event_num_max = -1
        event_num, in_event = 0, False

        # reads the lhe and looks into the events
        for line in ifile:
            if line[:1] == '#': continue
            if line.startswith('<scales'): continue

            if event_num_max > 0:
                if event_num > event_num_max: continue

            if line.startswith('<event>'):
                event_num += 1

                genp_ls = []
                weight = {}
                in_event = True
                continue

            if in_event:

                if not line.startswith('</event>'):
                    l0 = line.strip('\n')
                    
                    if l0.startswith('<wgt'):
                        l1 = l0.split()
                        weight[l1[1].split("=")[1].strip("'>")] = float(l1[2])
                        continue

                    if l0.startswith('<'): continue
                    if len(l0.split()) == 6: 
                        continue

                    genp_ls.append(l0)

                else:
                    #print("New Event")
                    ### event analysis

                # define the four momentum of the dilepton pair. 

                    h_p4 = ROOT.TLorentzVector(0, 0, 0, 0)
                    wp_p4 = ROOT.TLorentzVector(0, 0, 0, 0)
                    wm_p4 = ROOT.TLorentzVector(0, 0, 0, 0)
                    j1_p4 = ROOT.TLorentzVector(0, 0, 0, 0)
                    j2_p4 = ROOT.TLorentzVector(0, 0, 0, 0)

                    j1_flag = False

                    for p in genp_ls:

                # for each particle in an event: extract four momentum
                        i_p4 = ROOT.TLorentzVector(lhep_px(p), lhep_py(p), lhep_pz(p), lhep_E(p))
            
                        if lhep_status(p) == 1: 

                # e or mu
                            if lhep_pdgID(p) == 24: 
                                #print("Found wp")
                                wp_p4 += i_p4

                            if lhep_pdgID(p) == -24: 
                                #print("Found wm")
                                wm_p4 += i_p4

                            if abs(lhep_pdgID(p)) == 25: 
                                #print("Found higgs")
                                h_p4 += i_p4

                            if j1_flag == False and (abs(lhep_pdgID(p)) == 1 or abs(lhep_pdgID(p)) ==2 or abs(lhep_pdgID(p)) ==3 or abs(lhep_pdgID(p)) ==4 or abs(lhep_pdgID(p)) ==5 or abs(lhep_pdgID(p)) ==6): 
                                j1_p4 += i_p4
                                j1_flag = True
                                #print("Found j1")

                            elif abs(lhep_pdgID(p)) == 1 or abs(lhep_pdgID(p)) ==2 or abs(lhep_pdgID(p)) ==3 or abs(lhep_pdgID(p)) ==4 or abs(lhep_pdgID(p)) ==5 or abs(lhep_pdgID(p)) ==6: 
                                j2_p4 += i_p4
                                #print("Found j2")

                    jet_sum = j1_p4 + j2_p4
                    inv_mass = math.sqrt(jet_sum*jet_sum)
                    
                    #if inv_mass > 100:
                    #    continue
                    # for each event: store the observables in the histograms
                    h_higgs_pt[k].Fill(h_p4.Pt()) 
                    h_wp_pt[k].Fill(wp_p4.Pt()) 
                    h_wm_pt[k].Fill(wm_p4.Pt()) 
                    h_j1_pt[k].Fill(j1_p4.Pt()) 
                    h_j2_pt[k].Fill(j2_p4.Pt()) 
                    h_jet_invmass[k].Fill(inv_mass) 
                    h_jet_invmass_extended[k].Fill(inv_mass)
                    h_j1_eta[k].Fill(j1_p4.Eta()) 
                    h_j2_eta[k].Fill(j2_p4.Eta()) 
                    h_jet_etasep[k].Fill(abs(j1_p4.Eta()-j2_p4.Eta())) 

                        #print(abs(j1_p4.Eta()-j2_p4.Eta()))
                        

                    in_event = False
                    continue


        print("processed %i events" %event_num)

    draw_plots(h_higgs_pt, 'h_higgs_pt')
    draw_plots(h_wp_pt, 'h_wp_pt')
    draw_plots(h_wm_pt, 'h_wm_pt')
    draw_plots(h_j1_pt, 'h_j1_pt')
    draw_plots(h_j2_pt, 'h_j2_pt')
    draw_plots(h_jet_invmass, 'h_jet_invmass')
    draw_plots(h_jet_invmass_extended, 'h_jet_invmass_extended')
    draw_plots(h_j1_eta, 'h_j1_eta')
    draw_plots(h_j2_eta, 'h_j2_eta')
    draw_plots(h_jet_etasep, 'h_jet_etasep')

  