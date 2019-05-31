def getMT(doprint=True):
 #       List of genes
 #         List of mutations
 #           List of groups of variables
 #             Pair (variable list + range)
 #               List of variables
 MT = [] 
 #CACNA1C:
 MT.append([ [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -25.9 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19265197
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -27.0 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -37.3 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19265197
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -30.0 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-31.4, +7.0] ],               #http://www.ncbi.nlm.nih.gov/pubmed/21685391                                           
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.85, 1.45] ],
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-28.5, +16.3] ],
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], [0.72, 1.38] ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-38.5, +12.9] ],              #http://www.ncbi.nlm.nih.gov/pubmed/16157588                                           
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.46, 1.56] ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-27.8, +8.7] ],               #http://www.ncbi.nlm.nih.gov/pubmed/18836301                                           
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.89, 1.14] ],
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-19.1, +4.7] ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-11.2, +1.0] ],               #http://www.ncbi.nlm.nih.gov/pubmed/15299022                                           
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -3.1 ],
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 1.24 ] ] ])

 #CACNB2:
 MT.append([ [ [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -5.2 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/19358333
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 0.69 ] ], 
             [ [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 1.7] ],                        #http://www.ncbi.nlm.nih.gov/pubmed/7723731
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-4.9, 4.9] ],                 #http://www.ncbi.nlm.nih.gov/pubmed/19723630
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-5.1, 5.1] ],
               [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], [0.6, 1.68] ],
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], [0.6, 1.66] ] ],
             [ [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 1.26] ] ])                     #http://www.ncbi.nlm.nih.gov/pubmed/20025708

 #CACNA1D:
 MT.append([ [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -10.9 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], 0.73 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/21998310
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-3.0, 3.5] ],                 #(42A)
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 0.81 ],
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 1.25 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-10.6, 3.4] ],                #http://www.ncbi.nlm.nih.gov/pubmed/21998309 and
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.8, 1.12] ],                 #http://www.ncbi.nlm.nih.gov/pubmed/21998310
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-5.3, 1.2] ],                 #(43S)
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 0.66 ],
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 0.72 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], 6.6 ],                         #http://www.ncbi.nlm.nih.gov/pubmed/20951705 and
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.75, 1.19] ],                #http://www.ncbi.nlm.nih.gov/pubmed/21054386
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], [0.5, 1.12] ] ],               #(CaV1.3 KO)
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -9.8 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/25620733
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], 0.8 ],
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -15.4 ],
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 1.05 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-24.2, +6.1] ],               #http://www.ncbi.nlm.nih.gov/pubmed/23913004
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.7, 1.24] ],
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], -14.5 ],
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], [0.72, 1.28] ],
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 3.52 ] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -17.8 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/22760075
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], 0.81 ],
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], [0.77, 1.31] ] ] ])

 #CACNA1I:
 MT.append([ [ [ 'offma_Ca_LVAst', 1.3 ],                                         #http://www.ncbi.nlm.nih.gov/pubmed/15254077
               [ 'offha_Ca_LVAst', 1.6 ],
               [ ['taummin_Ca_LVAst', 'taumdiff_Ca_LVAst'], [0.87, 1.45] ],
               [ ['tauhmin_Ca_LVAst', 'tauhdiff_Ca_LVAst'], 0.8 ] ],
             [ [ 'offma_Ca_LVAst', -4.3 ],                                        #http://www.ncbi.nlm.nih.gov/pubmed/12080115
               [ 'sloma_Ca_LVAst', 1.14 ],
               [ 'offha_Ca_LVAst', -4.4 ],
               [ 'sloha_Ca_LVAst', [0.89, 1.04] ],
               [ ['taummin_Ca_LVAst', 'taumdiff_Ca_LVAst'], 0.53 ],
               [ ['tauhmin_Ca_LVAst', 'tauhdiff_Ca_LVAst'], 0.46 ] ] ])

 #CACNA1S:
 MT.append([ [ [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], 0.67 ] ],                      #http://www.ncbi.nlm.nih.gov/pubmed/20861472
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -30.02 ],                      #http://www.ncbi.nlm.nih.gov/pubmed/19134469
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], 0.62 ],
               [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], 0.49] ],
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -4.4 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/24240197
               [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], 0.95 ],
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], 20.6 ],
               [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], 0.97] ] ])

 #ATP2A:
 MT.append([ [ [ 'gamma_CaDynamics_E2', 0.6 ] ] ])                                #http://www.ncbi.nlm.nih.gov/pubmed/10970890
 #MT = [MT, { { { { 'Pup_SERCA', 0.66 ] ],                                        #http://www.ncbi.nlm.nih.gov/pubmed/9891028
 #            [ [ [ 'Pup_SERCA', [0.22, 2.31] ] ],                                #http://www.ncbi.nlm.nih.gov/pubmed/12975374
 #            [ [ [ 'Pup_SERCA', 0.2 ] ] ])                                       #http://www.ncbi.nlm.nih.gov/pubmed/12670936             
 
 #ATP2B:
 MT.append([ [ [ 'decay_CaDynamics_E2', 1.97 ] ],                                 #http://www.ncbi.nlm.nih.gov/pubmed/22789621 
             [ [ 'decay_CaDynamics_E2', 1.5 ],                                    #http://www.ncbi.nlm.nih.gov/pubmed/21232211
               [ 'minCai_CaDynamics_E2' , 1.4 ] ],
             [ [ 'decay_CaDynamics_E2', 4.45 ] ],                                 #http://www.ncbi.nlm.nih.gov/pubmed/17234811
             [ [ 'minCai_CaDynamics_E2', 1.1 ] ] ])                               #http://www.ncbi.nlm.nih.gov/pubmed/22047666

 #NRGN:
 MT.append([ [ [ 'gamma_CaDynamics_E2', 0.4 ] ] ])                                #http://www.ncbi.nlm.nih.gov/pubmed/15564582

 #SCN1A:
 MT.append([ [ [ 'offm_NaTa_t', -0.3 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/18632931
               [ 'offh_NaTa_t', 5 ],
               [ 'slom_NaTa_t', 1.15 ],
               [ 'sloh_NaTa_t', 1.23 ] ],
             [ [ 'offm_NaTa_t', 2.8 ],                                            #http://www.ncbi.nlm.nih.gov/pubmed/18632931
               [ 'offh_NaTa_t', 9.6 ],
               [ 'slom_NaTa_t', 0.984 ],
               [ 'sloh_NaTa_t', 1.042 ] ],
             [ [ 'offm_NaTa_t', -4.0 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/21864321                                           
               [ 'offh_NaTa_t', -5.8 ],
               [ 'slom_NaTa_t', 0.92 ],
               [ 'sloh_NaTa_t', 1.13 ],
               [ ['tauha_NaTa_t', 'tauhb_NaTa_t'], 1.47 ] ],
             [ [ 'offm_NaTa_t', -8.1 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/21864321                                           
               [ 'offh_NaTa_t', 2.2 ],
               [ 'slom_NaTa_t', 0.97 ],
               [ 'sloh_NaTa_t', 0.97 ],
               [ ['tauha_NaTa_t', 'tauhb_NaTa_t'], 1.59 ] ],
             [ [ 'offm_NaTa_t', 6.0 ],                                            #http://www.ncbi.nlm.nih.gov/pubmed/23398611                                           
               [ 'slom_NaTa_t', 1.16 ],
               [ ['tauha_NaTa_t', 'tauhb_NaTa_t'], 1.29 ] ],
             [ [ 'offm_NaTa_t', 10.0 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/16326807                                           
               [ 'offh_NaTa_t', -0.6 ],
               [ 'slom_NaTa_t', 1.15 ],
               [ 'sloh_NaTa_t', 1.14 ] ] ])

 #SCN9A:
 MT.append([ [ [ ['offh_Nap_Et2', 'offha_Nap_Et2', 'offhb_Nap_Et2'], 6.8 ] ],     #http://www.ncbi.nlm.nih.gov/pubmed/22136189
             [ [ ['offh_Nap_Et2', 'offha_Nap_Et2', 'offhb_Nap_Et2'], 3.5 ],       #http://www.ncbi.nlm.nih.gov/pubmed/18945915
               [ 'sloh_Nap_Et2', 0.55 ],
               [ 'offm_NaTa_t', -7.1 ],
               [ 'offh_NaTa_t', 17.0 ],
               [ 'sloh_NaTa_t', 0.69 ] ],
             [ [ 'offm_NaTa_t', -9.1 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/16392115
               [ 'offh_NaTa_t', 3.1 ] ],
             [ [ 'offm_NaTa_t', -7.6 ],                                           #http://www.ncbi.nlm.nih.gov/pubmed/15958509
               [ 'offh_NaTa_t', 4.3 ] ] ])
 #KCNS3:
 MT.append([ [ [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 2.0 ],  #http://www.ncbi.nlm.nih.gov/pubmed/10484328
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 2.5 ],
               [ 'sloh_K_Pst', 0.5 ] ] ])
 #KCNN3:
 MT.append([ [ [ 'offc_SK_E2', 0.86 ],                                            #http://www.ncbi.nlm.nih.gov/pubmed/14978258
               [ 'sloc_SK_E2', 1.24 ] ],
             [ [ 'sloc_SK_E2', 0.84 ] ] ])                                        #http://www.ncbi.nlm.nih.gov/pubmed/17167222
 #HCN1:
 MT.append([ [ [ ['offma_Ih', 'offmb_Ih'], -26.5 ],                               #http://www.ncbi.nlm.nih.gov/pubmed/17185333
               [ ['sloma_Ih', 'slomb_Ih'], 0.64 ] ],
             [ [ ['offma_Ih', 'offmb_Ih'], [-25.9, 17.7] ],                       #http://www.ncbi.nlm.nih.gov/pubmed/12668666
               [ ['sloma_Ih', 'slomb_Ih'], 0.6 ] ],
             [ [ ['offma_Ih', 'offmb_Ih'], 3.9 ],                #http://www.ncbi.nlm.nih.gov/pubmed/26578877  
               [ ['tauma_Ih', 'taumb_Ih'], 0.88 ] ] ])
 #KCNB1:
 MT.append([ [ [ 'offm_K_Pst', 5 ],                                               #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203K)
               [ 'offh_K_Pst', 3 ],
               [ 'slom_K_Pst', 1.11 ],
               [ 'sloh_K_Pst', 0.86 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 0.5 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 0.53 ] ],
             [ [ 'offm_K_Pst', 1 ],                                               #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203D)
               [ 'offh_K_Pst', -6 ],
               [ 'slom_K_Pst', 1.22 ],
               [ 'sloh_K_Pst', 1.0 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 0.89 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 1.13 ] ],
             [ [ 'offm_K_Pst', 6 ],                                               #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347K)
               [ 'offh_K_Pst', -8 ],
               [ 'slom_K_Pst', 1.33 ],
               [ 'sloh_K_Pst', 1.0 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 0.5 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 0.87 ] ],
             [ [ 'offm_K_Pst', -28 ],                                             #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347D)
               [ 'offh_K_Pst', -27 ],
               [ 'slom_K_Pst', 1.11 ],
               [ 'sloh_K_Pst', 0.71 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 1.13 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 2.27 ] ], 
             [ [ 'offm_K_Pst', 14 ],                                              #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (T203W)
               [ 'offh_K_Pst', -21 ],
               [ 'slom_K_Pst', 2.0 ],
               [ 'sloh_K_Pst', 1.0 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 0.39 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 1.2 ] ], 
             [ [ 'offm_K_Pst', -13 ],                                             #http://www.ncbi.nlm.nih.gov/pubmed/21455829 (S347W)
               [ 'offh_K_Pst', -13 ],
               [ 'slom_K_Pst', 1.33 ],
               [ 'sloh_K_Pst', 0.71 ],
               [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], 0.95 ],
               [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], 5.13 ] ] ])
 #CACNB2 reprise:
 MT.append([ [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], 3.00 ],                        #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N1 vs N4)
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], 3.48 ],
               [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], 1.01 ],                        #(tau_h_slow considered)
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 0.89 ] ],  
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -1.11 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N3 vs N4)
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], 5.14 ],
               [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], 0.6 ],                         #(tau_h_slow considered)
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 1.48 ] ], 
             [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], -1.88 ],                       #http://www.ncbi.nlm.nih.gov/pubmed/19723630 (N5 vs N4)
               [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], 2.69 ],
               [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], 0.6 ],                         #(tau_h_slow considered)
               [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], 1.35 ] ] ]) 

 MTgenes = ['CACNB2', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'HCN1', 'KCNB1', 'CACNB2'] 

 if doprint:
  print "MT loaded successfully:"
  print MT
 return MT


def getparamwiseMT():
 MT = [ [ ['offma_Ca_HVA', 'offmb_Ca_HVA'], [-37.3, 12.9] ],
        [ ['offha_Ca_HVA', 'offhb_Ca_HVA'], [-30.0, 20.6] ],
        [ ['sloma_Ca_HVA', 'slomb_Ca_HVA'], [0.46, 1.56] ],
        [ ['sloha_Ca_HVA', 'slohb_Ca_HVA'], [0.66, 1.38] ],
        [ ['tauma_Ca_HVA', 'taumb_Ca_HVA'], [0.49, 1.7] ],
        [ ['tauha_Ca_HVA', 'tauhb_Ca_HVA'], [0.5, 3.52] ],
        [ ['offma_Ca_LVAst'], [-4.3, 1.3] ],
        [ ['offha_Ca_LVAst'], [-4.4, 1.6] ],
        [ ['sloma_Ca_LVAst'], [1.14, 1.14] ],
        [ ['sloha_Ca_LVAst'], [0.89, 1.04] ],
        [ ['taummin_Ca_LVAst', 'taumdiff_Ca_LVAst'], [0.53, 1.45] ],
        [ ['tauhmin_Ca_LVAst', 'tauhdiff_Ca_LVAst'], [0.46, 0.8] ],
        [ ['gamma_CaDynamics_E2'], [0.4, 0.6] ],
        [ ['decay_CaDynamics_E2'], [1.5, 4.45] ],
        [ ['minCai_CaDynamics_E2'], [1.1, 1.4] ],
        [ ['offm_NaTa_t'], [-9.1, 10.0] ],
        [ ['offh_NaTa_t'], [-5.8, 17.0] ],
        [ ['slom_NaTa_t'], [0.92, 1.16] ],
        [ ['sloh_NaTa_t'], [0.69, 1.23] ],
        [ ['tauha_NaTa_t', 'tauhb_NaTa_t'], [1.29, 1.59] ],
        [ ['offh_Nap_Et2', 'offha_Nap_Et2', 'offhb_Nap_Et2'], [3.5, 6.8] ],
        [ ['sloh_Nap_Et2'], [0.55, 0.55] ],
        [ ['offm_K_Pst'], [-28.0, 14.0] ],
        [ ['offh_K_Pst'], [-27.0, 3.0] ],
        [ ['slom_K_Pst'], [1.11, 2.0] ],
        [ ['sloh_K_Pst'], [0.5, 1.0] ],
        [ ['taummin_K_Pst', 'taumdiff1_K_Pst', 'taumdiff2_K_Pst'], [0.39, 2.0] ],
        [ ['tauhmean_K_Pst', 'tauhdiff1_K_Pst', 'tauhdiff2_K_Pst'], [0.53, 2.5] ],
        [ ['offc_SK_E2'], [0.86, 0.86] ],
        [ ['sloc_SK_E2'], [0.84, 1.24] ],
        [ ['offma_Ih', 'offmb_Ih'], [-26.5, 17.7] ],
        [ ['sloma_Ih', 'slomb_Ih'], [0.6, 0.64] ],
        [ ['tauma_Ih', 'taumb_Ih'], [0.88, 0.88] ] ]
 print "paramwiseMT loaded successfully:"
 print MT
 return MT
          
def getdefvals(doprint=True):
 defVals = {'gamma_CaDynamics_E2': [0.000501, 0.000509], #if two values, the first is for somatic and second for apical segments
           'decay_CaDynamics_E2': [460.0, 122.0],
           'depth_CaDynamics_E2': 0.1,
           'minCai_CaDynamics_E2': 1e-4,
           'offma_Ca_HVA': -27.0,
           'offmb_Ca_HVA': -75.0,
           'offha_Ca_HVA': -13.0,
           'offhb_Ca_HVA': -15.0,
           'sloma_Ca_HVA': 3.8,
           'slomb_Ca_HVA': 17.0,
           'sloha_Ca_HVA': 50.0,
           'slohb_Ca_HVA': 28.0,
           'tauma_Ca_HVA': 1.0/0.055,
           'taumb_Ca_HVA': 1.0/0.94,
           'tauha_Ca_HVA': 1.0/0.000457,
           'tauhb_Ca_HVA': 1.0/0.0065,
           'offma_Ca_LVAst': -40.0,
           'offmt_Ca_LVAst': -35.0,
           'offha_Ca_LVAst': -90.0,
           'offht_Ca_LVAst': -50.0,
           'sloma_Ca_LVAst': 6.0,
           'slomt_Ca_LVAst': 5.0,
           'sloha_Ca_LVAst': 6.4,
           'sloht_Ca_LVAst': 7.0,
           'taummin_Ca_LVAst': 5.0,
           'taumdiff_Ca_LVAst': 20.0,
           'tauhmin_Ca_LVAst': 20.0,
           'tauhdiff_Ca_LVAst': 50.0,
           'ehcn_Ih': -45.0,
           'offma_Ih': -154.9,
           'sloma_Ih': 11.9,
           'tauma_Ih': 1.0/0.00643,
           'offmb_Ih': 0.0,
           'slomb_Ih': 33.1,
           'taumb_Ih': 1.0/0.193,
           'offma_Im': -35.0,
           'sloma_Im': 10.0,
           'tauma_Im': 1.0/3.3e-3,
           'offmb_Im': -35.0,
           'slomb_Im': 10.0,
           'taumb_Im': 1.0/3.3e-3,
           'offm_K_Pst': -11.0,
           'slom_K_Pst': 12.0,
           'offmt_K_Pst': -10.0,
           'slomt_K_Pst': 1.0/0.026,
           'taummin_K_Pst': 1.25,
           'taumdiff1_K_Pst': 175.03,
           'taumdiff2_K_Pst': 13.0,
           'offh_K_Pst': -64.0,
           'sloh_K_Pst': 11.0,
           'offht1_K_Pst': -65.0,
           'offht2_K_Pst': -85.0,
           'sloht_K_Pst': 48.0,
           'tauhmean_K_Pst': 360.0,
           'tauhdiff1_K_Pst': 1010.0,
           'tauhdiff2_K_Pst': 24.0,
           'offm_K_Tst': -10.0,
           'slom_K_Tst': 19.0,
           'offh_K_Tst': -76.0,
           'sloh_K_Tst': 10.0,
           'offmt_K_Tst': -81.0,
           'slomt_K_Tst': 59.0,
           'taummin_K_Tst': 0.34,
           'taumdiff_K_Tst': 0.92,
           'offht_K_Tst': -83.0,
           'sloht_K_Tst': 23.0,
           'tauhmin_K_Tst': 8.0,
           'tauhdiff_K_Tst': 49.0,
           'offm_NaTa_t': -38.0,
           'offh_NaTa_t': -66.0,
           'slom_NaTa_t': 6.0,
           'sloh_NaTa_t': 6.0,
           'tauma_NaTa_t': 1.0/0.182,
           'taumb_NaTa_t': 1.0/0.124,
           'tauha_NaTa_t': 1.0/0.015,
           'tauhb_NaTa_t': 1.0/0.015,
           'offm_Nap_Et2': -52.6,
           'slom_Nap_Et2': 4.6,
           'offma_Nap_Et2': -38.0,
           'offmb_Nap_Et2': -38.0,
           'sloma_Nap_Et2': 6.0,
           'slomb_Nap_Et2': 6.0,
           'tauma_Nap_Et2': 1.0/0.182,
           'taumb_Nap_Et2': 1.0/0.124,
           'taummax_Nap_Et2': 6.0,
           'offh_Nap_Et2': -48.8,
           'sloh_Nap_Et2': 10.0,
           'offha_Nap_Et2': -17.0,
           'offhb_Nap_Et2': -64.4,
           'sloha_Nap_Et2': 4.63,
           'slohb_Nap_Et2': 2.63,
           'tauha_Nap_Et2': 1.0/2.88e-6,
           'tauhb_Nap_Et2': 1.0/6.94e-6,
           'tauhmax_Nap_Et2': 1.0,
           'zTau_SK_E2': 1.0,
           'offc_SK_E2': 0.00043,
           'sloc_SK_E2': 4.8,
           'offma_SKv3_1': 18.7,
           'offmt_SKv3_1': -46.56,
           'sloma_SKv3_1': 9.7,
           'slomt_SKv3_1': 44.14,
           'taummax_SKv3_1': 4.0}          
 if doprint:
  print "defVals loaded successfully:"
  print defVals
 return defVals

def getgenenames():
 return ['CACNA1C', 'CACNB2', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'KCNN3', 'HCN1', 'KCNB1', 'CACNB2'] 
 #return ['CACNB2b', 'CACNA1D', 'CACNA1I', 'CACNA1S', 'ATP2A2', 'ATP2B2', 'NRGN', 'SCN1A', 'SCN9A', 'KCNS3', 'HCN1'] 

           
