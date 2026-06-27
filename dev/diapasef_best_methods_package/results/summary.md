# diaPASEF Fig.6 common-spectrum deconvolution test

Input matrices:

- Chromatogram matrix: `87 x 24`
- Mobilogram matrix: `61 x 24`

The same fragment columns are used in both matrices, so the shared spectrum matrix `A` can be estimated from
`argmin ||X - C A||_F^2 + ||Y - M A||_F^2`.

## MS-DIAL-like fixed-model result

  component assigned_site    purity S1_fraction S8_fraction
 S1_phospho    S1_phospho 0.7969819   0.7969819   0.2030181
 S8_phospho    S8_phospho 0.8100013   0.1899987   0.8100013

Top fragments:

  component rank                      fragment intensity
 S1_phospho    1  S1_phospho__obs1__mz654.3333 43664.413
 S1_phospho    2  S1_phospho__obs3__mz654.8325 29813.638
 S1_phospho    3  S8_phospho__obs7__mz780.8874 26884.114
 S1_phospho    4  S1_phospho__obs2__mz780.8820 26619.152
 S1_phospho    5  S8_phospho__obs8__mz781.3866 24204.536
 S1_phospho    6  S1_phospho__obs4__mz781.3812 24202.537
 S1_phospho    7 S1_phospho__obs5__mz1233.5555 20669.657
 S1_phospho    8  S1_phospho__obs6__mz781.8857 17264.465
 S8_phospho    1  S8_phospho__obs3__mz694.3114 17986.071
 S8_phospho    2 S8_phospho__obs10__mz788.8635 15458.326
 S8_phospho    3  S8_phospho__obs6__mz791.3811 14344.716
 S8_phospho    4  S1_phospho__obs6__mz781.8857 12769.216
 S8_phospho    5  S8_phospho__obs4__mz791.8857 11980.448
 S8_phospho    6  S8_phospho__obs5__mz694.8160 11783.953
 S8_phospho    7  S8_phospho__obs9__mz788.3589 11183.480
 S8_phospho    8  S8_phospho__obs1__mz786.7908  7081.814

## MS-DIAL-like single-peak chromatogram model result

Selected chromatogram model fragments:

       site                      fragment apex_index     score
 S1_phospho S1_phospho__obs5__mz1233.5555         64 14.993002
 S8_phospho  S8_phospho__obs3__mz694.3114         36  8.142516

  component assigned_site    purity S1_fraction S8_fraction
 S1_phospho    S1_phospho 0.7485794   0.7485794   0.2514206
 S8_phospho    S8_phospho 0.7213042   0.2786958   0.7213042

Top fragments:

  component rank                      fragment intensity
 S1_phospho    1  S1_phospho__obs1__mz654.3333 33677.689
 S1_phospho    2  S1_phospho__obs2__mz780.8820 24722.979
 S1_phospho    3  S8_phospho__obs7__mz780.8874 24682.340
 S1_phospho    4  S1_phospho__obs3__mz654.8325 22787.634
 S1_phospho    5  S1_phospho__obs4__mz781.3812 22304.214
 S1_phospho    6  S8_phospho__obs8__mz781.3866 22258.346
 S1_phospho    7  S1_phospho__obs6__mz781.8857 19292.549
 S1_phospho    8 S1_phospho__obs5__mz1233.5555 17657.468
 S8_phospho    1  S1_phospho__obs6__mz781.8857 12503.916
 S8_phospho    2  S8_phospho__obs3__mz694.3114 11403.886
 S8_phospho    3  S8_phospho__obs6__mz791.3811  8572.805
 S8_phospho    4 S8_phospho__obs10__mz788.8635  7915.525
 S8_phospho    5  S8_phospho__obs5__mz694.8160  7730.905
 S8_phospho    6  S8_phospho__obs4__mz791.8857  7339.028
 S8_phospho    7  S8_phospho__obs9__mz788.3589  7286.861
 S8_phospho    8  S8_phospho__obs1__mz786.7908  5119.906

## MS2DecR msdial_model RT-selected model result

Focused peak:

  left_rt  apex_rt right_rt
 18.20335 18.43304  18.9844

Triplet summary:

 model       mz apex_rt n_peaks
    M1       NA      NA      NA
    M2 1233.555  18.479       5
    M3       NA      NA      NA

Direct MS2DecR `msdial_deconv()` purity:

 component assigned_site   purity S1_fraction S8_fraction
        M2    S1_phospho 0.773735    0.773735    0.226265

Direct MS2DecR top fragments:

 component rank fragment intensity
        M2    1 1233.555  8202.280
        M2    2  654.333  6645.063
        M2    3  780.887  6146.554
        M2    4  780.882  6064.630
        M2    5  781.387  6021.859
        M2    6  781.381  5963.156
        M2    7  654.832  4835.210
        M2    8  781.886  2766.607

RT+IM common-A purity using the same RT models:

 component assigned_site   purity S1_fraction S8_fraction
        M2    S1_phospho 0.683422    0.683422    0.316578

Top fragments:

 component rank                      fragment intensity
        M2    1  S1_phospho__obs1__mz654.3333 13868.404
        M2    2  S1_phospho__obs2__mz780.8820 11987.961
        M2    3  S8_phospho__obs7__mz780.8874 11929.678
        M2    4  S8_phospho__obs8__mz781.3866 11106.565
        M2    5  S1_phospho__obs4__mz781.3812 11101.248
        M2    6 S1_phospho__obs5__mz1233.5555 10159.237
        M2    7  S1_phospho__obs3__mz654.8325  9399.441
        M2    8  S1_phospho__obs6__mz781.8857  8945.627

## MS2DecR msdial_model ideal_min 0.70 result

Triplet summary:

 model       mz  apex_rt n_peaks
    M1 788.3589 18.27225       1
    M2 654.3333 18.45602      13
    M3       NA       NA      NA

Direct MS2DecR `msdial_deconv()` purity:

 component assigned_site    purity S1_fraction S8_fraction
        M1    S8_phospho 0.9305168  0.06948317   0.9305168
        M2    S1_phospho 0.7723943  0.77239430   0.2276057

Direct MS2DecR top fragments:

 component rank fragment  intensity
        M1    1  788.359 2045.89035
        M1    2  788.864 1545.68198
        M1    3  654.333  310.87822
        M1    4  786.791  205.66580
        M1    5  747.860  171.97548
        M1    6  792.385  129.42923
        M1    7  694.816   78.92582
        M1    8  791.886   69.13557
        M2    1  654.333 8195.81998
        M2    2 1233.555 7420.41175
        M2    3  780.887 6832.61622
        M2    4  780.882 6765.68727
        M2    5  781.381 6143.68208
        M2    6  781.387 6094.49855
        M2    7  654.832 5343.59940
        M2    8  781.886 2823.21488

RT+IM common-A purity using the same RT models:

 component assigned_site    purity S1_fraction S8_fraction
        M1    S1_phospho 0.5405778   0.5405778   0.4594222
        M2    S1_phospho 0.7305829   0.7305829   0.2694171

Top fragments:

 component rank                      fragment intensity
        M1    1  S1_phospho__obs6__mz781.8857  7448.692
        M1    2  S1_phospho__obs1__mz654.3333  6968.770
        M1    3  S1_phospho__obs2__mz780.8820  5520.380
        M1    4  S8_phospho__obs7__mz780.8874  5385.234
        M1    5  S8_phospho__obs9__mz788.3589  4812.184
        M1    6  S1_phospho__obs4__mz781.3812  4571.293
        M1    7  S8_phospho__obs8__mz781.3866  4541.324
        M1    8  S1_phospho__obs3__mz654.8325  4271.753
        M2    1  S1_phospho__obs1__mz654.3333 14170.778
        M2    2  S8_phospho__obs7__mz780.8874 11024.569
        M2    3  S1_phospho__obs2__mz780.8820 11023.278
        M2    4  S1_phospho__obs4__mz781.3812  9840.393
        M2    5  S8_phospho__obs8__mz781.3866  9795.365
        M2    6  S1_phospho__obs3__mz654.8325  9468.663
        M2    7 S1_phospho__obs5__mz1233.5555  8324.721
        M2    8  S1_phospho__obs6__mz781.8857  7706.333

## MS2DecR RT ideal_min 0.70 + IM ideal_min 0.60 fwhm 3 result

RT triplet summary:

 model       mz  apex_rt n_peaks
    M1 788.3589 18.27225       1
    M2 654.3333 18.45602      13
    M3       NA       NA      NA

IM triplet summary:

 model       mz apex_rt n_peaks
    M1       NA      NA      NA
    M2 654.3333  0.9410       1
    M3 694.8160  0.9635       1

Model site alignment:

 model    rt_site    im_site
    M1 S8_phospho S8_phospho
    M2 S1_phospho S1_phospho

RT+IM common-A purity:

 component assigned_site    purity S1_fraction S8_fraction
        M1    S1_phospho 0.5145569   0.5145569   0.4854431
        M2    S1_phospho 0.7195227   0.7195227   0.2804773

Top fragments:

 component rank                      fragment intensity
        M1    1  S1_phospho__obs1__mz654.3333  3298.858
        M1    2  S1_phospho__obs2__mz780.8820  3112.037
        M1    3  S8_phospho__obs7__mz780.8874  3024.699
        M1    4  S1_phospho__obs4__mz781.3812  2829.395
        M1    5  S8_phospho__obs8__mz781.3866  2802.810
        M1    6  S1_phospho__obs6__mz781.8857  2521.515
        M1    7  S8_phospho__obs9__mz788.3589  2292.197
        M1    8 S8_phospho__obs10__mz788.8635  2261.620
        M2    1  S1_phospho__obs1__mz654.3333 10630.245
        M2    2  S8_phospho__obs7__mz780.8874  7714.774
        M2    3  S1_phospho__obs2__mz780.8820  7688.388
        M2    4  S1_phospho__obs3__mz654.8325  6960.298
        M2    5  S1_phospho__obs4__mz781.3812  6671.810
        M2    6  S8_phospho__obs8__mz781.3866  6653.817
        M2    7  S1_phospho__obs6__mz781.8857  6061.098
        M2    8 S1_phospho__obs5__mz1233.5555  4454.682

RT+IM window + baseline common-A purity:

 component assigned_site    purity S1_fraction S8_fraction
        M1    S8_phospho 0.9684297  0.03157025   0.9684297
        M2    S1_phospho 0.7705776  0.77057765   0.2294224

Top fragments:

 component rank                      fragment intensity
        M1    1 S8_phospho__obs10__mz788.8635 1449.8413
        M1    2  S8_phospho__obs9__mz788.3589 1289.7885
        M1    3  S8_phospho__obs3__mz694.3114  797.4860
        M1    4  S8_phospho__obs5__mz694.8160  766.0931
        M1    5 S8_phospho__obs12__mz792.3849  524.1222
        M1    6  S8_phospho__obs4__mz791.8857  461.1601
        M1    7  S8_phospho__obs6__mz791.3811  228.0517
        M1    8 S8_phospho__obs11__mz747.8598  190.3064
        M2    1  S1_phospho__obs1__mz654.3333 6666.7991
        M2    2  S8_phospho__obs7__mz780.8874 3477.5534
        M2    3  S1_phospho__obs2__mz780.8820 3379.4733
        M2    4  S1_phospho__obs3__mz654.8325 2938.6434
        M2    5  S1_phospho__obs6__mz781.8857 2594.7981
        M2    6 S1_phospho__obs5__mz1233.5555 1807.2011
        M2    7  S1_phospho__obs4__mz781.3812 1690.9972
        M2    8  S8_phospho__obs8__mz781.3866 1503.7779

## Joint MCR-ALS result

 component assigned_site    purity S1_fraction S8_fraction
      ALS1    S8_phospho 0.6327940   0.3672060   0.6327940
      ALS2    S1_phospho 0.7447384   0.7447384   0.2552616

Top fragments:

 component rank                      fragment intensity
      ALS1    1  S1_phospho__obs6__mz781.8857  28076.08
      ALS1    2 S8_phospho__obs10__mz788.8635  15157.31
      ALS1    3  S8_phospho__obs3__mz694.3114  15041.08
      ALS1    4  S8_phospho__obs6__mz791.3811  12571.30
      ALS1    5  S8_phospho__obs9__mz788.3589  11433.77
      ALS1    6  S8_phospho__obs4__mz791.8857  10759.28
      ALS1    7  S8_phospho__obs5__mz694.8160  10124.53
      ALS1    8  S1_phospho__obs2__mz780.8820   9835.50
      ALS2    1  S1_phospho__obs1__mz654.3333  49791.35
      ALS2    2  S1_phospho__obs2__mz780.8820  35448.41
      ALS2    3  S8_phospho__obs7__mz780.8874  35364.44
      ALS2    4  S1_phospho__obs3__mz654.8325  33343.70
      ALS2    5  S1_phospho__obs4__mz781.3812  32108.32
      ALS2    6  S8_phospho__obs8__mz781.3866  31995.39
      ALS2    7 S1_phospho__obs5__mz1233.5555  22893.45
      ALS2    8  S1_phospho__obs6__mz781.8857  22730.22

## ALS error

- initial: 58191.3457
- final: 33874.1774

## Model-initialized joint MCR-ALS result

 component assigned_site    purity S1_fraction S8_fraction
      ALS1    S1_phospho 0.7447384   0.7447384   0.2552616
      ALS2    S8_phospho 0.6327940   0.3672060   0.6327940

Top fragments:

 component rank                      fragment intensity
      ALS1    1  S1_phospho__obs1__mz654.3333  49791.35
      ALS1    2  S1_phospho__obs2__mz780.8820  35448.41
      ALS1    3  S8_phospho__obs7__mz780.8874  35364.44
      ALS1    4  S1_phospho__obs3__mz654.8325  33343.70
      ALS1    5  S1_phospho__obs4__mz781.3812  32108.32
      ALS1    6  S8_phospho__obs8__mz781.3866  31995.39
      ALS1    7 S1_phospho__obs5__mz1233.5555  22893.45
      ALS1    8  S1_phospho__obs6__mz781.8857  22730.22
      ALS2    1  S1_phospho__obs6__mz781.8857  28076.08
      ALS2    2 S8_phospho__obs10__mz788.8635  15157.31
      ALS2    3  S8_phospho__obs3__mz694.3114  15041.08
      ALS2    4  S8_phospho__obs6__mz791.3811  12571.30
      ALS2    5  S8_phospho__obs9__mz788.3589  11433.77
      ALS2    6  S8_phospho__obs4__mz791.8857  10759.28
      ALS2    7  S8_phospho__obs5__mz694.8160  10124.53
      ALS2    8  S1_phospho__obs2__mz780.8820   9835.50

## Model-initialized ALS error

- initial: 46949.6282
- final: 33874.1774

## Overcomplete peak-selected MCR-ALS result

Component selection table:

 component idx c_pass m_pass c_apex m_apex        score  spec_sum
      ALS1   1   TRUE  FALSE     36     45 1.441943e-06 107171.88
      ALS2   2   TRUE   TRUE     64     20 7.014532e-07 271077.71
      ALS3   3   TRUE  FALSE     80     47 9.883903e-08  47214.99
      ALS4   4   TRUE  FALSE     64     50 3.577519e-07  78691.07
      ALS5   5  FALSE   TRUE      7     35 1.124239e-07  67296.35

Selected component purity:

 component assigned_site    purity S1_fraction S8_fraction
      ALS1    S8_phospho 0.7728966   0.2271034   0.7728966
      ALS2    S1_phospho 0.7177634   0.7177634   0.2822366

Top fragments:

 component rank                     fragment intensity
      ALS1    1 S8_phospho__obs3__mz694.3114 14737.763
      ALS1    2 S8_phospho__obs6__mz791.3811 10482.628
      ALS1    3 S8_phospho__obs5__mz694.8160  9952.497
      ALS1    4 S8_phospho__obs4__mz791.8857  8743.758
      ALS1    5 S8_phospho__obs9__mz788.3589  7394.311
      ALS1    6 S8_phospho__obs1__mz786.7908  6321.994
      ALS1    7 S1_phospho__obs2__mz780.8820  6212.632
      ALS1    8 S8_phospho__obs7__mz780.8874  6092.559
      ALS2    1 S1_phospho__obs1__mz654.3333 51160.315
      ALS2    2 S1_phospho__obs3__mz654.8325 33571.284
      ALS2    3 S8_phospho__obs7__mz780.8874 30434.258
      ALS2    4 S1_phospho__obs2__mz780.8820 30414.782
      ALS2    5 S1_phospho__obs4__mz781.3812 26034.243
      ALS2    6 S8_phospho__obs8__mz781.3866 25787.256
      ALS2    7 S1_phospho__obs6__mz781.8857 16574.801
      ALS2    8 S1_phospho__obs8__mz655.3370 14380.899

## MCR-ALS post-selection strategy comparison

                                method S1_best_purity S8_best_purity
                    A_peak_select_only      0.7177634      0.7728966
     B_fixed_CM_window_baseline_A_only      0.7421356      0.7690327
             C_reALS_after_peak_select      0.7447384      0.6327940
 MSDIAL_like_window_baseline_reference      0.7705776      0.9684297
 mean_best_purity    error
        0.7453300 18627.09
        0.7555842       NA
        0.6887662 33874.18
        0.8695037       NA

### B. Fixed C/M, A-only window-baseline refit

 component assigned_site    purity S1_fraction S8_fraction
      ALS1    S8_phospho 0.7690327   0.2309673   0.7690327
      ALS2    S1_phospho 0.7421356   0.7421356   0.2578644

Top fragments:

 component rank                      fragment intensity
      ALS1    1  S8_phospho__obs3__mz694.3114 15321.120
      ALS1    2  S8_phospho__obs5__mz694.8160 10632.104
      ALS1    3  S8_phospho__obs6__mz791.3811  8923.276
      ALS1    4  S8_phospho__obs4__mz791.8857  7736.073
      ALS1    5  S1_phospho__obs2__mz780.8820  6814.568
      ALS1    6  S8_phospho__obs7__mz780.8874  6729.257
      ALS1    7  S8_phospho__obs8__mz781.3866  5930.125
      ALS1    8  S1_phospho__obs4__mz781.3812  5879.821
      ALS2    1  S1_phospho__obs1__mz654.3333 50388.092
      ALS2    2  S8_phospho__obs7__mz780.8874 37695.350
      ALS2    3  S1_phospho__obs2__mz780.8820 37225.701
      ALS2    4  S1_phospho__obs3__mz654.8325 34053.637
      ALS2    5  S8_phospho__obs8__mz781.3866 32448.001
      ALS2    6  S1_phospho__obs4__mz781.3812 32374.380
      ALS2    7  S1_phospho__obs6__mz781.8857 28344.440
      ALS2    8 S1_phospho__obs5__mz1233.5555 22912.868

### C. Re-ALS after peak selection

 component assigned_site    purity S1_fraction S8_fraction
      ALS1    S8_phospho 0.6327940   0.3672060   0.6327940
      ALS2    S1_phospho 0.7447384   0.7447384   0.2552616

Top fragments:

 component rank                      fragment intensity
      ALS1    1  S1_phospho__obs6__mz781.8857 28076.083
      ALS1    2 S8_phospho__obs10__mz788.8635 15157.308
      ALS1    3  S8_phospho__obs3__mz694.3114 15041.077
      ALS1    4  S8_phospho__obs6__mz791.3811 12571.304
      ALS1    5  S8_phospho__obs9__mz788.3589 11433.769
      ALS1    6  S8_phospho__obs4__mz791.8857 10759.280
      ALS1    7  S8_phospho__obs5__mz694.8160 10124.531
      ALS1    8  S1_phospho__obs2__mz780.8820  9835.499
      ALS2    1  S1_phospho__obs1__mz654.3333 49791.350
      ALS2    2  S1_phospho__obs2__mz780.8820 35448.409
      ALS2    3  S8_phospho__obs7__mz780.8874 35364.440
      ALS2    4  S1_phospho__obs3__mz654.8325 33343.696
      ALS2    5  S1_phospho__obs4__mz781.3812 32108.318
      ALS2    6  S8_phospho__obs8__mz781.3866 31995.395
      ALS2    7 S1_phospho__obs5__mz1233.5555 22893.452
      ALS2    8  S1_phospho__obs6__mz781.8857 22730.221

## Overcomplete ALS error

- initial: 52214.5740
- final: 18627.0872

## Summary Images

- `msdial_like_deconvolution_summary.png`
- `msdial_like_single_peak_chrom_deconvolution_summary.png`
- `ms2decr_msdial_model_deconvolution_summary.png`
- `ms2decr_msdial_ideal070_deconvolution_summary.png`
- `ms2decr_rt070_im060_deconvolution_summary.png`
- `ms2decr_rt070_im060_window_baseline_deconvolution_summary.png`
- `mcr_als_deconvolution_summary.png`
- `mcr_als_model_init_deconvolution_summary.png`
- `mcr_als_overcomplete_peak_selected_deconvolution_summary.png`
- `mcr_als_overcomplete_peak_selected_Aonly_deconvolution_summary.png`
- `mcr_als_overcomplete_peak_selected_reals_deconvolution_summary.png`
