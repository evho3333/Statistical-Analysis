
#To-Do:
    #Outlier Analysis??? Generalized ESD test
    #Tukey-Kramer p values??? Probably have to use Excel
    #Add IL-10:IL-12p40 ratio analysis

#On-Hold due to matplotlib compatibility and crashing issues:
    #Add file GUI option instead of typing file directory

import sys
import pandas

def LPS_PCR_posthoc_test (N, BBS, ten, thirty, hundred, threehundred, Mt, BBS_LPS, ten_LPS, thirty_LPS,
                          hundred_LPS, threehundred_LPS, Mt_LPS, BBS_LPS2, ten_LPS2, thirty_LPS2,
                          hundred_LPS2, threehundred_LPS2, Mt_LPS2, BBS_LPS3, ten_LPS3, thirty_LPS3,
                          hundred_LPS3, threehundred_LPS3, Mt_LPS3):
    import statistics
    import math
    import scipy.stats
    print("")
    print("Running - Function - LPS_PCR_posthoc_test...")

##    Listt = []
##
##    for a in N:
##        if a == ",":
##            continue
##        Listt.append (float(a))
##
##    N = sum (Listt)
    
    print ("")
    print ("1: Unprotected Fisher's LSD: doesn't correct for multiple comparisons")
    print ("")
    print ("2: Tukey's HSD: corrects for multiple comparisons and assumes every treatment has an equal n (no NaN values)")
    print ("")
    print ("3: Tukey-Kramer Method: corrects for multiple comparisons, but doesn't assume equal n for each treatment")
    print ("")
    print ("4: Dunnnett's Test: compares each bacterial dose to the control BBS group; assumes equal n for each treatment in critical value")
    print ("")
    answer = input ("Which Post-Hoc Test Would You Like to Run? (Enter The Number): ")

    #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/New T Critical Values.csv", header = None)

    crit_table = [12.7065, 4.3026, 3.1824, 2.7764, 2.5706, 2.4469, 2.3646, 2.306, 2.2621, 2.2282, 2.201, 2.1788, 2.1604,	2.1448,	2.1314,	2.1199,	2.1098,	2.1009,	2.093, 2.086, 2.0796, 2.0739, 2.0686, 2.0639, 2.0596, 2.0555, 2.0518, 2.0484, 2.0452, 2.0423, 2.0395, 2.0369, 2.0345, 2.0322, 2.0301, 2.0281, 2.0262, 2.0244, 2.0227, 2.0211, 2.0196, 2.0181, 2.0167, 2.0154, 2.0141, 2.0129, 2.0117, 2.0106, 2.0096, 2.0086, 2.0076, 2.0066, 2.0057, 2.0049, 2.0041, 2.0032, 2.0025, 2.0017, 2.001, 2.0003, 1.9996, 1.999, 1.9983, 1.9977, 1.9971, 1.9966, 1.996, 1.9955, 1.995, 1.9944, 1.9939, 1.9935, 1.993,1.9925, 1.9921, 1.9917, 1.9913, 1.9909, 1.9904, 1.9901, 1.9897, 1.9893, 1.9889, 1.9886, 1.9883, 1.9879, 1.9876, 1.9873, 1.987, 1.9867, 1.9864, 1.9861, 1.9858, 1.9855, 1.9852, 1.985, 1.9847, 1.9845, 1.9842, 1.984, 1.9837, 1.9835, 1.9833, 1.983, 1.9828, 1.9826, 1.9824, 1.9822, 1.982, 1.9818, 1.9816, 1.9814, 1.9812, 1.981, 1.9808, 1.9806, 1.9805, 1.9803, 1.9801, 1.9799]
    
    #################################################
  
    df_BBS = len (BBS) - 1
    df_ten = len (ten) - 1
    df_thirty = len (thirty) - 1
    df_hundred = len (hundred) - 1
    df_threehundred = len (threehundred) - 1
    df_Mt = len (Mt) - 1

    df_BBS_LPS = len (BBS_LPS) - 1
    df_ten_LPS = len (ten_LPS) - 1
    df_thirty_LPS = len (thirty_LPS) - 1
    df_hundred_LPS = len (hundred_LPS) - 1
    df_threehundred_LPS = len (threehundred_LPS) - 1
    df_Mt_LPS = len (Mt_LPS) - 1

    df_BBS_LPS2 = len (BBS_LPS2) - 1
    df_ten_LPS2 = len (ten_LPS2) - 1
    df_thirty_LPS2 = len (thirty_LPS2) - 1
    df_hundred_LPS2 = len (hundred_LPS2) - 1
    df_threehundred_LPS2 = len (threehundred_LPS2) - 1
    df_Mt_LPS2 = len (Mt_LPS2) - 1

    df_BBS_LPS3 = len (BBS_LPS3) - 1
    df_ten_LPS3 = len (ten_LPS3) - 1
    df_thirty_LPS3 = len (thirty_LPS3) - 1
    df_hundred_LPS3 = len (hundred_LPS3) - 1
    df_threehundred_LPS3 = len (threehundred_LPS3) - 1
    df_Mt_LPS3 = len (Mt_LPS3) - 1
    
    #################################################
    
    if len (BBS) == 1:
        BBS_var = 0
    else:
        BBS_var = statistics.variance (BBS)

    if len(ten) == 1:
        ten_var = 0
    else:
        ten_var = statistics.variance (ten)

    if len(thirty) == 1:
        thirty_var = 0
    else:
        thirty_var = statistics.variance (thirty)

    if len(hundred) == 1:
        hundred_var = 0
    else:
        hundred_var = statistics.variance (hundred)

    if len (threehundred) == 1:
        threehundred_var = 0
    else:
        threehundred_var = statistics.variance (threehundred)

    if len (Mt) == 1:
        Mt_var = 0
    else:
        Mt_var = statistics.variance (Mt)

    ############################

    if len (BBS_LPS) == 1:
        BBS_LPS_var = 0
    else:
        BBS_LPS_var = statistics.variance (BBS_LPS)

    if len(ten_LPS) == 1:
        ten_LPS_var = 0
    else:
        ten_LPS_var = statistics.variance (ten_LPS)

    if len(thirty_LPS) == 1:
        thirty_LPS_var = 0
    else:
        thirty_LPS_var = statistics.variance (thirty_LPS)

    if len(hundred_LPS) == 1:
        hundred_LPS_var = 0
    else:
        hundred_LPS_var = statistics.variance (hundred_LPS)

    if len (threehundred_LPS) == 1:
        threehundred_LPS_var = 0
    else:
        threehundred_LPS_var = statistics.variance (threehundred_LPS)

    if len (Mt_LPS) == 1:
        Mt_LPS_var = 0
    else:
        Mt_LPS_var = statistics.variance (Mt_LPS)

    #################################################

    if len (BBS_LPS2) == 1:
        BBS_LPS2_var = 0
    else:
        BBS_LPS2_var = statistics.variance (BBS_LPS2)

    if len(ten_LPS2) == 1:
        ten_LPS2_var = 0
    else:
        ten_LPS2_var = statistics.variance (ten_LPS2)

    if len(thirty_LPS2) == 1:
        thirty_LPS2_var = 0
    else:
        thirty_LPS2_var = statistics.variance (thirty_LPS2)

    if len(hundred_LPS2) == 1:
        hundred_LPS2_var = 0
    else:
        hundred_LPS2_var = statistics.variance (hundred_LPS2)

    if len (threehundred_LPS2) == 1:
        threehundred_LPS2_var = 0
    else:
        threehundred_LPS2_var = statistics.variance (threehundred_LPS2)

    if len (Mt_LPS2) == 1:
        Mt_LPS2_var = 0
    else:
        Mt_LPS2_var = statistics.variance (Mt_LPS2)

    #################################################

    if len (BBS_LPS3) == 1:
        BBS_LPS3_var = 0
    else:
        BBS_LPS3_var = statistics.variance (BBS_LPS3)

    if len(ten_LPS3) == 1:
        ten_LPS3_var = 0
    else:
        ten_LPS3_var = statistics.variance (ten_LPS3)

    if len(thirty_LPS3) == 1:
        thirty_LPS3_var = 0
    else:
        thirty_LPS3_var = statistics.variance (thirty_LPS3)

    if len(hundred_LPS3) == 1:
        hundred_LPS3_var = 0
    else:
        hundred_LPS3_var = statistics.variance (hundred_LPS3)

    if len (threehundred_LPS3) == 1:
        threehundred_LPS3_var = 0
    else:
        threehundred_LPS3_var = statistics.variance (threehundred_LPS3)

    if len (Mt_LPS3) == 1:
        Mt_LPS3_var = 0
    else:
        Mt_LPS3_var = statistics.variance (Mt_LPS3)

    #################################################

    SS_within_control = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var)
    df_within_control = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred
    df_within_control_actual = int (float(df_within_control) - 1)
    MSE_within_control = SS_within_control / df_within_control
    
    SS_within_LPS = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var)
    df_within_LPS = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS
    df_within_LPS_actual = int (float (df_within_LPS) - 1)
    MSE_within_LPS = SS_within_LPS / df_within_LPS
    
    SS_within_LPS2 = (df_BBS_LPS2 * BBS_LPS2_var) + (df_ten_LPS2 * ten_LPS2_var) + (df_thirty_LPS2 * thirty_LPS2_var) + (df_hundred_LPS2 * hundred_LPS2_var) + (df_threehundred_LPS2 * threehundred_LPS2_var)
    df_within_LPS2 = df_BBS_LPS2 + df_ten_LPS2 + df_thirty_LPS2 + df_hundred_LPS2 + df_threehundred_LPS2
    df_within_LPS2_actual = int (float (df_within_LPS2) - 1)
    MSE_within_LPS2 = SS_within_LPS2 / df_within_LPS2
    
    SS_within_LPS3 = (df_BBS_LPS3 * BBS_LPS3_var) + (df_ten_LPS3 * ten_LPS3_var) + (df_thirty_LPS3 * thirty_LPS3_var) + (df_hundred_LPS3 * hundred_LPS3_var) + (df_threehundred_LPS3 * threehundred_LPS3_var)
    df_within_LPS3 = df_BBS_LPS3 + df_ten_LPS3 + df_thirty_LPS3 + df_hundred_LPS3 + df_threehundred_LPS3
    df_within_LPS3_actual = int (float (df_within_LPS3) - 1)
    MSE_within_LPS3 = SS_within_LPS3 / df_within_LPS3

    SS_within_LPSvLPS2 = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var) + (df_BBS_LPS2 * BBS_LPS2_var) + (df_ten_LPS2 * ten_LPS2_var) + (df_thirty_LPS2 * thirty_LPS2_var) + (df_hundred_LPS2 * hundred_LPS2_var) + (df_threehundred_LPS2 * threehundred_LPS2_var)
    df_within_LPSvLPS2 = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS + df_BBS_LPS2 + df_ten_LPS2 + df_thirty_LPS2 + df_hundred_LPS2 + df_threehundred_LPS2
    df_within_LPSvLPS2_actual = int (float (df_within_LPSvLPS2) - 1)
    MSE_within_LPSvLPS2 = SS_within_LPSvLPS2 / df_within_LPSvLPS2

    SS_within_LPSvLPS3 = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var) + (df_BBS_LPS3 * BBS_LPS3_var) + (df_ten_LPS3 * ten_LPS3_var) + (df_thirty_LPS3 * thirty_LPS3_var) + (df_hundred_LPS3 * hundred_LPS3_var) + (df_threehundred_LPS3 * threehundred_LPS3_var)
    df_within_LPSvLPS3 = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS + df_BBS_LPS3 + df_ten_LPS3 + df_thirty_LPS3 + df_hundred_LPS3 + df_threehundred_LPS3
    df_within_LPSvLPS3_actual = int (float (df_within_LPSvLPS3) - 1)
    MSE_within_LPSvLPS3 = SS_within_LPSvLPS3 / df_within_LPSvLPS3

    SS_within_LPS2vLPS3 = (df_BBS_LPS2 * BBS_LPS2_var) + (df_ten_LPS2 * ten_LPS2_var) + (df_thirty_LPS2 * thirty_LPS2_var) + (df_hundred_LPS2 * hundred_LPS2_var) + (df_threehundred_LPS2 * threehundred_LPS2_var) + (df_BBS_LPS3 * BBS_LPS3_var) + (df_ten_LPS3 * ten_LPS3_var) + (df_thirty_LPS3 * thirty_LPS3_var) + (df_hundred_LPS3 * hundred_LPS3_var) + (df_threehundred_LPS3 * threehundred_LPS3_var)
    df_within_LPS2vLPS3 = df_BBS_LPS2 + df_ten_LPS2 + df_thirty_LPS2 + df_hundred_LPS2 + df_threehundred_LPS2 + df_BBS_LPS3 + df_ten_LPS3 + df_thirty_LPS3 + df_hundred_LPS3 + df_threehundred_LPS3
    df_within_LPS2vLPS3_actual = int (float (df_within_LPS2vLPS3) - 1)
    MSE_within_LPS2vLPS3 = SS_within_LPS2vLPS3 / df_within_LPS2vLPS3

    SS_within_CvLPS = SS_within_control + SS_within_LPS
    df_within_CvLPS = df_within_control + df_within_LPS
    df_within_CvLPS_actual = int (float(df_within_CvLPS)-1)
    MSE_within_CvLPS = SS_within_CvLPS / df_within_CvLPS

    SS_within_CvLPS2 = SS_within_control + SS_within_LPS2
    df_within_CvLPS2 = df_within_control + df_within_LPS2
    df_within_CvLPS2_actual = int (float(df_within_CvLPS2)-1)
    MSE_within_CvLPS2 = SS_within_CvLPS2 / df_within_CvLPS2

    SS_within_CvLPS3 = SS_within_control + SS_within_LPS3
    df_within_CvLPS3 = df_within_control + df_within_LPS3
    df_within_CvLPS3_actual = int (float(df_within_CvLPS3)-1)
    MSE_within_CvLPS3 = SS_within_CvLPS3 / df_within_CvLPS3

    SS_within_BBS_by_LPS = (df_BBS*BBS_var)+(df_BBS_LPS*BBS_LPS_var)+(df_BBS_LPS2*BBS_LPS2_var)+(df_BBS_LPS3*BBS_LPS3_var)
    df_within_BBS_by_LPS = df_BBS+df_BBS_LPS+df_BBS_LPS2+df_BBS_LPS3
    df_within_BBS_by_LPS_actual = int (float (df_within_BBS_by_LPS) - 1)
    MSE_within_BBS_by_LPS = SS_within_BBS_by_LPS / df_within_BBS_by_LPS

    SS_within_ten_by_LPS = (df_ten*ten_var)+(df_ten_LPS*ten_LPS_var)+(df_ten_LPS2*ten_LPS2_var)+(df_ten_LPS3*ten_LPS3_var)
    df_within_ten_by_LPS = df_ten+df_ten_LPS+df_ten_LPS2+df_ten_LPS3
    df_within_ten_by_LPS_actual = int (float (df_within_ten_by_LPS) - 1)
    MSE_within_ten_by_LPS = SS_within_ten_by_LPS / df_within_ten_by_LPS

    SS_within_thirty_by_LPS = (df_thirty*thirty_var)+(df_thirty_LPS*thirty_LPS_var)+(df_thirty_LPS2*thirty_LPS2_var)+(df_thirty_LPS3*thirty_LPS3_var)
    df_within_thirty_by_LPS = df_thirty+df_thirty_LPS+df_thirty_LPS2+df_thirty_LPS3
    df_within_thirty_by_LPS_actual = int (float (df_within_thirty_by_LPS) - 1)
    MSE_within_thirty_by_LPS = SS_within_thirty_by_LPS / df_within_thirty_by_LPS

    SS_within_hundred_by_LPS = (df_hundred*hundred_var)+(df_hundred_LPS*hundred_LPS_var)+(df_hundred_LPS2*hundred_LPS2_var)+(df_hundred_LPS3*hundred_LPS3_var)
    df_within_hundred_by_LPS = df_hundred+df_hundred_LPS+df_hundred_LPS2+df_hundred_LPS3
    df_within_hundred_by_LPS_actual = int (float (df_within_hundred_by_LPS) - 1)
    MSE_within_hundred_by_LPS = SS_within_hundred_by_LPS / df_within_hundred_by_LPS

    SS_within_threehundred_by_LPS = (df_threehundred*threehundred_var)+(df_threehundred_LPS*threehundred_LPS_var)+(df_threehundred_LPS2*threehundred_LPS2_var)+(df_threehundred_LPS3*threehundred_LPS3_var)
    df_within_threehundred_by_LPS = df_threehundred+df_threehundred_LPS+df_threehundred_LPS2+df_threehundred_LPS3
    df_within_threehundred_by_LPS_actual = int (float (df_within_threehundred_by_LPS) - 1)
    MSE_within_threehundred_by_LPS = SS_within_threehundred_by_LPS / df_within_threehundred_by_LPS

    SS_within_Mt_by_LPS = (df_Mt*Mt_var)+(df_Mt_LPS*Mt_LPS_var)+(df_Mt_LPS2*Mt_LPS2_var)+(df_Mt_LPS3*Mt_LPS3_var)
    df_within_Mt_by_LPS = df_Mt + df_Mt_LPS + df_Mt_LPS2 + df_Mt_LPS3
    df_within_Mt_by_LPS_actual = int(float(df_within_Mt_by_LPS)-1)
    MSE_within_Mt_by_LPS = SS_within_Mt_by_LPS / df_within_Mt_by_LPS



    Pooled_SS_within_all_LPS_for_vaccae = SS_within_control + SS_within_LPS + SS_within_LPS2 + SS_within_LPS3
    Pooled_df_within_all_LPS_for_vaccae = df_within_control + df_within_LPS + df_within_LPS2 + df_within_LPS3
    Pooled_df_within_all_LPS_for_vaccae_actual = None
    Pooled_MSE_within_all_LPS_for_vaccae = Pooled_SS_within_all_LPS_for_vaccae / Pooled_df_within_all_LPS_for_vaccae
    print("MSE for vaccae:", Pooled_MSE_within_all_LPS_for_vaccae) #verified: MSE for vaccae and MSE for LPS are the same value

    Pooled_SS_within_all_LPS_for_LPS = SS_within_BBS_by_LPS + SS_within_ten_by_LPS + SS_within_thirty_by_LPS + SS_within_hundred_by_LPS + SS_within_threehundred_by_LPS
    Pooled_df_within_all_LPS_for_LPS = df_within_BBS_by_LPS + df_within_ten_by_LPS + df_within_thirty_by_LPS + df_within_hundred_by_LPS + df_within_threehundred_by_LPS
    Pooled_df_within_all_LPS_for_LPS_actual = None
    Pooled_MSE_within_all_LPS_for_LPS = Pooled_SS_within_all_LPS_for_LPS / Pooled_df_within_all_LPS_for_LPS
    print("MSE for LPS:", Pooled_MSE_within_all_LPS_for_LPS) #verified: MSE for vaccae and MSE for LPS are the same value

    Pooled_SS_within_all_LPS_for_LPS_Mt = SS_within_BBS_by_LPS + SS_within_Mt_by_LPS
    Pooled_df_within_all_LPS_for_LPS_Mt = df_within_BBS_by_LPS + df_within_Mt_by_LPS
    Pooled_df_within_all_LPS_for_LPS_Mt_actual = None
    Pooled_MSE_within_all_LPS_for_LPS_Mt = Pooled_SS_within_all_LPS_for_LPS_Mt / Pooled_df_within_all_LPS_for_LPS_Mt
    print("MSE for TB LPS:", Pooled_MSE_within_all_LPS_for_LPS_Mt)
    
    #################################################
    
    if answer == "1":
    
        #Fisher LSD: any difference larger than the LSD is significant

        #t_crit_control = crit_table.iloc[df_within_control_actual, 0]
        #t_crit_LPS = crit_table.iloc[df_within_LPS_actual, 0]

        t_crit_control = crit_table[df_within_control_actual]
        t_crit_LPS = crit_table[df_within_LPS_actual]
        t_crit_LPS2 = crit_table[df_within_LPS2_actual]
        t_crit_LPS3 = crit_table[df_within_LPS3_actual]
        t_crit_LPSvLPS2 = crit_table[df_within_LPSvLPS2_actual]
        t_crit_LPSvLPS3 = crit_table[df_within_LPSvLPS3_actual]
        t_crit_LPS2vLPS3 = crit_table[df_within_LPS2vLPS3_actual]
        t_crit_CvLPS = crit_table[df_within_CvLPS_actual]
        t_crit_CvLPS2 = crit_table[df_within_CvLPS2_actual]
        t_crit_CvLPS3 = crit_table[df_within_CvLPS3_actual]

        #################################################
     
        LSD_BBS_ten = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        LSD_BBS_thirty = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        LSD_BBS_hundred = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        LSD_BBS_threehundred = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))
        LSD_BBS_Mt = abs (((sum (BBS) / len (BBS)) - (sum (Mt) / len (Mt))))

        LSD_BBS_ten_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        LSD_BBS_thirty_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        LSD_BBS_hundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        LSD_BBS_threehundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))
        LSD_BBS_Mt_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (Mt_LPS) / len (Mt_LPS))))

        LSD_BBS_ten_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (ten_LPS2) / len (ten_LPS2))))
        LSD_BBS_thirty_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        LSD_BBS_hundred_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        LSD_BBS_threehundred_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))
        LSD_BBS_Mt_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (Mt_LPS2) / len (Mt_LPS2))))

        LSD_BBS_ten_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (ten_LPS3) / len (ten_LPS3))))
        LSD_BBS_thirty_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        LSD_BBS_hundred_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        LSD_BBS_threehundred_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))
        LSD_BBS_Mt_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (Mt_LPS3) / len (Mt_LPS3))))

        LSD_BBS_LPSvLPS2 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        LSD_ten_LPSvLPS2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS2) / len (ten_LPS2))))
        LSD_thirty_LPSvLPS2 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        LSD_hundred_LPSvLPS2 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        LSD_threehundred_LPSvLPS2 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        LSD_BBS_LPSvLPS3 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        LSD_ten_LPSvLPS3 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS3) / len (ten_LPS3))))
        LSD_thirty_LPSvLPS3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        LSD_hundred_LPSvLPS3 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        LSD_threehundred_LPSvLPS3 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        LSD_BBS_LPS2vLPS3 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        LSD_ten_LPS2vLPS3 = abs (((sum (ten_LPS2) / len (ten_LPS2)) - (sum (ten_LPS3) / len (ten_LPS3))))
        LSD_thirty_LPS2vLPS3 = abs (((sum (thirty_LPS2) / len (thirty_LPS2)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        LSD_hundred_LPS2vLPS3 = abs (((sum (hundred_LPS2) / len (hundred_LPS2)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        LSD_threehundred_LPS2vLPS3 = abs (((sum (threehundred_LPS2) / len (threehundred_LPS2)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        LSD_BBS_CvLPS = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS) / len (BBS_LPS))))
        LSD_ten_CvLPS = abs (((sum (ten) / len (ten)) - (sum (ten_LPS) / len (ten_LPS))))
        LSD_thirty_CvLPS = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS) / len (thirty_LPS))))
        LSD_hundred_CvLPS = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS) / len (hundred_LPS))))
        LSD_threehundred_CvLPS = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        LSD_BBS_CvLPS2 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        LSD_ten_CvLPS2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS2) / len (ten_LPS2))))
        LSD_thirty_CvLPS2 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        LSD_hundred_CvLPS2 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        LSD_threehundred_CvLPS2 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        LSD_BBS_CvLPS3 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        LSD_ten_CvLPS3 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS3) / len (ten_LPS3))))
        LSD_thirty_CvLPS3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        LSD_hundred_CvLPS3 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        LSD_threehundred_CvLPS3 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        #################################################
        
        LSD_control_1 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(ten))))
        LSD_control_2 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(thirty))))
        LSD_control_3 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(hundred))))
        LSD_control_4 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(threehundred))))
        LSD_control_5 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(Mt))))

        LSD_LPS_1 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(ten_LPS))))
        LSD_LPS_2 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(thirty_LPS))))
        LSD_LPS_3 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(hundred_LPS))))
        LSD_LPS_4 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(threehundred_LPS))))
        LSD_LPS_5 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(Mt_LPS))))

        LSD_LPS2_1 = t_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ((1/len(BBS_LPS2)) + (1/len(ten_LPS2))))
        LSD_LPS2_2 = t_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ((1/len(BBS_LPS2)) + (1/len(thirty_LPS2))))
        LSD_LPS2_3 = t_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ((1/len(BBS_LPS2)) + (1/len(hundred_LPS2))))
        LSD_LPS2_4 = t_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ((1/len(BBS_LPS2)) + (1/len(threehundred_LPS2))))
        LSD_LPS2_5 = t_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ((1/len(BBS_LPS2)) + (1/len(Mt_LPS2))))

        LSD_LPS3_1 = t_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ((1/len(BBS_LPS3)) + (1/len(ten_LPS3))))
        LSD_LPS3_2 = t_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ((1/len(BBS_LPS3)) + (1/len(thirty_LPS3))))
        LSD_LPS3_3 = t_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ((1/len(BBS_LPS3)) + (1/len(hundred_LPS3))))
        LSD_LPS3_4 = t_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ((1/len(BBS_LPS3)) + (1/len(threehundred_LPS3))))
        LSD_LPS3_5 = t_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ((1/len(BBS_LPS3)) + (1/len(Mt_LPS3))))

        LSD_LPSvLPS2_1 = t_crit_LPSvLPS2 * math.sqrt (MSE_within_LPSvLPS2 * ((1/len(BBS_LPS)) + (1/len(BBS_LPS2))))
        LSD_LPSvLPS2_2 = t_crit_LPSvLPS2 * math.sqrt (MSE_within_LPSvLPS2 * ((1/len(ten_LPS)) + (1/len(ten_LPS2))))
        LSD_LPSvLPS2_3 = t_crit_LPSvLPS2 * math.sqrt (MSE_within_LPSvLPS2 * ((1/len(thirty_LPS)) + (1/len(thirty_LPS2))))
        LSD_LPSvLPS2_4 = t_crit_LPSvLPS2 * math.sqrt (MSE_within_LPSvLPS2 * ((1/len(hundred_LPS)) + (1/len(hundred_LPS2))))
        LSD_LPSvLPS2_5 = t_crit_LPSvLPS2 * math.sqrt (MSE_within_LPSvLPS2 * ((1/len(threehundred_LPS)) + (1/len(threehundred_LPS2))))

        LSD_LPSvLPS3_1 = t_crit_LPSvLPS3 * math.sqrt (MSE_within_LPSvLPS3 * ((1/len(BBS_LPS)) + (1/len(BBS_LPS3))))
        LSD_LPSvLPS3_2 = t_crit_LPSvLPS3 * math.sqrt (MSE_within_LPSvLPS3 * ((1/len(ten_LPS)) + (1/len(ten_LPS3))))
        LSD_LPSvLPS3_3 = t_crit_LPSvLPS3 * math.sqrt (MSE_within_LPSvLPS3 * ((1/len(thirty_LPS)) + (1/len(thirty_LPS3))))
        LSD_LPSvLPS3_4 = t_crit_LPSvLPS3 * math.sqrt (MSE_within_LPSvLPS3 * ((1/len(hundred_LPS)) + (1/len(hundred_LPS3))))
        LSD_LPSvLPS3_5 = t_crit_LPSvLPS3 * math.sqrt (MSE_within_LPSvLPS3 * ((1/len(threehundred_LPS)) + (1/len(threehundred_LPS3))))

        LSD_LPS2vLPS3_1 = t_crit_LPS2vLPS3 * math.sqrt (MSE_within_LPS2vLPS3 * ((1/len(BBS_LPS2)) + (1/len(BBS_LPS3))))
        LSD_LPS2vLPS3_2 = t_crit_LPS2vLPS3 * math.sqrt (MSE_within_LPS2vLPS3 * ((1/len(ten_LPS2)) + (1/len(ten_LPS3))))
        LSD_LPS2vLPS3_3 = t_crit_LPS2vLPS3 * math.sqrt (MSE_within_LPS2vLPS3 * ((1/len(thirty_LPS2)) + (1/len(thirty_LPS3))))
        LSD_LPS2vLPS3_4 = t_crit_LPS2vLPS3 * math.sqrt (MSE_within_LPS2vLPS3 * ((1/len(hundred_LPS2)) + (1/len(hundred_LPS3))))
        LSD_LPS2vLPS3_5 = t_crit_LPS2vLPS3 * math.sqrt (MSE_within_LPS2vLPS3 * ((1/len(threehundred_LPS2)) + (1/len(threehundred_LPS3))))

        LSD_CvLPS_1 = t_crit_CvLPS * math.sqrt (MSE_within_CvLPS * ((1/len(BBS)) + (1/len(BBS_LPS))))
        LSD_CvLPS_2 = t_crit_CvLPS * math.sqrt (MSE_within_CvLPS * ((1/len(ten)) + (1/len(ten_LPS))))
        LSD_CvLPS_3 = t_crit_CvLPS * math.sqrt (MSE_within_CvLPS * ((1/len(thirty)) + (1/len(thirty_LPS))))
        LSD_CvLPS_4 = t_crit_CvLPS * math.sqrt (MSE_within_CvLPS * ((1/len(hundred)) + (1/len(hundred_LPS))))
        LSD_CvLPS_5 = t_crit_CvLPS * math.sqrt (MSE_within_CvLPS * ((1/len(threehundred)) + (1/len(threehundred_LPS))))

        LSD_CvLPS2_1 = t_crit_CvLPS2 * math.sqrt (MSE_within_CvLPS2 * ((1/len(BBS)) + (1/len(BBS_LPS2))))
        LSD_CvLPS2_2 = t_crit_CvLPS2 * math.sqrt (MSE_within_CvLPS2 * ((1/len(ten)) + (1/len(ten_LPS2))))
        LSD_CvLPS2_3 = t_crit_CvLPS2 * math.sqrt (MSE_within_CvLPS2 * ((1/len(thirty)) + (1/len(thirty_LPS2))))
        LSD_CvLPS2_4 = t_crit_CvLPS2 * math.sqrt (MSE_within_CvLPS2 * ((1/len(hundred)) + (1/len(hundred_LPS2))))
        LSD_CvLPS2_5 = t_crit_CvLPS2 * math.sqrt (MSE_within_CvLPS2 * ((1/len(threehundred)) + (1/len(threehundred_LPS2))))

        LSD_CvLPS3_1 = t_crit_CvLPS3 * math.sqrt (MSE_within_CvLPS3 * ((1/len(BBS)) + (1/len(BBS_LPS3))))
        LSD_CvLPS3_2 = t_crit_CvLPS3 * math.sqrt (MSE_within_CvLPS3 * ((1/len(ten)) + (1/len(ten_LPS3))))
        LSD_CvLPS3_3 = t_crit_CvLPS3 * math.sqrt (MSE_within_CvLPS3 * ((1/len(thirty)) + (1/len(thirty_LPS3))))
        LSD_CvLPS3_4 = t_crit_CvLPS3 * math.sqrt (MSE_within_CvLPS3 * ((1/len(hundred)) + (1/len(hundred_LPS3))))
        LSD_CvLPS3_5 = t_crit_CvLPS3 * math.sqrt (MSE_within_CvLPS3 * ((1/len(threehundred)) + (1/len(threehundred_LPS3))))

        #################################################
        
        SED_control_1 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        SED_LPS2_1 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
        SED_LPS2_2 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
        SED_LPS2_3 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
        SED_LPS2_4 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

        SED_LPS3_1 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
        SED_LPS3_2 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
        SED_LPS3_3 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
        SED_LPS3_4 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

        SED_LPSvLPS2_1 = math.sqrt(MSE_within_LPSvLPS2 * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS2)) ))
        SED_LPSvLPS2_2 = math.sqrt(MSE_within_LPSvLPS2 * ( (1/len(ten_LPS)) + (1/len(ten_LPS2)) ))
        SED_LPSvLPS2_3 = math.sqrt(MSE_within_LPSvLPS2 * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS2)) ))
        SED_LPSvLPS2_4 = math.sqrt(MSE_within_LPSvLPS2 * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS2)) ))
        SED_LPSvLPS2_5 = math.sqrt(MSE_within_LPSvLPS2 * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS2)) ))

        SED_LPSvLPS3_1 = math.sqrt(MSE_within_LPSvLPS3 * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS3)) ))
        SED_LPSvLPS3_2 = math.sqrt(MSE_within_LPSvLPS3 * ( (1/len(ten_LPS)) + (1/len(ten_LPS3)) ))
        SED_LPSvLPS3_3 = math.sqrt(MSE_within_LPSvLPS3 * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS3)) ))
        SED_LPSvLPS3_4 = math.sqrt(MSE_within_LPSvLPS3 * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS3)) ))
        SED_LPSvLPS3_5 = math.sqrt(MSE_within_LPSvLPS3 * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS3)) ))

        SED_LPS2vLPS3_1 = math.sqrt(MSE_within_LPS2vLPS3 * ( (1/len(BBS_LPS2)) + (1/len(BBS_LPS3)) ))
        SED_LPS2vLPS3_2 = math.sqrt(MSE_within_LPS2vLPS3 * ( (1/len(ten_LPS2)) + (1/len(ten_LPS3)) ))
        SED_LPS2vLPS3_3 = math.sqrt(MSE_within_LPS2vLPS3 * ( (1/len(thirty_LPS2)) + (1/len(thirty_LPS3)) ))
        SED_LPS2vLPS3_4 = math.sqrt(MSE_within_LPS2vLPS3 * ( (1/len(hundred_LPS2)) + (1/len(hundred_LPS3)) ))
        SED_LPS2vLPS3_5 = math.sqrt(MSE_within_LPS2vLPS3 * ( (1/len(threehundred_LPS2)) + (1/len(threehundred_LPS3)) ))

        SED_CvLPS_1 = math.sqrt(MSE_within_CvLPS * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
        SED_CvLPS_2 = math.sqrt(MSE_within_CvLPS * ( (1/len(ten)) + (1/len(ten_LPS)) ))
        SED_CvLPS_3 = math.sqrt(MSE_within_CvLPS * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
        SED_CvLPS_4 = math.sqrt(MSE_within_CvLPS * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
        SED_CvLPS_5 = math.sqrt(MSE_within_CvLPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))

        SED_CvLPS2_1 = math.sqrt(MSE_within_CvLPS2 * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
        SED_CvLPS2_2 = math.sqrt(MSE_within_CvLPS2 * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
        SED_CvLPS2_3 = math.sqrt(MSE_within_CvLPS2 * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
        SED_CvLPS2_4 = math.sqrt(MSE_within_CvLPS2 * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
        SED_CvLPS2_5 = math.sqrt(MSE_within_CvLPS2 * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))

        SED_CvLPS3_1 = math.sqrt(MSE_within_CvLPS3 * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))
        SED_CvLPS3_2 = math.sqrt(MSE_within_CvLPS3 * ( (1/len(ten)) + (1/len(ten_LPS3)) ))
        SED_CvLPS3_3 = math.sqrt(MSE_within_CvLPS3 * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))
        SED_CvLPS3_4 = math.sqrt(MSE_within_CvLPS3 * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))
        SED_CvLPS3_5 = math.sqrt(MSE_within_CvLPS3 * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

        #################################################

        tscore_control_1 = LSD_BBS_ten / SED_control_1
        tscore_control_2 = LSD_BBS_thirty / SED_control_2
        tscore_control_3 = LSD_BBS_hundred / SED_control_3
        tscore_control_4 = LSD_BBS_threehundred / SED_control_4

        tscore_LPS_1 = LSD_BBS_ten_LPS / SED_LPS_1
        tscore_LPS_2 = LSD_BBS_thirty_LPS / SED_LPS_2
        tscore_LPS_3 = LSD_BBS_hundred_LPS / SED_LPS_3
        tscore_LPS_4 = LSD_BBS_threehundred_LPS / SED_LPS_4

        tscore_LPS2_1 = LSD_BBS_ten_LPS2 / SED_LPS2_1
        tscore_LPS2_2 = LSD_BBS_thirty_LPS2 / SED_LPS2_2
        tscore_LPS2_3 = LSD_BBS_hundred_LPS2 / SED_LPS2_3
        tscore_LPS2_4 = LSD_BBS_threehundred_LPS2 / SED_LPS2_4

        tscore_LPS3_1 = LSD_BBS_ten_LPS3 / SED_LPS3_1
        tscore_LPS3_2 = LSD_BBS_thirty_LPS3 / SED_LPS3_2
        tscore_LPS3_3 = LSD_BBS_hundred_LPS3 / SED_LPS3_3
        tscore_LPS3_4 = LSD_BBS_threehundred_LPS3 / SED_LPS3_4

        tscore_LPSvLPS2_1 = LSD_BBS_LPSvLPS2 / SED_LPSvLPS2_1
        tscore_LPSvLPS2_2 = LSD_ten_LPSvLPS2 / SED_LPSvLPS2_2
        tscore_LPSvLPS2_3 = LSD_thirty_LPSvLPS2 / SED_LPSvLPS2_3
        tscore_LPSvLPS2_4 = LSD_hundred_LPSvLPS2 / SED_LPSvLPS2_4
        tscore_LPSvLPS2_5 = LSD_threehundred_LPSvLPS2 / SED_LPSvLPS2_5

        tscore_LPSvLPS3_1 = LSD_BBS_LPSvLPS3 / SED_LPSvLPS3_1
        tscore_LPSvLPS3_2 = LSD_ten_LPSvLPS3 / SED_LPSvLPS3_2
        tscore_LPSvLPS3_3 = LSD_thirty_LPSvLPS3 / SED_LPSvLPS3_3
        tscore_LPSvLPS3_4 = LSD_hundred_LPSvLPS3 / SED_LPSvLPS3_4
        tscore_LPSvLPS3_5 = LSD_threehundred_LPSvLPS3 / SED_LPSvLPS3_5

        tscore_LPS2vLPS3_1 = LSD_BBS_LPS2vLPS3 / SED_LPS2vLPS3_1
        tscore_LPS2vLPS3_2 = LSD_ten_LPS2vLPS3 / SED_LPS2vLPS3_2
        tscore_LPS2vLPS3_3 = LSD_thirty_LPS2vLPS3 / SED_LPS2vLPS3_3
        tscore_LPS2vLPS3_4 = LSD_hundred_LPS2vLPS3 / SED_LPS2vLPS3_4
        tscore_LPS2vLPS3_5 = LSD_threehundred_LPS2vLPS3 / SED_LPS2vLPS3_5

        tscore_CvLPS_1 = LSD_BBS_CvLPS / SED_CvLPS_1
        tscore_CvLPS_2 = LSD_BBS_CvLPS / SED_CvLPS_2
        tscore_CvLPS_3 = LSD_BBS_CvLPS / SED_CvLPS_3
        tscore_CvLPS_4 = LSD_BBS_CvLPS / SED_CvLPS_4
        tscore_CvLPS_5 = LSD_BBS_CvLPS / SED_CvLPS_5

        tscore_CvLPS2_1 = LSD_BBS_CvLPS2 / SED_CvLPS2_1
        tscore_CvLPS2_2 = LSD_BBS_CvLPS2 / SED_CvLPS2_2
        tscore_CvLPS2_3 = LSD_BBS_CvLPS2 / SED_CvLPS2_3
        tscore_CvLPS2_4 = LSD_BBS_CvLPS2 / SED_CvLPS2_4
        tscore_CvLPS2_5 = LSD_BBS_CvLPS2 / SED_CvLPS2_5

        tscore_CvLPS3_1 = LSD_BBS_CvLPS3 / SED_CvLPS3_1
        tscore_CvLPS3_2 = LSD_BBS_CvLPS3 / SED_CvLPS3_2
        tscore_CvLPS3_3 = LSD_BBS_CvLPS3 / SED_CvLPS3_3
        tscore_CvLPS3_4 = LSD_BBS_CvLPS3 / SED_CvLPS3_4
        tscore_CvLPS3_5 = LSD_BBS_CvLPS3 / SED_CvLPS3_5

        #################################################

        pvalue_control_1 = 2 * scipy.stats.t.sf (tscore_control_1, df_within_control)
        pvalue_control_2 = 2 * scipy.stats.t.sf (tscore_control_2, df_within_control)
        pvalue_control_3 = 2 * scipy.stats.t.sf (tscore_control_3, df_within_control)
        pvalue_control_4 = 2 * scipy.stats.t.sf (tscore_control_4, df_within_control)

        pvalue_LPS_1 = 2 * scipy.stats.t.sf (tscore_LPS_1, df_within_LPS)
        pvalue_LPS_2 = 2 * scipy.stats.t.sf (tscore_LPS_2, df_within_LPS)
        pvalue_LPS_3 = 2 * scipy.stats.t.sf (tscore_LPS_3, df_within_LPS)
        pvalue_LPS_4 = 2 * scipy.stats.t.sf (tscore_LPS_4, df_within_LPS)

        pvalue_LPS2_1 = 2 * scipy.stats.t.sf (tscore_LPS2_1, df_within_LPS2)
        pvalue_LPS2_2 = 2 * scipy.stats.t.sf (tscore_LPS2_2, df_within_LPS2)
        pvalue_LPS2_3 = 2 * scipy.stats.t.sf (tscore_LPS2_3, df_within_LPS2)
        pvalue_LPS2_4 = 2 * scipy.stats.t.sf (tscore_LPS2_4, df_within_LPS2)

        pvalue_LPS3_1 = 2 * scipy.stats.t.sf (tscore_LPS3_1, df_within_LPS3)
        pvalue_LPS3_2 = 2 * scipy.stats.t.sf (tscore_LPS3_2, df_within_LPS3)
        pvalue_LPS3_3 = 2 * scipy.stats.t.sf (tscore_LPS3_3, df_within_LPS3)
        pvalue_LPS3_4 = 2 * scipy.stats.t.sf (tscore_LPS3_4, df_within_LPS3)

        pvalue_LPSvLPS2_1 = 2 * scipy.stats.t.sf (tscore_LPSvLPS2_1, df_within_LPSvLPS2)
        pvalue_LPSvLPS2_2 = 2 * scipy.stats.t.sf (tscore_LPSvLPS2_2, df_within_LPSvLPS2)
        pvalue_LPSvLPS2_3 = 2 * scipy.stats.t.sf (tscore_LPSvLPS2_3, df_within_LPSvLPS2)
        pvalue_LPSvLPS2_4 = 2 * scipy.stats.t.sf (tscore_LPSvLPS2_4, df_within_LPSvLPS2)
        pvalue_LPSvLPS2_5 = 2 * scipy.stats.t.sf (tscore_LPSvLPS2_5, df_within_LPSvLPS2)

        pvalue_LPSvLPS3_1 = 2 * scipy.stats.t.sf (tscore_LPSvLPS3_1, df_within_LPSvLPS3)
        pvalue_LPSvLPS3_2 = 2 * scipy.stats.t.sf (tscore_LPSvLPS3_2, df_within_LPSvLPS3)
        pvalue_LPSvLPS3_3 = 2 * scipy.stats.t.sf (tscore_LPSvLPS3_3, df_within_LPSvLPS3)
        pvalue_LPSvLPS3_4 = 2 * scipy.stats.t.sf (tscore_LPSvLPS3_4, df_within_LPSvLPS3)
        pvalue_LPSvLPS3_5 = 2 * scipy.stats.t.sf (tscore_LPSvLPS3_5, df_within_LPSvLPS3)

        pvalue_LPS2vLPS3_1 = 2 * scipy.stats.t.sf (tscore_LPS2vLPS3_1, df_within_LPS2vLPS3)
        pvalue_LPS2vLPS3_2 = 2 * scipy.stats.t.sf (tscore_LPS2vLPS3_2, df_within_LPS2vLPS3)
        pvalue_LPS2vLPS3_3 = 2 * scipy.stats.t.sf (tscore_LPS2vLPS3_3, df_within_LPS2vLPS3)
        pvalue_LPS2vLPS3_4 = 2 * scipy.stats.t.sf (tscore_LPS2vLPS3_4, df_within_LPS2vLPS3)
        pvalue_LPS2vLPS3_5 = 2 * scipy.stats.t.sf (tscore_LPS2vLPS3_5, df_within_LPS2vLPS3)

        pvalue_CvLPS_1 = 2 * scipy.stats.t.sf (tscore_CvLPS_1, df_within_CvLPS)
        pvalue_CvLPS_2 = 2 * scipy.stats.t.sf (tscore_CvLPS_2, df_within_CvLPS)
        pvalue_CvLPS_3 = 2 * scipy.stats.t.sf (tscore_CvLPS_3, df_within_CvLPS)
        pvalue_CvLPS_4 = 2 * scipy.stats.t.sf (tscore_CvLPS_4, df_within_CvLPS)
        pvalue_CvLPS_5 = 2 * scipy.stats.t.sf (tscore_CvLPS_5, df_within_CvLPS)

        pvalue_CvLPS2_1 = 2 * scipy.stats.t.sf (tscore_CvLPS2_1, df_within_CvLPS2)
        pvalue_CvLPS2_2 = 2 * scipy.stats.t.sf (tscore_CvLPS2_2, df_within_CvLPS2)
        pvalue_CvLPS2_3 = 2 * scipy.stats.t.sf (tscore_CvLPS2_3, df_within_CvLPS2)
        pvalue_CvLPS2_4 = 2 * scipy.stats.t.sf (tscore_CvLPS2_4, df_within_CvLPS2)
        pvalue_CvLPS2_5 = 2 * scipy.stats.t.sf (tscore_CvLPS2_5, df_within_CvLPS2)

        pvalue_CvLPS3_1 = 2 * scipy.stats.t.sf (tscore_CvLPS3_1, df_within_CvLPS3)
        pvalue_CvLPS3_2 = 2 * scipy.stats.t.sf (tscore_CvLPS3_2, df_within_CvLPS3)
        pvalue_CvLPS3_3 = 2 * scipy.stats.t.sf (tscore_CvLPS3_3, df_within_CvLPS3)
        pvalue_CvLPS3_4 = 2 * scipy.stats.t.sf (tscore_CvLPS3_4, df_within_CvLPS3)
        pvalue_CvLPS3_5 = 2 * scipy.stats.t.sf (tscore_CvLPS3_5, df_within_CvLPS3)

        print ("========================================================")
        print ("")
        print ("Your Fisher's LSD Results Are:")
        print ("")
        
        if LSD_BBS_ten > LSD_control_1:
            print ("BBS vs. 10:  *  | p =", pvalue_control_1)
        else:
            print ("BBS vs. 10:  ns | p =", pvalue_control_1)

        if LSD_BBS_thirty > LSD_control_2:
            print ("BBS vs. 30:  *  | p =", pvalue_control_2)
        else:
            print ("BBS vs. 30:  ns | p =", pvalue_control_2)

        if LSD_BBS_hundred > LSD_control_3:
            print ("BBS vs. 100: *  | p =", pvalue_control_3)
        else:
            print ("BBS vs. 100: ns | p =", pvalue_control_3)

        if LSD_BBS_threehundred > LSD_control_4:
            print ("BBS vs. 300: *  | p =", pvalue_control_4)
        else:
            print ("BBS vs. 300: ns | p =", pvalue_control_4)

        print ("")

        print ("////////////////////////////////////////////////////////")

        print ("")
        
        if LSD_BBS_ten_LPS > LSD_LPS_1:
            print ("1st Dose BBS LPS vs. 10 LPS:  *  | p =", pvalue_LPS_1)
        else:
            print ("1st Dose BBS LPS vs. 10 LPS:  ns | p =", pvalue_LPS_1)

        if LSD_BBS_thirty_LPS > LSD_LPS_2:
            print ("1st Dose BBS LPS vs. 30 LPS:  *  | p =", pvalue_LPS_2)
        else:
            print ("1st Dose BBS LPS vs. 30 LPS:  ns | p =", pvalue_LPS_2)

        if LSD_BBS_hundred_LPS > LSD_LPS_3:
            print ("1st Dose BBS LPS vs. 100 LPS: *  | p =", pvalue_LPS_3)
        else:
            print ("1st Dose BBS LPS vs. 100 LPS: ns | p =", pvalue_LPS_3)

        if LSD_BBS_threehundred_LPS > LSD_LPS_4:
            print ("1st Dose BBS LPS vs. 300 LPS: *  | p =", pvalue_LPS_4)
        else:
            print ("1st Dose BBS LPS vs. 300 LPS: ns | p =", pvalue_LPS_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        
        if LSD_BBS_ten_LPS2 > LSD_LPS2_1:
            print ("2nd Dose BBS LPS vs. 10 LPS:  *  | p =", pvalue_LPS2_1)
        else:
            print ("2nd Dose BBS LPS vs. 10 LPS:  ns | p =", pvalue_LPS2_1)

        if LSD_BBS_thirty_LPS2 > LSD_LPS2_2:
            print ("2nd Dose BBS LPS vs. 30 LPS:  *  | p =", pvalue_LPS2_2)
        else:
            print ("2nd Dose BBS LPS vs. 30 LPS:  ns | p =", pvalue_LPS2_2)

        if LSD_BBS_hundred_LPS2 > LSD_LPS2_3:
            print ("2nd Dose BBS LPS vs. 100 LPS: *  | p =", pvalue_LPS2_3)
        else:
            print ("2nd Dose BBS LPS vs. 100 LPS: ns | p =", pvalue_LPS2_3)

        if LSD_BBS_threehundred_LPS2 > LSD_LPS2_4:
            print ("2nd Dose BBS LPS vs. 300 LPS: *  | p =", pvalue_LPS2_4)
        else:
            print ("2nd Dose BBS LPS vs. 300 LPS: ns | p =", pvalue_LPS2_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        
        if LSD_BBS_ten_LPS3 > LSD_LPS3_1:
            print ("3rd Dose BBS LPS vs. 10 LPS:  *  | p =", pvalue_LPS3_1)
        else:
            print ("3rd Dose BBS LPS vs. 10 LPS:  ns | p =", pvalue_LPS3_1)

        if LSD_BBS_thirty_LPS3 > LSD_LPS3_2:
            print ("3rd Dose BBS LPS vs. 30 LPS:  *  | p =", pvalue_LPS3_2)
        else:
            print ("3rd Dose BBS LPS vs. 30 LPS:  ns | p =", pvalue_LPS2_2)

        if LSD_BBS_hundred_LPS3 > LSD_LPS3_3:
            print ("3rd Dose BBS LPS vs. 100 LPS: *  | p =", pvalue_LPS3_3)
        else:
            print ("3rd Dose BBS LPS vs. 100 LPS: ns | p =", pvalue_LPS2_3)

        if LSD_BBS_threehundred_LPS3 > LSD_LPS3_4:
            print ("3rd Dose BBS LPS vs. 300 LPS: *  | p =", pvalue_LPS3_4)
        else:
            print ("3rd Dose BBS LPS vs. 300 LPS: ns | p =", pvalue_LPS3_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if LSD_BBS_CvLPS > LSD_CvLPS_1:
            print ("BBS vs. 1st Dose BBS LPS:  * | p =", pvalue_CvLPS_1)
        else:
            print ("BBS vs. 1st Dose BBS LPS: ns | p =", pvalue_CvLPS_1)

        if LSD_BBS_CvLPS > LSD_CvLPS_2:
            print ("10 vs. 1st Dose 10 LPS:    * | p =", pvalue_CvLPS_2)
        else:
            print ("10 vs. 1st Dose 10 LPS:   ns | p =", pvalue_CvLPS_2)

        if LSD_BBS_CvLPS > LSD_CvLPS_3:
            print ("30 vs. 1st Dose 30 LPS:    * | p =", pvalue_CvLPS_3)
        else:
            print ("30 vs. 1st Dose 30 LPS:   ns | p =", pvalue_CvLPS_3)

        if LSD_BBS_CvLPS > LSD_CvLPS_4:
            print ("100 vs. 1st Dose 100 LPS:  * | p =", pvalue_CvLPS_4)
        else:
            print ("100 vs. 1st Dose 100 LPS: ns | p =", pvalue_CvLPS_4)

        if LSD_BBS_CvLPS > LSD_CvLPS_5:
            print ("300 vs. 1st Dose 300 LPS:  * | p =", pvalue_CvLPS_5)
        else:
            print ("300 vs. 1st Dose 300 LPS: ns | p =", pvalue_CvLPS_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if LSD_BBS_CvLPS2 > LSD_CvLPS2_1:
            print ("BBS vs. 2nd Dose BBS LPS:  * | p =", pvalue_CvLPS2_1)
        else:
            print ("BBS vs. 2nd Dose BBS LPS: ns | p =", pvalue_CvLPS2_1)

        if LSD_BBS_CvLPS2 > LSD_CvLPS2_2:
            print ("10 vs. 2nd Dose 10 LPS:    * | p =", pvalue_CvLPS2_2)
        else:
            print ("10 vs. 2nd Dose 10 LPS:   ns | p =", pvalue_CvLPS2_2)

        if LSD_BBS_CvLPS2 > LSD_CvLPS2_3:
            print ("30 vs. 2nd Dose 30 LPS:    * | p =", pvalue_CvLPS2_3)
        else:
            print ("30 vs. 2nd Dose 30 LPS:   ns | p =", pvalue_CvLPS2_3)

        if LSD_BBS_CvLPS2 > LSD_CvLPS2_4:
            print ("100 vs. 2nd Dose 100 LPS:  * | p =", pvalue_CvLPS2_4)
        else:
            print ("100 vs. 2nd Dose 100 LPS: ns | p =", pvalue_CvLPS2_4)

        if LSD_BBS_CvLPS2 > LSD_CvLPS2_5:
            print ("300 vs. 2nd Dose 300 LPS:  * | p =", pvalue_CvLPS2_5)
        else:
            print ("300 vs. 2nd Dose 300 LPS: ns | p =", pvalue_CvLPS2_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if LSD_BBS_CvLPS3 > LSD_CvLPS3_1:
            print ("BBS vs. 3rd Dose BBS LPS:  * | p =", pvalue_CvLPS3_1)
        else:
            print ("BBS vs. 3rd Dose BBS LPS: ns | p =", pvalue_CvLPS3_1)

        if LSD_BBS_CvLPS3 > LSD_CvLPS3_2:
            print ("10 vs. 3rd Dose 10 LPS:    * | p =", pvalue_CvLPS3_2)
        else:
            print ("10 vs. 3rd Dose 10 LPS:   ns | p =", pvalue_CvLPS3_2)

        if LSD_BBS_CvLPS3 > LSD_CvLPS3_3:
            print ("30 vs. 3rd Dose 30 LPS:    * | p =", pvalue_CvLPS3_3)
        else:
            print ("30 vs. 3rd Dose 30 LPS:   ns | p =", pvalue_CvLPS3_3)

        if LSD_BBS_CvLPS3 > LSD_CvLPS3_4:
            print ("100 vs. 3rd Dose 100 LPS:  * | p =", pvalue_CvLPS3_4)
        else:
            print ("100 vs. 3rd Dose 100 LPS: ns | p =", pvalue_CvLPS3_4)

        if LSD_BBS_CvLPS3 > LSD_CvLPS3_5:
            print ("300 vs. 3rd Dose 300 LPS:  * | p =", pvalue_CvLPS3_5)
        else:
            print ("300 vs. 3rd Dose 300 LPS: ns | p =", pvalue_CvLPS3_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")
        
        if LSD_BBS_LPSvLPS2 > LSD_LPSvLPS2_1:
            print ("1st Dose BBS LPS vs. 2nd Dose BBS LPS:  *  | p =", pvalue_LPSvLPS2_1)
        else:
            print ("1st Dose BBS LPS vs. 2nd Dose BBS LPS:  ns | p =", pvalue_LPSvLPS2_1)
            
        if LSD_BBS_LPSvLPS2 > LSD_LPSvLPS2_2:
            print ("1st Dose 10 LPS vs. 2nd Dose 10 LPS:    *  | p =", pvalue_LPSvLPS2_2)
        else:
            print ("1st Dose 10 LPS vs. 2nd Dose 10 LPS:    ns | p =", pvalue_LPSvLPS2_2)

        if LSD_BBS_LPSvLPS2 > LSD_LPSvLPS2_3:
            print ("1st Dose 30 LPS vs. 2nd Dose 30 LPS:    *  | p =", pvalue_LPSvLPS2_3)
        else:
            print ("1st Dose 30 LPS vs. 2nd Dose 30 LPS:    ns | p =", pvalue_LPSvLPS2_3)

        if LSD_BBS_LPSvLPS2 > LSD_LPSvLPS2_4:
            print ("1st Dose 100 LPS vs. 2nd Dose 100 LPS:  *  | p =", pvalue_LPSvLPS2_4)
        else:
            print ("1st Dose 100 LPS vs. 2nd Dose 100 LPS:  ns | p =", pvalue_LPSvLPS2_4)

        if LSD_BBS_LPSvLPS2 > LSD_LPSvLPS2_5:
            print ("1st Dose 300 LPS vs. 2nd Dose 300 LPS:  *  | p =", pvalue_LPSvLPS2_5)
        else:
            print ("1st Dose 300 LPS vs. 2nd Dose 300 LPS:  ns | p =", pvalue_LPSvLPS2_5)    

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if LSD_BBS_LPSvLPS3 > LSD_LPSvLPS3_1:
            print ("1st Dose BBS LPS vs. 3rd Dose BBS LPS:  *  | p =", pvalue_LPSvLPS3_1)
        else:
            print ("1st Dose BBS LPS vs. 3rd Dose BBS LPS:  ns | p =", pvalue_LPSvLPS3_1)
            
        if LSD_BBS_LPSvLPS3 > LSD_LPSvLPS3_2:
            print ("1st Dose 10 LPS vs. 3rd Dose 10 LPS:    *  | p =", pvalue_LPSvLPS3_2)
        else:
            print ("1st Dose 10 LPS vs. 3rd Dose 10 LPS:    ns | p =", pvalue_LPSvLPS3_2)

        if LSD_BBS_LPSvLPS3 > LSD_LPSvLPS3_3:
            print ("1st Dose 30 LPS vs. 3rd Dose 30 LPS:    *  | p =", pvalue_LPSvLPS3_3)
        else:
            print ("1st Dose 30 LPS vs. 3rd Dose 30 LPS:    ns | p =", pvalue_LPSvLPS3_3)

        if LSD_BBS_LPSvLPS3 > LSD_LPSvLPS3_4:
            print ("1st Dose 100 LPS vs. 3rd Dose 100 LPS:  *  | p =", pvalue_LPSvLPS3_4)
        else:
            print ("1st Dose 100 LPS vs. 3rd Dose 100 LPS:  ns | p =", pvalue_LPSvLPS3_4)

        if LSD_BBS_LPSvLPS3 > LSD_LPSvLPS3_5:
            print ("1st Dose 300 LPS vs. 3rd Dose 300 LPS:  *  | p =", pvalue_LPSvLPS3_5)
        else:
            print ("1st Dose 300 LPS vs. 3rd Dose 300 LPS:  ns | p =", pvalue_LPSvLPS3_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if LSD_BBS_LPS2vLPS3 > LSD_LPS2vLPS3_1:
            print ("2nd Dose BBS LPS vs. 3rd Dose BBS LPS:  *  | p =", pvalue_LPS2vLPS3_1)
        else:
            print ("2nd Dose BBS LPS vs. 3rd Dose BBS LPS:  ns | p =", pvalue_LPS2vLPS3_1)
            
        if LSD_BBS_LPS2vLPS3 > LSD_LPS2vLPS3_2:
            print ("2nd Dose 10 LPS vs. 3rd Dose 10 LPS:    *  | p =", pvalue_LPS2vLPS3_2)
        else:
            print ("2nd Dose 10 LPS vs. 3rd Dose 10 LPS:    ns | p =", pvalue_LPS2vLPS3_2)

        if LSD_BBS_LPS2vLPS3 > LSD_LPS2vLPS3_3:
            print ("2nd Dose 30 LPS vs. 3rd Dose 30 LPS:    *  | p =", pvalue_LPS2vLPS3_3)
        else:
            print ("2nd Dose 30 LPS vs. 3rd Dose 30 LPS:    ns | p =", pvalue_LPS2vLPS3_3)

        if LSD_BBS_LPS2vLPS3 > LSD_LPS2vLPS3_4:
            print ("2nd Dose 100 LPS vs. 3rd Dose 100 LPS:  *  | p =", pvalue_LPS2vLPS3_4)
        else:
            print ("2nd Dose 100 LPS vs. 3rd Dose 100 LPS:  ns | p =", pvalue_LPS2vLPS3_4)

        if LSD_BBS_LPS2vLPS3 > LSD_LPS2vLPS3_5:
            print ("2nd Dose 300 LPS vs. 3rd Dose 300 LPS:  *  | p =", pvalue_LPS2vLPS3_5)
        else:
            print ("2nd Dose 300 LPS vs. 3rd Dose 300 LPS:  ns | p =", pvalue_LPS2vLPS3_5)

        print ("")
        print ("========================================================")
        print ("")
        
        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("You must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")

    elif answer == "2":

        #Tukey's HSD: if the HSD is higher than Q critical value, il est significant

        from statsmodels.stats.libqsturng import psturng
        
        #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/Q critical values.csv", header = None)
        crit_table = [37.082,	10.881,	7.502,	6.287,	5.673,	5.305,	5.06,	4.886,	4.755,	4.654,	4.574,	4.508,	4.453,	4.407,	4.367,	4.333,	4.303,	4.276,	4.253,	4.232,	4.213,	4.196,	4.18,	4.166,	4.153,	4.141,	4.13,	4.12,	4.111,	4.102,	4.094,	4.086,	4.079,	4.072,	4.066,	4.06,	4.054,	4.049,	4.044,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887]
        crit_table2 = [49.071,13.988,9.462,7.826,6.995,6.493,6.158,5.918,5.738,5.598,5.486,5.395,5.318,5.253,5.198,5.15,5.108,5.071,5.037,5.008,4.981,4.957,4.935,4.915,4.897,4.88,4.864,4.85,4.837,4.824,4.812,4.802,4.791,4.782,4.773,4.764,4.756,4.749,4.741,4.735,4.735,4.735,4.735,4.735,4.735,4.735,4.735,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517]

        #################################################

        HSD_BBS_ten = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_thirty = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_hundred = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_threehundred = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred)))) / (math.sqrt (MSE_within_control / float(N)))

        HSD_BBS_ten_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_thirty_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_hundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_threehundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))

        HSD_BBS_ten_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (ten_LPS2) / len (ten_LPS2)))) / (math.sqrt (MSE_within_LPS2 / float(N)))
        HSD_BBS_thirty_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (thirty_LPS2) / len (thirty_LPS2)))) / (math.sqrt (MSE_within_LPS2 / float(N)))
        HSD_BBS_hundred_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (hundred_LPS2) / len (hundred_LPS2)))) / (math.sqrt (MSE_within_LPS2 / float(N)))
        HSD_BBS_threehundred_LPS2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (threehundred_LPS2) / len (threehundred_LPS2)))) / (math.sqrt (MSE_within_LPS2 / float(N)))

        HSD_BBS_ten_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (ten_LPS3) / len (ten_LPS3)))) / (math.sqrt (MSE_within_LPS3 / float(N)))
        HSD_BBS_thirty_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (thirty_LPS3) / len (thirty_LPS3)))) / (math.sqrt (MSE_within_LPS3 / float(N)))
        HSD_BBS_hundred_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (hundred_LPS3) / len (hundred_LPS3)))) / (math.sqrt (MSE_within_LPS3 / float(N)))
        HSD_BBS_threehundred_LPS3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (threehundred_LPS3) / len (threehundred_LPS3)))) / (math.sqrt (MSE_within_LPS3 / float(N)))

        HSD_LPSvLPS2_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS2) / len (BBS_LPS2)))) / (math.sqrt (MSE_within_LPSvLPS2 / float(N)))
        HSD_LPSvLPS2_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS2) / len (ten_LPS2)))) / (math.sqrt (MSE_within_LPSvLPS2 / float(N)))
        HSD_LPSvLPS2_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS2) / len (thirty_LPS2)))) / (math.sqrt (MSE_within_LPSvLPS2 / float(N)))
        HSD_LPSvLPS2_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS2) / len (hundred_LPS2)))) / (math.sqrt (MSE_within_LPSvLPS2 / float(N)))
        HSD_LPSvLPS2_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS2) / len (threehundred_LPS2)))) / (math.sqrt (MSE_within_LPSvLPS2 / float(N)))

        HSD_LPSvLPS3_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS3) / len (BBS_LPS3)))) / (math.sqrt (MSE_within_LPSvLPS3 / float(N)))
        HSD_LPSvLPS3_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS3) / len (ten_LPS3)))) / (math.sqrt (MSE_within_LPSvLPS3 / float(N)))
        HSD_LPSvLPS3_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS3) / len (thirty_LPS3)))) / (math.sqrt (MSE_within_LPSvLPS3 / float(N)))
        HSD_LPSvLPS3_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS3) / len (hundred_LPS3)))) / (math.sqrt (MSE_within_LPSvLPS3 / float(N)))
        HSD_LPSvLPS3_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS3) / len (threehundred_LPS3)))) / (math.sqrt (MSE_within_LPSvLPS3 / float(N)))

        HSD_LPS2vLPS3_1 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (BBS_LPS3) / len (BBS_LPS3)))) / (math.sqrt (MSE_within_LPS2vLPS3 / float(N)))
        HSD_LPS2vLPS3_2 = abs (((sum (ten_LPS2) / len (ten_LPS2)) - (sum (ten_LPS3) / len (ten_LPS3)))) / (math.sqrt (MSE_within_LPS2vLPS3 / float(N)))
        HSD_LPS2vLPS3_3 = abs (((sum (thirty_LPS2) / len (thirty_LPS2)) - (sum (thirty_LPS3) / len (thirty_LPS3)))) / (math.sqrt (MSE_within_LPS2vLPS3 / float(N)))
        HSD_LPS2vLPS3_4 = abs (((sum (hundred_LPS2) / len (hundred_LPS2)) - (sum (hundred_LPS3) / len (hundred_LPS3)))) / (math.sqrt (MSE_within_LPS2vLPS3 / float(N)))
        HSD_LPS2vLPS3_5 = abs (((sum (threehundred_LPS2) / len (threehundred_LPS2)) - (sum (threehundred_LPS3) / len (threehundred_LPS3)))) / (math.sqrt (MSE_within_LPS2vLPS3 / float(N)))

        HSD_CvLPS_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS) / len (BBS_LPS)))) / (math.sqrt (MSE_within_CvLPS / float(N)))
        HSD_CvLPS_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS) / len (ten_LPS)))) / (math.sqrt (MSE_within_CvLPS / float(N)))
        HSD_CvLPS_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS) / len (thirty_LPS)))) / (math.sqrt (MSE_within_CvLPS / float(N)))
        HSD_CvLPS_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS) / len (hundred_LPS)))) / (math.sqrt (MSE_within_CvLPS / float(N)))
        HSD_CvLPS_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS) / len (threehundred_LPS)))) / (math.sqrt (MSE_within_CvLPS / float(N)))

        HSD_CvLPS2_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS2) / len (BBS_LPS2)))) / (math.sqrt (MSE_within_CvLPS2 / float(N)))
        HSD_CvLPS2_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS2) / len (ten_LPS2)))) / (math.sqrt (MSE_within_CvLPS2 / float(N)))
        HSD_CvLPS2_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS2) / len (thirty_LPS2)))) / (math.sqrt (MSE_within_CvLPS2 / float(N)))
        HSD_CvLPS2_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS2) / len (hundred_LPS2)))) / (math.sqrt (MSE_within_CvLPS2 / float(N)))
        HSD_CvLPS2_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS2) / len (threehundred_LPS2)))) / (math.sqrt (MSE_within_CvLPS2 / float(N)))

        HSD_CvLPS3_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS3) / len (BBS_LPS3)))) / (math.sqrt (MSE_within_CvLPS3 / float(N)))
        HSD_CvLPS3_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS3) / len (ten_LPS3)))) / (math.sqrt (MSE_within_CvLPS3 / float(N)))
        HSD_CvLPS3_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS3) / len (thirty_LPS3)))) / (math.sqrt (MSE_within_CvLPS3 / float(N)))
        HSD_CvLPS3_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS3) / len (hundred_LPS3)))) / (math.sqrt (MSE_within_CvLPS3 / float(N)))
        HSD_CvLPS3_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS3) / len (threehundred_LPS3)))) / (math.sqrt (MSE_within_CvLPS3 / float(N)))
        
        #################################################

        #HSD_control = crit_table.iloc[int (float(df_within_control)-1), 0]
        #HSD_LPS = crit_table.iloc[int (float(df_within_LPS)-1), 0]

        HSD_control = crit_table[int (float(df_within_control)-1)]
        HSD_LPS = crit_table[int (float(df_within_LPS)-1)]
        HSD_LPS2 = crit_table[int (float(df_within_LPS2)-1)]
        HSD_LPS3 = crit_table[int (float(df_within_LPS3)-1)]
        HSD_LPSvLPS2 = crit_table2[int (float(df_within_LPSvLPS2)-1)]
        HSD_LPSvLPS3 = crit_table2[int (float(df_within_LPSvLPS3)-1)]
        HSD_LPS2vLPS3 = crit_table2[int (float(df_within_LPS2vLPS3)-1)]
        HSD_CvLPS = crit_table2[int (float(df_within_CvLPS)-1)]
        HSD_CvLPS2 = crit_table2[int (float(df_within_CvLPS2)-1)]
        HSD_CvLPS3 = crit_table2[int (float(df_within_CvLPS3)-1)]

        #################################################

        MD_control_1 = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        MD_control_2 = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        MD_control_3 = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        MD_control_4 = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))

        MD_LPS_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_LPS_2 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_LPS_3 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_LPS_4 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        MD_LPS2_1 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_LPS2_2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_LPS2_3 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_LPS2_4 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_LPS3_1 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPS3_2 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPS3_3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPS3_4 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_LPSvLPS2_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        MD_LPSvLPS2_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_LPSvLPS2_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_LPSvLPS2_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_LPSvLPS2_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_LPSvLPS3_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_LPSvLPS3_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPSvLPS3_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPSvLPS3_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPSvLPS3_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_LPS2vLPS3_1 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_LPS2vLPS3_2 = abs (((sum (ten_LPS2) / len (ten_LPS2)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPS2vLPS3_3 = abs (((sum (thirty_LPS2) / len (thirty_LPS2)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPS2vLPS3_4 = abs (((sum (hundred_LPS2) / len (hundred_LPS2)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPS2vLPS3_5 = abs (((sum (threehundred_LPS2) / len (threehundred_LPS2)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_CvLPS_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS) / len (BBS_LPS))))
        MD_CvLPS_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_CvLPS_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_CvLPS_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_CvLPS_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        MD_CvLPS2_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        MD_CvLPS2_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_CvLPS2_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_CvLPS2_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_CvLPS2_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_CvLPS3_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_CvLPS3_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_CvLPS3_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_CvLPS3_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_CvLPS3_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        #################################################

        SED_control_1 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        SED_LPS2_1 = math.sqrt (MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
        SED_LPS2_2 = math.sqrt (MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
        SED_LPS2_3 = math.sqrt (MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
        SED_LPS2_4 = math.sqrt (MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

        SED_LPS3_1 = math.sqrt (MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
        SED_LPS3_2 = math.sqrt (MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
        SED_LPS3_3 = math.sqrt (MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
        SED_LPS3_4 = math.sqrt (MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

        SED_LPSvLPS2_1 = math.sqrt (MSE_within_LPSvLPS2 * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS2)) ))
        SED_LPSvLPS2_2 = math.sqrt (MSE_within_LPSvLPS2 * ( (1/len(ten_LPS)) + (1/len(ten_LPS2)) ))
        SED_LPSvLPS2_3 = math.sqrt (MSE_within_LPSvLPS2 * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS2)) ))
        SED_LPSvLPS2_4 = math.sqrt (MSE_within_LPSvLPS2 * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS2)) ))
        SED_LPSvLPS2_5 = math.sqrt (MSE_within_LPSvLPS2 * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS2)) ))

        SED_LPSvLPS3_1 = math.sqrt (MSE_within_LPSvLPS3 * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS3)) ))
        SED_LPSvLPS3_2 = math.sqrt (MSE_within_LPSvLPS3 * ( (1/len(ten_LPS)) + (1/len(ten_LPS3)) ))
        SED_LPSvLPS3_3 = math.sqrt (MSE_within_LPSvLPS3 * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS3)) ))
        SED_LPSvLPS3_4 = math.sqrt (MSE_within_LPSvLPS3 * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS3)) ))
        SED_LPSvLPS3_5 = math.sqrt (MSE_within_LPSvLPS3 * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS3)) ))
        
        SED_LPS2vLPS3_1 = math.sqrt (MSE_within_LPS2vLPS3 * ( (1/len(BBS_LPS2)) + (1/len(BBS_LPS3)) ))
        SED_LPS2vLPS3_2 = math.sqrt (MSE_within_LPS2vLPS3 * ( (1/len(ten_LPS2)) + (1/len(ten_LPS3)) ))
        SED_LPS2vLPS3_3 = math.sqrt (MSE_within_LPS2vLPS3 * ( (1/len(thirty_LPS2)) + (1/len(thirty_LPS3)) ))
        SED_LPS2vLPS3_4 = math.sqrt (MSE_within_LPS2vLPS3 * ( (1/len(hundred_LPS2)) + (1/len(hundred_LPS3)) ))
        SED_LPS2vLPS3_5 = math.sqrt (MSE_within_LPS2vLPS3 * ( (1/len(threehundred_LPS2)) + (1/len(threehundred_LPS3)) ))

        SED_CvLPS_1 = math.sqrt (MSE_within_CvLPS * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
        SED_CvLPS_2 = math.sqrt (MSE_within_CvLPS * ( (1/len(ten)) + (1/len(ten_LPS)) ))
        SED_CvLPS_3 = math.sqrt (MSE_within_CvLPS * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
        SED_CvLPS_4 = math.sqrt (MSE_within_CvLPS * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
        SED_CvLPS_5 = math.sqrt (MSE_within_CvLPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))

        SED_CvLPS2_1 = math.sqrt (MSE_within_CvLPS2 * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
        SED_CvLPS2_2 = math.sqrt (MSE_within_CvLPS2 * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
        SED_CvLPS2_3 = math.sqrt (MSE_within_CvLPS2 * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
        SED_CvLPS2_4 = math.sqrt (MSE_within_CvLPS2 * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
        SED_CvLPS2_5 = math.sqrt (MSE_within_CvLPS2 * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))

        SED_CvLPS3_1 = math.sqrt (MSE_within_CvLPS3 * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))
        SED_CvLPS3_2 = math.sqrt (MSE_within_CvLPS3 * ( (1/len(ten)) + (1/len(ten_LPS3)) ))
        SED_CvLPS3_3 = math.sqrt (MSE_within_CvLPS3 * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))
        SED_CvLPS3_4 = math.sqrt (MSE_within_CvLPS3 * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))
        SED_CvLPS3_5 = math.sqrt (MSE_within_CvLPS3 * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

        #################################################

        q_control_1 = MD_control_1 / (SED_control_1 * 1/math.sqrt(2))
        q_control_2 = MD_control_2 / (SED_control_2 * 1/math.sqrt(2))
        q_control_3 = MD_control_3 / (SED_control_3 * 1/math.sqrt(2))
        q_control_4 = MD_control_4 / (SED_control_4 * 1/math.sqrt(2))

        q_LPS_1 = MD_LPS_1 / (SED_LPS_1 * 1/math.sqrt(2))
        q_LPS_2 = MD_LPS_2 / (SED_LPS_2 * 1/math.sqrt(2))
        q_LPS_3 = MD_LPS_3 / (SED_LPS_3 * 1/math.sqrt(2))
        q_LPS_4 = MD_LPS_4 / (SED_LPS_4 * 1/math.sqrt(2))

        q_LPS2_1 = MD_LPS2_1 / (SED_LPS2_1 * 1/math.sqrt(2))
        q_LPS2_2 = MD_LPS2_2 / (SED_LPS2_2 * 1/math.sqrt(2))
        q_LPS2_3 = MD_LPS2_3 / (SED_LPS2_3 * 1/math.sqrt(2))
        q_LPS2_4 = MD_LPS2_4 / (SED_LPS2_4 * 1/math.sqrt(2))

        q_LPS3_1 = MD_LPS3_1 / (SED_LPS3_1 * 1/math.sqrt(2))
        q_LPS3_2 = MD_LPS3_2 / (SED_LPS3_2 * 1/math.sqrt(2))
        q_LPS3_3 = MD_LPS3_3 / (SED_LPS3_3 * 1/math.sqrt(2))
        q_LPS3_4 = MD_LPS3_4 / (SED_LPS3_4 * 1/math.sqrt(2))

        q_LPSvLPS2_1 = MD_LPSvLPS2_1 / (SED_LPSvLPS2_1 * 1/math.sqrt(2))
        q_LPSvLPS2_2 = MD_LPSvLPS2_2 / (SED_LPSvLPS2_2 * 1/math.sqrt(2))
        q_LPSvLPS2_3 = MD_LPSvLPS2_3 / (SED_LPSvLPS2_3 * 1/math.sqrt(2))
        q_LPSvLPS2_4 = MD_LPSvLPS2_4 / (SED_LPSvLPS2_4 * 1/math.sqrt(2))
        q_LPSvLPS2_5 = MD_LPSvLPS2_5 / (SED_LPSvLPS2_5 * 1/math.sqrt(2))

        q_LPSvLPS3_1 = MD_LPSvLPS3_1 / (SED_LPSvLPS3_1 * 1/math.sqrt(2))
        q_LPSvLPS3_2 = MD_LPSvLPS3_2 / (SED_LPSvLPS3_2 * 1/math.sqrt(2))
        q_LPSvLPS3_3 = MD_LPSvLPS3_3 / (SED_LPSvLPS3_3 * 1/math.sqrt(2))
        q_LPSvLPS3_4 = MD_LPSvLPS3_4 / (SED_LPSvLPS3_4 * 1/math.sqrt(2))
        q_LPSvLPS3_5 = MD_LPSvLPS3_5 / (SED_LPSvLPS3_5 * 1/math.sqrt(2))

        q_LPS2vLPS3_1 = MD_LPS2vLPS3_1 / (SED_LPS2vLPS3_1 * 1/math.sqrt(2))
        q_LPS2vLPS3_2 = MD_LPS2vLPS3_2 / (SED_LPS2vLPS3_2 * 1/math.sqrt(2))
        q_LPS2vLPS3_3 = MD_LPS2vLPS3_3 / (SED_LPS2vLPS3_3 * 1/math.sqrt(2))
        q_LPS2vLPS3_4 = MD_LPS2vLPS3_4 / (SED_LPS2vLPS3_4 * 1/math.sqrt(2))
        q_LPS2vLPS3_5 = MD_LPS2vLPS3_5 / (SED_LPS2vLPS3_5 * 1/math.sqrt(2))

        q_CvLPS_1 = MD_CvLPS_1 / (SED_CvLPS_1 * 1/math.sqrt(2))
        q_CvLPS_2 = MD_CvLPS_2 / (SED_CvLPS_2 * 1/math.sqrt(2))
        q_CvLPS_3 = MD_CvLPS_3 / (SED_CvLPS_3 * 1/math.sqrt(2))
        q_CvLPS_4 = MD_CvLPS_4 / (SED_CvLPS_4 * 1/math.sqrt(2))
        q_CvLPS_5 = MD_CvLPS_5 / (SED_CvLPS_5 * 1/math.sqrt(2))

        q_CvLPS2_1 = MD_CvLPS2_1 / (SED_CvLPS2_1 * 1/math.sqrt(2))
        q_CvLPS2_2 = MD_CvLPS2_2 / (SED_CvLPS2_2 * 1/math.sqrt(2))
        q_CvLPS2_3 = MD_CvLPS2_3 / (SED_CvLPS2_3 * 1/math.sqrt(2))
        q_CvLPS2_4 = MD_CvLPS2_4 / (SED_CvLPS2_4 * 1/math.sqrt(2))
        q_CvLPS2_5 = MD_CvLPS2_5 / (SED_CvLPS2_5 * 1/math.sqrt(2))

        q_CvLPS3_1 = MD_CvLPS3_1 / (SED_CvLPS3_1 * 1/math.sqrt(2))
        q_CvLPS3_2 = MD_CvLPS3_2 / (SED_CvLPS3_2 * 1/math.sqrt(2))
        q_CvLPS3_3 = MD_CvLPS3_3 / (SED_CvLPS3_3 * 1/math.sqrt(2))
        q_CvLPS3_4 = MD_CvLPS3_4 / (SED_CvLPS3_4 * 1/math.sqrt(2))
        q_CvLPS3_5 = MD_CvLPS3_5 / (SED_CvLPS3_5 * 1/math.sqrt(2))
        
        #################################################

        p_control_1 = psturng (q_control_1, 5, df_within_control)
        p_control_2 = psturng (q_control_2, 5, df_within_control)
        p_control_3 = psturng (q_control_3, 5, df_within_control)
        p_control_4 = psturng (q_control_4, 5, df_within_control)
        
        p_LPS_1 = psturng (q_LPS_1, 5, df_within_LPS)
        p_LPS_2 = psturng (q_LPS_2, 5, df_within_LPS)
        p_LPS_3 = psturng (q_LPS_3, 5, df_within_LPS)
        p_LPS_4 = psturng (q_LPS_4, 5, df_within_LPS)

        p_LPS2_1 = psturng (q_LPS2_1, 5, df_within_LPS2)
        p_LPS2_2 = psturng (q_LPS2_2, 5, df_within_LPS2)
        p_LPS2_3 = psturng (q_LPS2_3, 5, df_within_LPS2)
        p_LPS2_4 = psturng (q_LPS2_4, 5, df_within_LPS2)

        p_LPS3_1 = psturng (q_LPS3_1, 5, df_within_LPS3)
        p_LPS3_2 = psturng (q_LPS3_2, 5, df_within_LPS3)
        p_LPS3_3 = psturng (q_LPS3_3, 5, df_within_LPS3)
        p_LPS3_4 = psturng (q_LPS3_4, 5, df_within_LPS3)

        p_LPSvLPS2_1 = psturng (q_LPSvLPS2_1, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_2 = psturng (q_LPSvLPS2_2, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_3 = psturng (q_LPSvLPS2_3, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_4 = psturng (q_LPSvLPS2_4, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_5 = psturng (q_LPSvLPS2_5, 10, df_within_LPSvLPS2)

        p_LPSvLPS3_1 = psturng (q_LPSvLPS3_1, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_2 = psturng (q_LPSvLPS3_2, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_3 = psturng (q_LPSvLPS3_3, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_4 = psturng (q_LPSvLPS3_4, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_5 = psturng (q_LPSvLPS3_5, 10, df_within_LPSvLPS3)

        p_LPS2vLPS3_1 = psturng (q_LPS2vLPS3_1, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_2 = psturng (q_LPS2vLPS3_2, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_3 = psturng (q_LPS2vLPS3_3, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_4 = psturng (q_LPS2vLPS3_4, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_5 = psturng (q_LPS2vLPS3_5, 10, df_within_LPS2vLPS3)

        p_CvLPS_1 = psturng (q_CvLPS_1, 10, df_within_CvLPS)
        p_CvLPS_2 = psturng (q_CvLPS_2, 10, df_within_CvLPS)
        p_CvLPS_3 = psturng (q_CvLPS_3, 10, df_within_CvLPS)
        p_CvLPS_4 = psturng (q_CvLPS_4, 10, df_within_CvLPS)
        p_CvLPS_5 = psturng (q_CvLPS_5, 10, df_within_CvLPS)

        p_CvLPS2_1 = psturng (q_CvLPS2_1, 10, df_within_CvLPS2)
        p_CvLPS2_2 = psturng (q_CvLPS2_2, 10, df_within_CvLPS2)
        p_CvLPS2_3 = psturng (q_CvLPS2_3, 10, df_within_CvLPS2)
        p_CvLPS2_4 = psturng (q_CvLPS2_4, 10, df_within_CvLPS2)
        p_CvLPS2_5 = psturng (q_CvLPS2_5, 10, df_within_CvLPS2)

        p_CvLPS3_1 = psturng (q_CvLPS3_1, 10, df_within_CvLPS3)
        p_CvLPS3_2 = psturng (q_CvLPS3_2, 10, df_within_CvLPS3)
        p_CvLPS3_3 = psturng (q_CvLPS3_3, 10, df_within_CvLPS3)
        p_CvLPS3_4 = psturng (q_CvLPS3_4, 10, df_within_CvLPS3)
        p_CvLPS3_5 = psturng (q_CvLPS3_5, 10, df_within_CvLPS3)
        
        #################################################

        print ("========================================================")
        print ("")
        print ("Your Tukey's HSD Results Are:")
        print ("")

        if HSD_BBS_ten > HSD_control:
            print ("BBS vs. 10:  *  | p =", p_control_1)
        else:
            print ("BBS vs. 10:  ns | p =", p_control_1)

        if HSD_BBS_thirty > HSD_control:
            print ("BBS vs. 30:  *  | p =", p_control_2)
        else:
            print ("BBS vs. 30:  ns | p =", p_control_2)

        if HSD_BBS_hundred > HSD_control:
            print ("BBS vs. 100: *  | p =", p_control_3)
        else:
            print ("BBS vs. 100: ns | p =", p_control_3)

        if HSD_BBS_threehundred > HSD_control:
            print ("BBS vs. 300: *  | p =", p_control_4)
        else:
            print ("BBS vs. 300: ns | p =", p_control_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_BBS_ten_LPS > HSD_LPS:
            print ("1st Dose BBS LPS vs. 10 LPS:  *  | p =", p_LPS_1)
        else:
            print ("1st Dose BBS LPS vs. 10 LPS:  ns | p =", p_LPS_1)
            
        if HSD_BBS_thirty_LPS > HSD_LPS:
            print ("1st Dose BBS LPS vs. 30 LPS:  *  | p =", p_LPS_2)
        else:
            print ("1st Dose BBS LPS vs. 30 LPS:  ns | p =", p_LPS_2)

        if HSD_BBS_hundred_LPS > HSD_LPS:
            print ("1st Dose BBS LPS vs. 100 LPS: *  | p =", p_LPS_3)
        else:
            print ("1st Dose BBS LPS vs. 100 LPS: ns | p =", p_LPS_3)

        if HSD_BBS_threehundred_LPS > HSD_LPS:
            print ("1st Dose BBS LPS vs. 300 LPS: *  | p =", p_LPS_4)
        else:
            print ("1st Dose BBS LPS vs. 300 LPS: ns | p =", p_LPS_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_BBS_ten_LPS2 > HSD_LPS2:
            print ("2nd Dose BBS LPS vs. 10 LPS:  *  | p =", p_LPS2_1)
        else:
            print ("2nd Dose BBS LPS vs. 10 LPS:  ns | p =", p_LPS2_1)
            
        if HSD_BBS_thirty_LPS2 > HSD_LPS2:
            print ("2nd Dose BBS LPS vs. 30 LPS:  *  | p =", p_LPS2_2)
        else:
            print ("2nd Dose BBS LPS vs. 30 LPS:  ns | p =", p_LPS2_2)

        if HSD_BBS_hundred_LPS2 > HSD_LPS2:
            print ("2nd Dose BBS LPS vs. 100 LPS: *  | p =", p_LPS2_3)
        else:
            print ("2nd Dose BBS LPS vs. 100 LPS: ns | p =", p_LPS2_3)

        if HSD_BBS_threehundred_LPS2 > HSD_LPS2:
            print ("2nd Dose BBS LPS vs. 300 LPS: *  | p =", p_LPS2_4)
        else:
            print ("2nd Dose BBS LPS vs. 300 LPS: ns | p =", p_LPS2_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_BBS_ten_LPS3 > HSD_LPS3:
            print ("3rd Dose BBS LPS vs. 10 LPS:  *  | p =", p_LPS3_1)
        else:
            print ("3rd Dose BBS LPS vs. 10 LPS:  ns | p =", p_LPS3_1)
            
        if HSD_BBS_thirty_LPS3 > HSD_LPS3:
            print ("3rd Dose BBS LPS vs. 30 LPS:  *  | p =", p_LPS3_2)
        else:
            print ("3rd Dose BBS LPS vs. 30 LPS:  ns | p =", p_LPS3_2)

        if HSD_BBS_hundred_LPS3 > HSD_LPS3:
            print ("3rd Dose BBS LPS vs. 100 LPS: *  | p =", p_LPS3_3)
        else:
            print ("3rd Dose BBS LPS vs. 100 LPS: ns | p =", p_LPS3_3)

        if HSD_BBS_threehundred_LPS3 > HSD_LPS3:
            print ("3rd Dose BBS LPS vs. 300 LPS: *  | p =", p_LPS3_4)
        else:
            print ("3rd Dose BBS LPS vs. 300 LPS: ns | p =", p_LPS3_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_CvLPS_1 > HSD_CvLPS:
            print ("Control BBS vs. 1st Dose BBS LPS:  * | p =", p_CvLPS_1)
        else:
            print ("Control BBS vs. 1st Dose BBS LPS: ns | p =", p_CvLPS_1)

        if HSD_CvLPS_2 > HSD_CvLPS:
            print ("Control 10 vs. 1st Dose 10 LPS:    * | p =", p_CvLPS_2)
        else:
            print ("Control 10 vs. 1st Dose 10 LPS:   ns | p =", p_CvLPS_2)

        if HSD_CvLPS_3 > HSD_CvLPS:
            print ("Control 30 vs. 1st Dose 30 LPS:    * | p =", p_CvLPS_3)
        else:
            print ("Control 30 vs. 1st Dose 30 LPS:   ns | p =", p_CvLPS_3)

        if HSD_CvLPS_4 > HSD_CvLPS:
            print ("Control 100 vs. 1st Dose 100 LPS:  * | p =", p_CvLPS_4)
        else:
            print ("Control 100 vs. 1st Dose 100 LPS: ns | p =", p_CvLPS_4)

        if HSD_CvLPS_5 > HSD_CvLPS:
            print ("Control 300 vs. 1st Dose 300 LPS:  * | p =", p_CvLPS_5)
        else:
            print ("Control 300 vs. 1st Dose 300 LPS: ns | p =", p_CvLPS_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_CvLPS2_1 > HSD_CvLPS2:
            print ("Control BBS vs. 2nd Dose BBS LPS:  * | p =", p_CvLPS2_1)
        else:
            print ("Control BBS vs. 2nd Dose BBS LPS: ns | p =", p_CvLPS2_1)

        if HSD_CvLPS2_2 > HSD_CvLPS2:
            print ("Control 10 vs. 2nd Dose 10 LPS:    * | p =", p_CvLPS2_2)
        else:
            print ("Control 10 vs. 2nd Dose 10 LPS:   ns | p =", p_CvLPS2_2)

        if HSD_CvLPS2_3 > HSD_CvLPS2:
            print ("Control 30 vs. 2nd Dose 30 LPS:    * | p =", p_CvLPS2_3)
        else:
            print ("Control 30 vs. 2nd Dose 30 LPS:   ns | p =", p_CvLPS2_3)

        if HSD_CvLPS2_4 > HSD_CvLPS2:
            print ("Control 100 vs. 2nd Dose 100 LPS:  * | p =", p_CvLPS2_4)
        else:
            print ("Control 100 vs. 2nd Dose 100 LPS: ns | p =", p_CvLPS2_4)

        if HSD_CvLPS2_5 > HSD_CvLPS2:
            print ("Control 300 vs. 2nd Dose 300 LPS:  * | p =", p_CvLPS2_5)
        else:
            print ("Control 300 vs. 2nd Dose 300 LPS: ns | p =", p_CvLPS2_5)        

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_CvLPS3_1 > HSD_CvLPS3:
            print ("Control BBS vs. 3rd Dose BBS LPS:  * | p =", p_CvLPS3_1)
        else:
            print ("Control BBS vs. 3rd Dose BBS LPS: ns | p =", p_CvLPS3_1)

        if HSD_CvLPS3_2 > HSD_CvLPS3:
            print ("Control 10 vs. 3rd Dose 10 LPS:    * | p =", p_CvLPS3_2)
        else:
            print ("Control 10 vs. 3rd Dose 10 LPS:   ns | p =", p_CvLPS3_2)

        if HSD_CvLPS3_3 > HSD_CvLPS3:
            print ("Control 30 vs. 3rd Dose 30 LPS:    * | p =", p_CvLPS3_3)
        else:
            print ("Control 30 vs. 3rd Dose 30 LPS:   ns | p =", p_CvLPS3_3)

        if HSD_CvLPS3_4 > HSD_CvLPS3:
            print ("Control 100 vs. 3rd Dose 100 LPS:  * | p =", p_CvLPS3_4)
        else:
            print ("Control 100 vs. 3rd Dose 100 LPS: ns | p =", p_CvLPS3_4)

        if HSD_CvLPS3_5 > HSD_CvLPS3:
            print ("Control 300 vs. 3rd Dose 300 LPS:  * | p =", p_CvLPS3_5)
        else:
            print ("Control 300 vs. 3rd Dose 300 LPS: ns | p =", p_CvLPS3_5)   

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_LPSvLPS2_1 > HSD_LPSvLPS2:
            print ("1st Dose BBS LPS vs. 2nd Dose BBS LPS:  * | p =", p_LPSvLPS2_1)
        else:
            print ("1st Dose BBS LPS vs. 2nd Dose BBS LPS: ns | p =", p_LPSvLPS2_1)

        if HSD_LPSvLPS2_2 > HSD_LPSvLPS2:
            print ("1st Dose 10 LPS vs. 2nd Dose 10 LPS:    * | p =", p_LPSvLPS2_2)
        else:
            print ("1st Dose 10 LPS vs. 2nd Dose 10 LPS:   ns | p =", p_LPSvLPS2_2)

        if HSD_LPSvLPS2_3 > HSD_LPSvLPS2:
            print ("1st Dose 30 LPS vs. 2nd Dose 30 LPS:    * | p =", p_LPSvLPS2_3)
        else:
            print ("1st Dose 30 LPS vs. 2nd Dose 30 LPS:   ns | p =", p_LPSvLPS2_3)

        if HSD_LPSvLPS2_4 > HSD_LPSvLPS2:
            print ("1st Dose 100 LPS vs. 2nd Dose 100 LPS:  * | p =", p_LPSvLPS2_4)
        else:
            print ("1st Dose 100 LPS vs. 2nd Dose 100 LPS: ns | p =", p_LPSvLPS2_4)

        if HSD_LPSvLPS2_5 > HSD_LPSvLPS2:
            print ("1st Dose 300 LPS vs. 2nd Dose 300 LPS:  * | p =", p_LPSvLPS2_5)
        else:
            print ("1st Dose 300 LPS vs. 2nd Dose 300 LPS: ns | p =", p_LPSvLPS2_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_LPSvLPS3_1 > HSD_LPSvLPS3:
            print ("1st Dose BBS LPS vs. 3rd Dose BBS LPS:  * | p =", p_LPSvLPS3_1)
        else:
            print ("1st Dose BBS LPS vs. 3rd Dose BBS LPS: ns | p =", p_LPSvLPS3_1)

        if HSD_LPSvLPS3_2 > HSD_LPSvLPS3:
            print ("1st Dose 10 LPS vs. 3rd Dose 10 LPS:    * | p =", p_LPSvLPS3_2)
        else:
            print ("1st Dose 10 LPS vs. 3rd Dose 10 LPS:   ns | p =", p_LPSvLPS3_2)

        if HSD_LPSvLPS3_3 > HSD_LPSvLPS3:
            print ("1st Dose 30 LPS vs. 3rd Dose 30 LPS:    * | p =", p_LPSvLPS3_3)
        else:
            print ("1st Dose 30 LPS vs. 3rd Dose 30 LPS:   ns | p =", p_LPSvLPS3_3)

        if HSD_LPSvLPS3_4 > HSD_LPSvLPS3:
            print ("1st Dose 100 LPS vs. 3rd Dose 100 LPS:  * | p =", p_LPSvLPS3_4)
        else:
            print ("1st Dose 100 LPS vs. 3rd Dose 100 LPS: ns | p =", p_LPSvLPS3_4)

        if HSD_LPSvLPS3_5 > HSD_LPSvLPS3:
            print ("1st Dose 300 LPS vs. 3rd Dose 300 LPS:  * | p =", p_LPSvLPS3_5)
        else:
            print ("1st Dose 300 LPS vs. 3rd Dose 300 LPS: ns | p =", p_LPSvLPS3_5)        

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_LPS2vLPS3_1 > HSD_LPS2vLPS3:
            print ("2nd Dose BBS LPS vs. 3rd Dose BBS LPS:  * | p =", p_LPS2vLPS3_1)
        else:
            print ("2nd Dose BBS LPS vs. 3rd Dose BBS LPS: ns | p =", p_LPS2vLPS3_1)

        if HSD_LPS2vLPS3_2 > HSD_LPS2vLPS3:
            print ("2nd Dose 10 LPS vs. 3rd Dose 10 LPS:    * | p =", p_LPS2vLPS3_2)
        else:
            print ("2nd Dose 10 LPS vs. 3rd Dose 10 LPS:   ns | p =", p_LPS2vLPS3_2)

        if HSD_LPS2vLPS3_3 > HSD_LPS2vLPS3:
            print ("2nd Dose 30 LPS vs. 3rd Dose 30 LPS:    * | p =", p_LPS2vLPS3_3)
        else:
            print ("2nd Dose 30 LPS vs. 3rd Dose 30 LPS:   ns | p =", p_LPS2vLPS3_3)

        if HSD_LPS2vLPS3_4 > HSD_LPS2vLPS3:
            print ("2nd Dose 100 LPS vs. 3rd Dose 100 LPS:  * | p =", p_LPS2vLPS3_4)
        else:
            print ("2nd Dose 100 LPS vs. 3rd Dose 100 LPS: ns | p =", p_LPS2vLPS3_4)

        if HSD_LPS2vLPS3_5 > HSD_LPS2vLPS3:
            print ("2nd Dose 300 LPS vs. 3rd Dose 300 LPS:  * | p =", p_LPS2vLPS3_5)
        else:
            print ("2nd Dose 300 LPS vs. 3rd Dose 300 LPS: ns | p =", p_LPS2vLPS3_5)   

        print ("")
        print ("========================================================")
        print ("")

        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")
    
    elif answer == "3":

        #Tukey-Kramer: if the MD is greater than the CD, significant

        from statsmodels.stats.libqsturng import psturng
        
        #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/Q critical values.csv", header = None)
        crit_table = [37.082,	10.881,	7.502,	6.287,	5.673,	5.305,	5.06,	4.886,	4.755,	4.654,	4.574,	4.508,	4.453,	4.407,	4.367,	4.333,	4.303,	4.276,	4.253,	4.232,	4.213,	4.196,	4.18,	4.166,	4.153,	4.141,	4.13,	4.12,	4.111,	4.102,	4.094,	4.086,	4.079,	4.072,	4.066,	4.06,	4.054,	4.049,	4.044,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887]
        crit_table2 = [49.071,13.988,9.462,7.826,6.995,6.493,6.158,5.918,5.738,5.598,5.486,5.395,5.318,5.253,5.198,5.15,5.108,5.071,5.037,5.008,4.981,4.957,4.935,4.915,4.897,4.88,4.864,4.85,4.837,4.824,4.812,4.802,4.791,4.782,4.773,4.764,4.756,4.749,4.741,4.735,4.735,4.735,4.735,4.735,4.735,4.735,4.735,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.69,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.646,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.603,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.56,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517,4.517]

        #################################################

        #HSD_control = crit_table.iloc[int (float(df_within_control)-1), 0]
        #HSD_LPS = crit_table.iloc[int (float(df_within_LPS)-1), 0]

        HSD_control = crit_table[int (float(df_within_control)-1)]
        HSD_LPS = crit_table[int (float(df_within_LPS)-1)]
        HSD_LPS2 = crit_table[int (float(df_within_LPS2)-1)]
        HSD_LPS3 = crit_table[int (float(df_within_LPS3)-1)]
        HSD_LPSvLPS2 = crit_table2[int (float(df_within_LPSvLPS2)-1)]
        HSD_LPSvLPS3 = crit_table2[int (float(df_within_LPSvLPS3)-1)]
        HSD_LPS2vLPS3 = crit_table2[int (float(df_within_LPS2vLPS3)-1)]
        HSD_CvLPS = crit_table2[int (float(df_within_CvLPS)-1)]
        HSD_CvLPS2 = crit_table2[int (float(df_within_CvLPS2)-1)]
        HSD_CvLPS3 = crit_table2[int (float(df_within_CvLPS3)-1)]

        #################################################

        CD_control_1 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(ten)) ))
        CD_control_2 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(thirty)) ))
        CD_control_3 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(hundred)) ))
        CD_control_4 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(threehundred)) ))

        CD_LPS_1 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        CD_LPS_2 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        CD_LPS_3 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        CD_LPS_4 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        CD_LPS2_1 = HSD_LPS2 * math.sqrt ( (MSE_within_LPS2 / 2) * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
        CD_LPS2_2 = HSD_LPS2 * math.sqrt ( (MSE_within_LPS2 / 2) * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
        CD_LPS2_3 = HSD_LPS2 * math.sqrt ( (MSE_within_LPS2 / 2) * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
        CD_LPS2_4 = HSD_LPS2 * math.sqrt ( (MSE_within_LPS2 / 2) * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

        CD_LPS3_1 = HSD_LPS3 * math.sqrt ( (MSE_within_LPS3 / 2) * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
        CD_LPS3_2 = HSD_LPS3 * math.sqrt ( (MSE_within_LPS3 / 2) * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
        CD_LPS3_3 = HSD_LPS3 * math.sqrt ( (MSE_within_LPS3 / 2) * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
        CD_LPS3_4 = HSD_LPS3 * math.sqrt ( (MSE_within_LPS3 / 2) * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

        CD_LPSvLPS2_1 = HSD_LPSvLPS2 * math.sqrt ( (MSE_within_LPSvLPS2 / 2) * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS2)) ))
        CD_LPSvLPS2_2 = HSD_LPSvLPS2 * math.sqrt ( (MSE_within_LPSvLPS2 / 2) * ( (1/len(ten_LPS)) + (1/len(ten_LPS2)) ))
        CD_LPSvLPS2_3 = HSD_LPSvLPS2 * math.sqrt ( (MSE_within_LPSvLPS2 / 2) * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS2)) ))
        CD_LPSvLPS2_4 = HSD_LPSvLPS2 * math.sqrt ( (MSE_within_LPSvLPS2 / 2) * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS2)) ))
        CD_LPSvLPS2_5 = HSD_LPSvLPS2 * math.sqrt ( (MSE_within_LPSvLPS2 / 2) * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS2)) ))

        CD_LPSvLPS3_1 = HSD_LPSvLPS3 * math.sqrt ( (MSE_within_LPSvLPS3 / 2) * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS3)) ))
        CD_LPSvLPS3_2 = HSD_LPSvLPS3 * math.sqrt ( (MSE_within_LPSvLPS3 / 2) * ( (1/len(ten_LPS)) + (1/len(ten_LPS3)) ))
        CD_LPSvLPS3_3 = HSD_LPSvLPS3 * math.sqrt ( (MSE_within_LPSvLPS3 / 2) * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS3)) ))
        CD_LPSvLPS3_4 = HSD_LPSvLPS3 * math.sqrt ( (MSE_within_LPSvLPS3 / 2) * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS3)) ))
        CD_LPSvLPS3_5 = HSD_LPSvLPS3 * math.sqrt ( (MSE_within_LPSvLPS3 / 2) * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS3)) ))

        CD_LPS2vLPS3_1 = HSD_LPS2vLPS3 * math.sqrt ( (MSE_within_LPS2vLPS3 / 2) * ( (1/len(BBS_LPS2)) + (1/len(BBS_LPS3)) ))
        CD_LPS2vLPS3_2 = HSD_LPS2vLPS3 * math.sqrt ( (MSE_within_LPS2vLPS3 / 2) * ( (1/len(ten_LPS2)) + (1/len(ten_LPS3)) ))
        CD_LPS2vLPS3_3 = HSD_LPS2vLPS3 * math.sqrt ( (MSE_within_LPS2vLPS3 / 2) * ( (1/len(thirty_LPS2)) + (1/len(thirty_LPS3)) ))
        CD_LPS2vLPS3_4 = HSD_LPS2vLPS3 * math.sqrt ( (MSE_within_LPS2vLPS3 / 2) * ( (1/len(hundred_LPS2)) + (1/len(hundred_LPS3)) ))
        CD_LPS2vLPS3_5 = HSD_LPS2vLPS3 * math.sqrt ( (MSE_within_LPS2vLPS3 / 2) * ( (1/len(threehundred_LPS2)) + (1/len(threehundred_LPS3)) ))

        CD_CvLPS_1 = HSD_CvLPS * math.sqrt ( (MSE_within_CvLPS / 2) * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
        CD_CvLPS_2 = HSD_CvLPS * math.sqrt ( (MSE_within_CvLPS / 2) * ( (1/len(ten)) + (1/len(ten_LPS)) ))
        CD_CvLPS_3 = HSD_CvLPS * math.sqrt ( (MSE_within_CvLPS / 2) * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
        CD_CvLPS_4 = HSD_CvLPS * math.sqrt ( (MSE_within_CvLPS / 2) * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
        CD_CvLPS_5 = HSD_CvLPS * math.sqrt ( (MSE_within_CvLPS / 2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))

        CD_CvLPS2_1 = HSD_CvLPS2 * math.sqrt ( (MSE_within_CvLPS2 / 2) * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
        CD_CvLPS2_2 = HSD_CvLPS2 * math.sqrt ( (MSE_within_CvLPS2 / 2) * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
        CD_CvLPS2_3 = HSD_CvLPS2 * math.sqrt ( (MSE_within_CvLPS2 / 2) * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
        CD_CvLPS2_4 = HSD_CvLPS2 * math.sqrt ( (MSE_within_CvLPS2 / 2) * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
        CD_CvLPS2_5 = HSD_CvLPS2 * math.sqrt ( (MSE_within_CvLPS2 / 2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))

        CD_CvLPS3_1 = HSD_CvLPS3 * math.sqrt ( (MSE_within_CvLPS3 / 2) * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))
        CD_CvLPS3_2 = HSD_CvLPS3 * math.sqrt ( (MSE_within_CvLPS3 / 2) * ( (1/len(ten)) + (1/len(ten_LPS3)) ))
        CD_CvLPS3_3 = HSD_CvLPS3 * math.sqrt ( (MSE_within_CvLPS3 / 2) * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))
        CD_CvLPS3_4 = HSD_CvLPS3 * math.sqrt ( (MSE_within_CvLPS3 / 2) * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))
        CD_CvLPS3_5 = HSD_CvLPS3 * math.sqrt ( (MSE_within_CvLPS3 / 2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

        #################################################

        MD_control_1 = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        MD_control_2 = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        MD_control_3 = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        MD_control_4 = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))

        MD_LPS_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_LPS_2 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_LPS_3 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_LPS_4 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        MD_LPS2_1 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_LPS2_2 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_LPS2_3 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_LPS2_4 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_LPS3_1 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPS3_2 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPS3_3 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPS3_4 = abs (((sum (BBS_LPS3) / len (BBS_LPS3)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_LPSvLPS2_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        MD_LPSvLPS2_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_LPSvLPS2_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_LPSvLPS2_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_LPSvLPS2_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_LPSvLPS3_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_LPSvLPS3_2 = abs (((sum (ten_LPS) / len (ten_LPS)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPSvLPS3_3 = abs (((sum (thirty_LPS) / len (thirty_LPS)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPSvLPS3_4 = abs (((sum (hundred_LPS) / len (hundred_LPS)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPSvLPS3_5 = abs (((sum (threehundred_LPS) / len (threehundred_LPS)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_LPS2vLPS3_1 = abs (((sum (BBS_LPS2) / len (BBS_LPS2)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_LPS2vLPS3_2 = abs (((sum (ten_LPS2) / len (ten_LPS2)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_LPS2vLPS3_3 = abs (((sum (thirty_LPS2) / len (thirty_LPS2)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_LPS2vLPS3_4 = abs (((sum (hundred_LPS2) / len (hundred_LPS2)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_LPS2vLPS3_5 = abs (((sum (threehundred_LPS2) / len (threehundred_LPS2)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        MD_CvLPS_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS) / len (BBS_LPS))))
        MD_CvLPS_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_CvLPS_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_CvLPS_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_CvLPS_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        MD_CvLPS2_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS2) / len (BBS_LPS2))))
        MD_CvLPS2_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS2) / len (ten_LPS2))))
        MD_CvLPS2_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS2) / len (thirty_LPS2))))
        MD_CvLPS2_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS2) / len (hundred_LPS2))))
        MD_CvLPS2_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS2) / len (threehundred_LPS2))))

        MD_CvLPS3_1 = abs (((sum (BBS) / len (BBS)) - (sum (BBS_LPS3) / len (BBS_LPS3))))
        MD_CvLPS3_2 = abs (((sum (ten) / len (ten)) - (sum (ten_LPS3) / len (ten_LPS3))))
        MD_CvLPS3_3 = abs (((sum (thirty) / len (thirty)) - (sum (thirty_LPS3) / len (thirty_LPS3))))
        MD_CvLPS3_4 = abs (((sum (hundred) / len (hundred)) - (sum (hundred_LPS3) / len (hundred_LPS3))))
        MD_CvLPS3_5 = abs (((sum (threehundred) / len (threehundred)) - (sum (threehundred_LPS3) / len (threehundred_LPS3))))

        #################################################

        SED_control_1 = math.sqrt ((MSE_within_control/2) * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt ((MSE_within_control/2) * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt ((MSE_within_control/2) * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt ((MSE_within_control/2) * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt ((MSE_within_LPS/2) * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt ((MSE_within_LPS/2) * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt ((MSE_within_LPS/2) * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt ((MSE_within_LPS/2) * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        SED_LPS2_1 = math.sqrt ((MSE_within_LPS2/2) * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
        SED_LPS2_2 = math.sqrt ((MSE_within_LPS2/2) * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
        SED_LPS2_3 = math.sqrt ((MSE_within_LPS2/2) * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
        SED_LPS2_4 = math.sqrt ((MSE_within_LPS2/2) * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

        SED_LPS3_1 = math.sqrt ((MSE_within_LPS3/2) * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
        SED_LPS3_2 = math.sqrt ((MSE_within_LPS3/2) * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
        SED_LPS3_3 = math.sqrt ((MSE_within_LPS3/2) * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
        SED_LPS3_4 = math.sqrt ((MSE_within_LPS3/2) * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

        SED_LPSvLPS2_1 = math.sqrt ((MSE_within_LPSvLPS2/2) * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS2)) ))
        SED_LPSvLPS2_2 = math.sqrt ((MSE_within_LPSvLPS2/2) * ( (1/len(ten_LPS)) + (1/len(ten_LPS2)) ))
        SED_LPSvLPS2_3 = math.sqrt ((MSE_within_LPSvLPS2/2) * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS2)) ))
        SED_LPSvLPS2_4 = math.sqrt ((MSE_within_LPSvLPS2/2) * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS2)) ))
        SED_LPSvLPS2_5 = math.sqrt ((MSE_within_LPSvLPS2/2) * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS2)) ))

        SED_LPSvLPS3_1 = math.sqrt ((MSE_within_LPSvLPS3/2) * ( (1/len(BBS_LPS)) + (1/len(BBS_LPS3)) ))
        SED_LPSvLPS3_2 = math.sqrt ((MSE_within_LPSvLPS3/2) * ( (1/len(ten_LPS)) + (1/len(ten_LPS3)) ))
        SED_LPSvLPS3_3 = math.sqrt ((MSE_within_LPSvLPS3/2) * ( (1/len(thirty_LPS)) + (1/len(thirty_LPS3)) ))
        SED_LPSvLPS3_4 = math.sqrt ((MSE_within_LPSvLPS3/2) * ( (1/len(hundred_LPS)) + (1/len(hundred_LPS3)) ))
        SED_LPSvLPS3_5 = math.sqrt ((MSE_within_LPSvLPS3/2) * ( (1/len(threehundred_LPS)) + (1/len(threehundred_LPS3)) ))
        
        SED_LPS2vLPS3_1 = math.sqrt ((MSE_within_LPS2vLPS3/2) * ( (1/len(BBS_LPS2)) + (1/len(BBS_LPS3)) ))
        SED_LPS2vLPS3_2 = math.sqrt ((MSE_within_LPS2vLPS3/2) * ( (1/len(ten_LPS2)) + (1/len(ten_LPS3)) ))
        SED_LPS2vLPS3_3 = math.sqrt ((MSE_within_LPS2vLPS3/2) * ( (1/len(thirty_LPS2)) + (1/len(thirty_LPS3)) ))
        SED_LPS2vLPS3_4 = math.sqrt ((MSE_within_LPS2vLPS3/2) * ( (1/len(hundred_LPS2)) + (1/len(hundred_LPS3)) ))
        SED_LPS2vLPS3_5 = math.sqrt ((MSE_within_LPS2vLPS3/2) * ( (1/len(threehundred_LPS2)) + (1/len(threehundred_LPS3)) ))

        SED_CvLPS_1 = math.sqrt ((MSE_within_CvLPS/2) * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
        SED_CvLPS_2 = math.sqrt ((MSE_within_CvLPS/2) * ( (1/len(ten)) + (1/len(ten_LPS)) ))
        SED_CvLPS_3 = math.sqrt ((MSE_within_CvLPS/2) * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
        SED_CvLPS_4 = math.sqrt ((MSE_within_CvLPS/2) * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
        SED_CvLPS_5 = math.sqrt ((MSE_within_CvLPS/2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))

        SED_CvLPS2_1 = math.sqrt ((MSE_within_CvLPS2/2) * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
        SED_CvLPS2_2 = math.sqrt ((MSE_within_CvLPS2/2) * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
        SED_CvLPS2_3 = math.sqrt ((MSE_within_CvLPS2/2) * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
        SED_CvLPS2_4 = math.sqrt ((MSE_within_CvLPS2/2) * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
        SED_CvLPS2_5 = math.sqrt ((MSE_within_CvLPS2/2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))

        SED_CvLPS3_1 = math.sqrt ((MSE_within_CvLPS3/2) * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))
        SED_CvLPS3_2 = math.sqrt ((MSE_within_CvLPS3/2) * ( (1/len(ten)) + (1/len(ten_LPS3)) ))
        SED_CvLPS3_3 = math.sqrt ((MSE_within_CvLPS3/2) * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))
        SED_CvLPS3_4 = math.sqrt ((MSE_within_CvLPS3/2) * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))
        SED_CvLPS3_5 = math.sqrt ((MSE_within_CvLPS3/2) * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

        #################################################

        q_control_1 = MD_control_1 / (SED_control_1 * 1/math.sqrt(2))
        q_control_2 = MD_control_2 / (SED_control_2 * 1/math.sqrt(2))
        q_control_3 = MD_control_3 / (SED_control_3 * 1/math.sqrt(2))
        q_control_4 = MD_control_4 / (SED_control_4 * 1/math.sqrt(2))

        q_LPS_1 = MD_LPS_1 / (SED_LPS_1 * 1/math.sqrt(2))
        q_LPS_2 = MD_LPS_2 / (SED_LPS_2 * 1/math.sqrt(2))
        q_LPS_3 = MD_LPS_3 / (SED_LPS_3 * 1/math.sqrt(2))
        q_LPS_4 = MD_LPS_4 / (SED_LPS_4 * 1/math.sqrt(2))

        q_LPS2_1 = MD_LPS2_1 / (SED_LPS2_1 * 1/math.sqrt(2))
        q_LPS2_2 = MD_LPS2_2 / (SED_LPS2_2 * 1/math.sqrt(2))
        q_LPS2_3 = MD_LPS2_3 / (SED_LPS2_3 * 1/math.sqrt(2))
        q_LPS2_4 = MD_LPS2_4 / (SED_LPS2_4 * 1/math.sqrt(2))

        q_LPS3_1 = MD_LPS3_1 / (SED_LPS3_1 * 1/math.sqrt(2))
        q_LPS3_2 = MD_LPS3_2 / (SED_LPS3_2 * 1/math.sqrt(2))
        q_LPS3_3 = MD_LPS3_3 / (SED_LPS3_3 * 1/math.sqrt(2))
        q_LPS3_4 = MD_LPS3_4 / (SED_LPS3_4 * 1/math.sqrt(2))

        q_LPSvLPS2_1 = MD_LPSvLPS2_1 / (SED_LPSvLPS2_1 * 1/math.sqrt(2))
        q_LPSvLPS2_2 = MD_LPSvLPS2_2 / (SED_LPSvLPS2_2 * 1/math.sqrt(2))
        q_LPSvLPS2_3 = MD_LPSvLPS2_3 / (SED_LPSvLPS2_3 * 1/math.sqrt(2))
        q_LPSvLPS2_4 = MD_LPSvLPS2_4 / (SED_LPSvLPS2_4 * 1/math.sqrt(2))
        q_LPSvLPS2_5 = MD_LPSvLPS2_5 / (SED_LPSvLPS2_5 * 1/math.sqrt(2))

        q_LPSvLPS3_1 = MD_LPSvLPS3_1 / (SED_LPSvLPS3_1 * 1/math.sqrt(2))
        q_LPSvLPS3_2 = MD_LPSvLPS3_2 / (SED_LPSvLPS3_2 * 1/math.sqrt(2))
        q_LPSvLPS3_3 = MD_LPSvLPS3_3 / (SED_LPSvLPS3_3 * 1/math.sqrt(2))
        q_LPSvLPS3_4 = MD_LPSvLPS3_4 / (SED_LPSvLPS3_4 * 1/math.sqrt(2))
        q_LPSvLPS3_5 = MD_LPSvLPS3_5 / (SED_LPSvLPS3_5 * 1/math.sqrt(2))

        q_LPS2vLPS3_1 = MD_LPS2vLPS3_1 / (SED_LPS2vLPS3_1 * 1/math.sqrt(2))
        q_LPS2vLPS3_2 = MD_LPS2vLPS3_2 / (SED_LPS2vLPS3_2 * 1/math.sqrt(2))
        q_LPS2vLPS3_3 = MD_LPS2vLPS3_3 / (SED_LPS2vLPS3_3 * 1/math.sqrt(2))
        q_LPS2vLPS3_4 = MD_LPS2vLPS3_4 / (SED_LPS2vLPS3_4 * 1/math.sqrt(2))
        q_LPS2vLPS3_5 = MD_LPS2vLPS3_5 / (SED_LPS2vLPS3_5 * 1/math.sqrt(2))

        q_CvLPS_1 = MD_CvLPS_1 / (SED_CvLPS_1 * 1/math.sqrt(2))
        q_CvLPS_2 = MD_CvLPS_2 / (SED_CvLPS_2 * 1/math.sqrt(2))
        q_CvLPS_3 = MD_CvLPS_3 / (SED_CvLPS_3 * 1/math.sqrt(2))
        q_CvLPS_4 = MD_CvLPS_4 / (SED_CvLPS_4 * 1/math.sqrt(2))
        q_CvLPS_5 = MD_CvLPS_5 / (SED_CvLPS_5 * 1/math.sqrt(2))

        q_CvLPS2_1 = MD_CvLPS2_1 / (SED_CvLPS2_1 * 1/math.sqrt(2))
        q_CvLPS2_2 = MD_CvLPS2_2 / (SED_CvLPS2_2 * 1/math.sqrt(2))
        q_CvLPS2_3 = MD_CvLPS2_3 / (SED_CvLPS2_3 * 1/math.sqrt(2))
        q_CvLPS2_4 = MD_CvLPS2_4 / (SED_CvLPS2_4 * 1/math.sqrt(2))
        q_CvLPS2_5 = MD_CvLPS2_5 / (SED_CvLPS2_5 * 1/math.sqrt(2))

        q_CvLPS3_1 = MD_CvLPS3_1 / (SED_CvLPS3_1 * 1/math.sqrt(2))
        q_CvLPS3_2 = MD_CvLPS3_2 / (SED_CvLPS3_2 * 1/math.sqrt(2))
        q_CvLPS3_3 = MD_CvLPS3_3 / (SED_CvLPS3_3 * 1/math.sqrt(2))
        q_CvLPS3_4 = MD_CvLPS3_4 / (SED_CvLPS3_4 * 1/math.sqrt(2))
        q_CvLPS3_5 = MD_CvLPS3_5 / (SED_CvLPS3_5 * 1/math.sqrt(2))
        
        #################################################

        p_control_1 = psturng (q_control_1, 5, df_within_control)
        p_control_2 = psturng (q_control_2, 5, df_within_control)
        p_control_3 = psturng (q_control_3, 5, df_within_control)
        p_control_4 = psturng (q_control_4, 5, df_within_control)
        
        p_LPS_1 = psturng (q_LPS_1, 5, df_within_LPS)
        p_LPS_2 = psturng (q_LPS_2, 5, df_within_LPS)
        p_LPS_3 = psturng (q_LPS_3, 5, df_within_LPS)
        p_LPS_4 = psturng (q_LPS_4, 5, df_within_LPS)

        p_LPS2_1 = psturng (q_LPS2_1, 5, df_within_LPS2)
        p_LPS2_2 = psturng (q_LPS2_2, 5, df_within_LPS2)
        p_LPS2_3 = psturng (q_LPS2_3, 5, df_within_LPS2)
        p_LPS2_4 = psturng (q_LPS2_4, 5, df_within_LPS2)

        p_LPS3_1 = psturng (q_LPS3_1, 5, df_within_LPS3)
        p_LPS3_2 = psturng (q_LPS3_2, 5, df_within_LPS3)
        p_LPS3_3 = psturng (q_LPS3_3, 5, df_within_LPS3)
        p_LPS3_4 = psturng (q_LPS3_4, 5, df_within_LPS3)

        p_LPSvLPS2_1 = psturng (q_LPSvLPS2_1, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_2 = psturng (q_LPSvLPS2_2, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_3 = psturng (q_LPSvLPS2_3, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_4 = psturng (q_LPSvLPS2_4, 10, df_within_LPSvLPS2)
        p_LPSvLPS2_5 = psturng (q_LPSvLPS2_5, 10, df_within_LPSvLPS2)

        p_LPSvLPS3_1 = psturng (q_LPSvLPS3_1, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_2 = psturng (q_LPSvLPS3_2, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_3 = psturng (q_LPSvLPS3_3, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_4 = psturng (q_LPSvLPS3_4, 10, df_within_LPSvLPS3)
        p_LPSvLPS3_5 = psturng (q_LPSvLPS3_5, 10, df_within_LPSvLPS3)

        p_LPS2vLPS3_1 = psturng (q_LPS2vLPS3_1, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_2 = psturng (q_LPS2vLPS3_2, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_3 = psturng (q_LPS2vLPS3_3, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_4 = psturng (q_LPS2vLPS3_4, 10, df_within_LPS2vLPS3)
        p_LPS2vLPS3_5 = psturng (q_LPS2vLPS3_5, 10, df_within_LPS2vLPS3)

        p_CvLPS_1 = psturng (q_CvLPS_1, 10, df_within_CvLPS)
        p_CvLPS_2 = psturng (q_CvLPS_2, 10, df_within_CvLPS)
        p_CvLPS_3 = psturng (q_CvLPS_3, 10, df_within_CvLPS)
        p_CvLPS_4 = psturng (q_CvLPS_4, 10, df_within_CvLPS)
        p_CvLPS_5 = psturng (q_CvLPS_5, 10, df_within_CvLPS)

        p_CvLPS2_1 = psturng (q_CvLPS2_1, 10, df_within_CvLPS2)
        p_CvLPS2_2 = psturng (q_CvLPS2_2, 10, df_within_CvLPS2)
        p_CvLPS2_3 = psturng (q_CvLPS2_3, 10, df_within_CvLPS2)
        p_CvLPS2_4 = psturng (q_CvLPS2_4, 10, df_within_CvLPS2)
        p_CvLPS2_5 = psturng (q_CvLPS2_5, 10, df_within_CvLPS2)

        p_CvLPS3_1 = psturng (q_CvLPS3_1, 10, df_within_CvLPS3)
        p_CvLPS3_2 = psturng (q_CvLPS3_2, 10, df_within_CvLPS3)
        p_CvLPS3_3 = psturng (q_CvLPS3_3, 10, df_within_CvLPS3)
        p_CvLPS3_4 = psturng (q_CvLPS3_4, 10, df_within_CvLPS3)
        p_CvLPS3_5 = psturng (q_CvLPS3_5, 10, df_within_CvLPS3)
        
        #################################################

        print ("========================================================")
        print ("")
        print ("Your Tukey-Kramer Results Are:")
        print ("")

        if MD_control_1 > CD_control_1:
            print ("BBS vs. 10:  *  | p=") #, p_control_1)
        else:
            print ("BBS vs. 10:  ns | p=")#, p_control_1)

        if MD_control_2 > CD_control_2:
            print ("BBS vs. 30:  *  | p=")#, p_control_2)
        else:
            print ("BBS vs. 30:  ns | p=")#, p_control_2)

        if MD_control_3 > CD_control_3:
            print ("BBS vs. 100: *  | p=")#, p_control_3)
        else:
            print ("BBS vs. 100: ns | p=")#, p_control_3)

        if MD_control_4 > CD_control_4:
            print ("BBS vs. 300: *  | p=")#, p_control_4)
        else:
            print ("BBS vs. 300: ns | p=")#, p_control_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS_1 > CD_LPS_1:
            print ("BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS_1)
        else:
            print ("BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS_1)
            
        if MD_LPS_2 > CD_LPS_2:
            print ("BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS_2)
        else:
            print ("BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS_2)

        if MD_LPS_3 > CD_LPS_3:
            print ("BBS LPS vs. 100 LPS: *  | p=")#, p_LPS_3)
        else:
            print ("BBS LPS vs. 100 LPS: ns | p=")#, p_LPS_3)

        if MD_LPS_4 > CD_LPS_4:
            print ("BBS LPS vs. 300 LPS: *  | p=")#, p_LPS_4)
        else:
            print ("BBS LPS vs. 300 LPS: ns | p=")#, p_LPS_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS2_1 > CD_LPS2_1:
            print ("2nd Dose BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS2_1)
        else:
            print ("2nd Dose BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS2_1)
            
        if MD_LPS2_2 > CD_LPS2_2:
            print ("2nd Dose BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS2_2)
        else:
            print ("2nd Dose BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS2_2)

        if MD_LPS2_3 > CD_LPS2_3:
            print ("2nd Dose BBS LPS vs. 100 LPS: *  | p=")#, p_LPS2_3)
        else:
            print ("2nd Dose BBS LPS vs. 100 LPS: ns | p=")#, p_LPS2_3)

        if MD_LPS2_4 > CD_LPS2_4:
            print ("2nd Dose BBS LPS vs. 300 LPS: *  | p=")#, p_LPS2_4)
        else:
            print ("2nd Dose BBS LPS vs. 300 LPS: ns | p=")#, p_LPS2_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS3_1 > CD_LPS3_1:
            print ("3rd Dose BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS3_1)
        else:
            print ("3rd Dose BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS3_1)
            
        if MD_LPS3_2 > CD_LPS3_2:
            print ("3rd Dose BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS3_2)
        else:
            print ("3rd Dose BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS3_2)

        if MD_LPS3_3 > CD_LPS3_3:
            print ("3rd Dose BBS LPS vs. 100 LPS: *  | p=")#, p_LPS3_3)
        else:
            print ("3rd Dose BBS LPS vs. 100 LPS: ns | p=")#, p_LPS3_3)

        if MD_LPS3_4 > CD_LPS3_4:
            print ("3rd Dose BBS LPS vs. 300 LPS: *  | p=")#, p_LPS3_4)
        else:
            print ("3rd Dose BBS LPS vs. 300 LPS: ns | p=")#, p_LPS3_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_CvLPS_1 > CD_CvLPS_1:
            print ("Control BBS vs. 1st Dose LPS BBS: *  | p=")#, p_CvLPS_1)
        else:
            print ("Control BBS vs. 1st Dose LPS BBS: ns | p=")#, p_CvLPS_1)

        if MD_CvLPS_2 > CD_CvLPS_2:
            print ("Control 10 vs. 1st Dose LPS 10:   *  | p=")#, p_CvLPS_2)
        else:
            print ("Control 10 vs. 1st Dose LPS 10:   ns | p=")#, p_CvLPS_2)

        if MD_CvLPS_3 > CD_CvLPS_3:
            print ("Control 30 vs. 1st Dose LPS 30:   *  | p=")#, p_CvLPS_3)
        else:
            print ("Control 30 vs. 1st Dose LPS 30:   ns | p=")#, p_CvLPS_3)

        if MD_CvLPS_4 > CD_CvLPS_4:
            print ("Control 100 vs. 1st Dose LPS 100: *  | p=")#, p_CvLPS_4)
        else:
            print ("Control 100 vs. 1st Dose LPS 100: ns | p=")#, p_CvLPS_4)

        if MD_CvLPS_5 > CD_CvLPS_5:
            print ("Control 300 vs. 1st Dose LPS 300: *  | p=")#, p_CvLPS_5)
        else:
            print ("Control 300 vs. 1st Dose LPS 300: ns | p=")#, p_CvLPS_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_CvLPS2_1 > CD_CvLPS2_1:
            print ("Control BBS vs. 2nd Dose LPS BBS: *  | p=")#, p_CvLPS2_1)
        else:
            print ("Control BBS vs. 2nd Dose LPS BBS: ns | p=")#, p_CvLPS2_1)

        if MD_CvLPS2_2 > CD_CvLPS2_2:
            print ("Control 10 vs. 2nd Dose LPS 10:   *  | p=")#, p_CvLPS2_2)
        else:
            print ("Control 10 vs. 2nd Dose LPS 10:   ns | p=")#, p_CvLPS2_2)

        if MD_CvLPS2_3 > CD_CvLPS2_3:
            print ("Control 30 vs. 2nd Dose LPS 30:   *  | p=")#, p_CvLPS2_3)
        else:
            print ("Control 30 vs. 2nd Dose LPS 30:   ns | p=")#, p_CvLPS2_3)

        if MD_CvLPS2_4 > CD_CvLPS2_4:
            print ("Control 100 vs. 2nd Dose LPS 100: *  | p=")#, p_CvLPS2_4)
        else:
            print ("Control 100 vs. 2nd Dose LPS 100: ns | p=")#, p_CvLPS2_4)

        if MD_CvLPS2_5 > CD_CvLPS2_5:
            print ("Control 300 vs. 2nd Dose LPS 300: *  | p=")#, p_CvLPS2_5)
        else:
            print ("Control 300 vs. 2nd Dose LPS 300: ns | p=")#, p_CvLPS2_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_CvLPS3_1 > CD_CvLPS3_1:
            print ("Control BBS vs. 3rd Dose LPS BBS: *  | p=")#, p_CvLPS3_1)
        else:
            print ("Control BBS vs. 3rd Dose LPS BBS: ns | p=")#, p_CvLPS3_1)

        if MD_CvLPS3_2 > CD_CvLPS3_2:
            print ("Control 10 vs. 3rd Dose LPS 10:   *  | p=")#, p_CvLPS3_2)
        else:
            print ("Control 10 vs. 3rd Dose LPS 10:   ns | p=")#, p_CvLPS3_2)

        if MD_CvLPS3_3 > CD_CvLPS3_3:
            print ("Control 30 vs. 3rd Dose LPS 30:   *  | p=")#, p_CvLPS3_3)
        else:
            print ("Control 30 vs. 3rd Dose LPS 30:   ns | p=")#, p_CvLPS3_3)

        if MD_CvLPS3_4 > CD_CvLPS3_4:
            print ("Control 100 vs. 3rd Dose LPS 100: *  | p=")#, p_CvLPS3_4)
        else:
            print ("Control 100 vs. 3rd Dose LPS 100: ns | p=")#, p_CvLPS3_4)

        if MD_CvLPS3_5 > CD_CvLPS3_5:
            print ("Control 300 vs. 3rd Dose LPS 300: *  | p=")#, p_CvLPS3_5)
        else:
            print ("Control 300 vs. 3rd Dose LPS 300: ns | p=")#, p_CvLPS3_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPSvLPS2_1 > CD_LPSvLPS2_1:
            print ("1st Dose LPS BBS vs. 2nd Dose LPS BBS: *  | p=")#, p_LPSvLPS2_1)
        else:
            print ("1st Dose LPS BBS vs. 2nd Dose LPS BBS: ns | p=")#, p_LPSvLPS2_1)

        if MD_LPSvLPS2_2 > CD_LPSvLPS2_2:
            print ("1st Dose LPS 10 vs. 2nd Dose LPS 10:   *  | p=")#, p_LPSvLPS2_2)
        else:
            print ("1st Dose LPS 10 vs. 2nd Dose LPS 10:   ns | p=")#, p_LPSvLPS2_2)

        if MD_LPSvLPS2_3 > CD_LPSvLPS2_3:
            print ("1st Dose LPS 30 vs. 2nd Dose LPS 30:   *  | p=")#, p_LPSvLPS2_3)
        else:
            print ("1st Dose LPS 30 vs. 2nd Dose LPS 30:   ns | p=")#, p_LPSvLPS2_3)

        if MD_LPSvLPS2_4 > CD_LPSvLPS2_4:
            print ("1st Dose LPS 100 vs. 2nd Dose LPS 100: *  | p=")#, p_LPSvLPS2_4)
        else:
            print ("1st Dose LPS 100 vs. 2nd Dose LPS 100: ns | p=")#, p_LPSvLPS2_4)

        if MD_LPSvLPS2_5 > CD_LPSvLPS2_5:
            print ("1st Dose LPS 300 vs. 2nd Dose LPS 300: *  | p=")#, p_LPSvLPS2_5)
        else:
            print ("1st Dose LPS 300 vs. 2nd Dose LPS 300: ns | p=")#, p_LPSvLPS2_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPSvLPS3_1 > CD_LPSvLPS3_1:
            print ("1st Dose LPS BBS vs. 3rd Dose LPS BBS: *  | p=")#, p_LPSvLPS3_1)
        else:
            print ("1st Dose LPS BBS vs. 3rd Dose LPS BBS: ns | p=")#, p_LPSvLPS3_1)

        if MD_LPSvLPS3_2 > CD_LPSvLPS3_2:
            print ("1st Dose LPS 10 vs. 3rd Dose LPS 10:   *  | p=")#, p_LPSvLPS3_2)
        else:
            print ("1st Dose LPS 10 vs. 3rd Dose LPS 10:   ns | p=")#, p_LPSvLPS3_2)

        if MD_LPSvLPS3_3 > CD_LPSvLPS3_3:
            print ("1st Dose LPS 30 vs. 3rd Dose LPS 30:   *  | p=")#, p_LPSvLPS3_3)
        else:
            print ("1st Dose LPS 30 vs. 3rd Dose LPS 30:   ns | p=")#, p_LPSvLPS3_3)

        if MD_LPSvLPS3_4 > CD_LPSvLPS3_4:
            print ("1st Dose LPS 100 vs. 3rd Dose LPS 100: *  | p=")#, p_LPSvLPS3_4)
        else:
            print ("1st Dose LPS 100 vs. 3rd Dose LPS 100: ns | p=")#, p_LPSvLPS3_4)

        if MD_LPSvLPS3_5 > CD_LPSvLPS3_5:
            print ("1st Dose LPS 300 vs. 3rd Dose LPS 300: *  | p=")#, p_LPSvLPS3_5)
        else:
            print ("1st Dose LPS 300 vs. 3rd Dose LPS 300: ns | p=")#, p_LPSvLPS3_5)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS2vLPS3_1 > CD_LPS2vLPS3_1:
            print ("2nd Dose LPS BBS vs. 3rd Dose LPS BBS: *  | p=")#, p_LPS2vLPS3_1)
        else:
            print ("2nd Dose LPS BBS vs. 3rd Dose LPS BBS: ns | p=")#, p_LPS2vLPS3_1)

        if MD_LPS2vLPS3_2 > CD_LPS2vLPS3_2:
            print ("2nd Dose LPS 10 vs. 3rd Dose LPS 10:   *  | p=")#, p_LPS2vLPS3_2)
        else:
            print ("2nd Dose LPS 10 vs. 3rd Dose LPS 10:   ns | p=")#, p_LPS2vLPS3_2)

        if MD_LPS2vLPS3_3 > CD_LPS2vLPS3_3:
            print ("2nd Dose LPS 30 vs. 3rd Dose LPS 30:   *  | p=")#, p_LPS2vLPS3_3)
        else:
            print ("2nd Dose LPS 30 vs. 3rd Dose LPS 30:   ns | p=")#, p_LPS2vLPS3_3)

        if MD_LPS2vLPS3_4 > CD_LPS2vLPS3_4:
            print ("2nd Dose LPS 100 vs. 3rd Dose LPS 100: *  | p=")#, p_LPS2vLPS3_4)
        else:
            print ("2nd Dose LPS 100 vs. 3rd Dose LPS 100: ns | p=")#, p_LPS2vLPS3_4)

        if MD_LPS2vLPS3_5 > CD_LPS2vLPS3_5:
            print ("2nd Dose LPS 300 vs. 3rd Dose LPS 300: *  | p=")#, p_LPS2vLPS3_5)
        else:
            print ("2nd Dose LPS 300 vs. 3rd Dose LPS 300: ns | p=")#, p_LPS2vLPS3_5)

        print ("")
        print ("*** p values are coming in a later version ***")

        print ("")
        print ("========================================================")
        print ("")
        
        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")


    if answer == "4":

        #k = 5, df = index up to df = 100, alpha level = 0.05; TABLE FOR VACCAE COMPARISONS
        crit_table = [21.85, 6.513, 4.538, 3.832, 3.476, 3.263, 3.123, 3.023, 2.948, 2.890, 2.845, 2.807, 2.776, 2.75, 2.727, 2.708, 2.691, 2.676, 2.663, 2.651, 2.64, 2.631, 2.622, 2.614, 2.607, 2.6, 2.594, 2.588, 2.583, 2.578, 2.574, 2.569, 2.565, 2.561, 2.558, 2.555, 2.551, 2.548, 2.546, 2.543, 2.543, 2.543, 2.543, 2.543, 2.543, 2.543, 2.543, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.526, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.508, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491, 2.491]

        #k = 4, df = index up to df = 100, alpha level = 0.05; TABLE FOR LPS COMPARISONS
        crit_table_LPS = [20.034,6.065,4.263,3.618,3.293,3.099,2.971,2.880,2.812,2.759,2.717,2.683,2.654,2.631,2.610,2.592,2.577,2.563,2.551,2.540,2.531,2.522,2.514,2.507,2.5,2.494,2.488,2.483,2.478,2.474,2.469,2.466,2.462,2.458,2.455,2.452,2.449,2.447,2.444,2.441,2.411,2.411,2.411,2.411,2.411,2.411,2.411,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.426,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.410,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394,2.394]

        d_crit_control = crit_table[df_within_control]
        d_crit_LPS = crit_table[df_within_LPS]
        d_crit_LPS2 = crit_table[df_within_LPS2]
        d_crit_LPS3 = crit_table[df_within_LPS3]
        d_crit_BBS_by_LPS = crit_table_LPS[df_within_BBS_by_LPS] 
        d_crit_ten_by_LPS = crit_table_LPS[df_within_ten_by_LPS] 
        d_crit_thirty_by_LPS = crit_table_LPS[df_within_thirty_by_LPS] 
        d_crit_hundred_by_LPS = crit_table_LPS[df_within_hundred_by_LPS] 
        d_crit_threehundred_by_LPS = crit_table_LPS[df_within_threehundred_by_LPS] 
        d_crit_Mt_by_LPS = crit_table_LPS[df_within_Mt_by_LPS] 

        try:
            d_crit_for_vaccae = crit_table[Pooled_df_within_all_LPS_for_vaccae]
        except:
            d_crit_for_vaccae = crit_table[-1]

        try:
            d_crit_for_LPS = crit_table_LPS[Pooled_df_within_all_LPS_for_LPS]
        except:
            d_crit_for_LPS = crit_table_LPS[-1]

        try:
            d_crit_for_LPS_Mt = crit_table_LPS[Pooled_df_within_all_LPS_for_LPS_Mt]
        except:
            d_crit_for_LPS_Mt = crit_table_LPS[-1]

        #################################################

        MD_control_1 = abs ( (sum(ten)/len(ten)) - (sum(BBS)/len(BBS)) )
        MD_control_2 = abs ( (sum(thirty)/len(thirty)) - (sum(BBS)/len(BBS)) )
        MD_control_3 = abs ( (sum(hundred)/len(hundred)) - (sum(BBS)/len(BBS)) )
        MD_control_4 = abs ( (sum(threehundred)/len(threehundred)) - (sum(BBS)/len(BBS)) )

        MD_LPS_1 = abs ( (sum(ten_LPS)/len(ten_LPS)) - (sum(BBS_LPS)/len(BBS_LPS)) )
        MD_LPS_2 = abs ( (sum(thirty_LPS)/len(thirty_LPS)) - (sum(BBS_LPS)/len(BBS_LPS)) )
        MD_LPS_3 = abs ( (sum(hundred_LPS)/len(hundred_LPS)) - (sum(BBS_LPS)/len(BBS_LPS)) )
        MD_LPS_4 = abs ( (sum(threehundred_LPS)/len(threehundred_LPS)) - (sum(BBS_LPS)/len(BBS_LPS)) )

        MD_LPS2_1 = abs ( (sum(ten_LPS2)/len(ten_LPS2)) - (sum(BBS_LPS2)/len(BBS_LPS2)) )
        MD_LPS2_2 = abs ( (sum(thirty_LPS2)/len(thirty_LPS2)) - (sum(BBS_LPS2)/len(BBS_LPS2)) )
        MD_LPS2_3 = abs ( (sum(hundred_LPS2)/len(hundred_LPS2)) - (sum(BBS_LPS2)/len(BBS_LPS2)) )
        MD_LPS2_4 = abs ( (sum(threehundred_LPS2)/len(threehundred_LPS2)) - (sum(BBS_LPS2)/len(BBS_LPS2)) )

        MD_LPS3_1 = abs ( (sum(ten_LPS3)/len(ten_LPS3)) - (sum(BBS_LPS3)/len(BBS_LPS3)) )
        MD_LPS3_2 = abs ( (sum(thirty_LPS3)/len(thirty_LPS3)) - (sum(BBS_LPS3)/len(BBS_LPS3)) )
        MD_LPS3_3 = abs ( (sum(hundred_LPS3)/len(hundred_LPS3)) - (sum(BBS_LPS3)/len(BBS_LPS3)) )
        MD_LPS3_4 = abs ( (sum(threehundred_LPS3)/len(threehundred_LPS3)) - (sum(BBS_LPS3)/len(BBS_LPS3)) )

        MD_BBS_by_LPS_1 = abs ( (sum(BBS_LPS)/len(BBS_LPS)) - (sum(BBS)/len(BBS)) )
        MD_BBS_by_LPS_2 = abs ( (sum(BBS_LPS2)/len(BBS_LPS2)) - (sum(BBS)/len(BBS)) )
        MD_BBS_by_LPS_3 = abs ( (sum(BBS_LPS3)/len(BBS_LPS3)) - (sum(BBS)/len(BBS)) )

        MD_ten_by_LPS_1 = abs ( (sum(ten_LPS)/len(ten_LPS)) - (sum(ten)/len(ten)) )
        MD_ten_by_LPS_2 = abs ( (sum(ten_LPS2)/len(ten_LPS2)) - (sum(ten)/len(ten)) )
        MD_ten_by_LPS_3 = abs ( (sum(ten_LPS3)/len(ten_LPS3)) - (sum(ten)/len(ten)) )

        MD_thirty_by_LPS_1 = abs ( (sum(thirty_LPS)/len(thirty_LPS)) - (sum(thirty)/len(thirty)) )
        MD_thirty_by_LPS_2 = abs ( (sum(thirty_LPS2)/len(thirty_LPS2)) - (sum(thirty)/len(thirty)) )
        MD_thirty_by_LPS_3 = abs ( (sum(thirty_LPS3)/len(thirty_LPS3)) - (sum(thirty)/len(thirty)) )

        MD_hundred_by_LPS_1 = abs ( (sum(hundred_LPS)/len(hundred_LPS)) - (sum(hundred)/len(hundred)) )
        MD_hundred_by_LPS_2 = abs ( (sum(hundred_LPS2)/len(hundred_LPS2)) - (sum(hundred)/len(hundred)) )
        MD_hundred_by_LPS_3 = abs ( (sum(hundred_LPS3)/len(hundred_LPS3)) - (sum(hundred)/len(hundred)) )

        MD_threehundred_by_LPS_1 = abs ( (sum(threehundred_LPS)/len(threehundred_LPS)) - (sum(threehundred)/len(threehundred)) )
        MD_threehundred_by_LPS_2 = abs ( (sum(threehundred_LPS2)/len(threehundred_LPS2)) - (sum(threehundred)/len(threehundred)) )
        MD_threehundred_by_LPS_3 = abs ( (sum(threehundred_LPS3)/len(threehundred_LPS3)) - (sum(threehundred)/len(threehundred)) )

        MD_Mt_by_LPS_1 = abs ( (sum(Mt_LPS)/len(Mt_LPS)) - (sum(Mt)/len(Mt)) )
        MD_Mt_by_LPS_2 = abs ( (sum(Mt_LPS2)/len(Mt_LPS2)) - (sum(Mt)/len(Mt)) )
        MD_Mt_by_LPS_3 = abs ( (sum(Mt_LPS3)/len(Mt_LPS3)) - (sum(Mt)/len(Mt)) )

        #################################################

        d_control_1 = d_crit_control * math.sqrt (MSE_within_control * ( (1/len(ten)) + (1/len(BBS)) ))
        d_control_2 = d_crit_control * math.sqrt (MSE_within_control * ( (1/len(thirty)) + (1/len(BBS)) ))
        d_control_3 = d_crit_control * math.sqrt (MSE_within_control * ( (1/len(hundred)) + (1/len(BBS)) ))
        d_control_4 = d_crit_control * math.sqrt (MSE_within_control * ( (1/len(threehundred)) + (1/len(BBS)) ))

        d_LPS_1 = d_crit_LPS * math.sqrt (MSE_within_LPS * ( (1/len(ten_LPS)) + (1/len(BBS_LPS)) ))
        d_LPS_2 = d_crit_LPS * math.sqrt (MSE_within_LPS * ( (1/len(thirty_LPS)) + (1/len(BBS_LPS)) ))
        d_LPS_3 = d_crit_LPS * math.sqrt (MSE_within_LPS * ( (1/len(hundred_LPS)) + (1/len(BBS_LPS)) ))
        d_LPS_4 = d_crit_LPS * math.sqrt (MSE_within_LPS * ( (1/len(threehundred_LPS)) + (1/len(BBS_LPS)) ))

        d_LPS2_1 = d_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ( (1/len(ten_LPS2)) + (1/len(BBS_LPS2)) ))
        d_LPS2_2 = d_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ( (1/len(thirty_LPS2)) + (1/len(BBS_LPS2)) ))
        d_LPS2_3 = d_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ( (1/len(hundred_LPS2)) + (1/len(BBS_LPS2)) ))
        d_LPS2_4 = d_crit_LPS2 * math.sqrt (MSE_within_LPS2 * ( (1/len(threehundred_LPS2)) + (1/len(BBS_LPS2)) ))

        d_LPS3_1 = d_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ( (1/len(ten_LPS3)) + (1/len(BBS_LPS3)) ))
        d_LPS3_2 = d_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ( (1/len(thirty_LPS3)) + (1/len(BBS_LPS3)) ))
        d_LPS3_3 = d_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ( (1/len(hundred_LPS3)) + (1/len(BBS_LPS3)) ))
        d_LPS3_4 = d_crit_LPS3 * math.sqrt (MSE_within_LPS3 * ( (1/len(threehundred_LPS3)) + (1/len(BBS_LPS3)) ))

        d_BBS_by_LPS_1 = d_crit_BBS_by_LPS * math.sqrt (MSE_within_BBS_by_LPS * ( (1/len(BBS_LPS)) + (1/len(BBS)) ))
        d_BBS_by_LPS_2 = d_crit_BBS_by_LPS * math.sqrt (MSE_within_BBS_by_LPS * ( (1/len(BBS_LPS2)) + (1/len(BBS)) ))
        d_BBS_by_LPS_3 = d_crit_BBS_by_LPS * math.sqrt (MSE_within_BBS_by_LPS * ( (1/len(BBS_LPS3)) + (1/len(BBS)) ))

        d_ten_by_LPS_1 = d_crit_ten_by_LPS * math.sqrt (MSE_within_ten_by_LPS * ( (1/len(ten_LPS)) + (1/len(ten)) ))
        d_ten_by_LPS_2 = d_crit_ten_by_LPS * math.sqrt (MSE_within_ten_by_LPS * ( (1/len(ten_LPS2)) + (1/len(ten)) ))
        d_ten_by_LPS_3 = d_crit_ten_by_LPS * math.sqrt (MSE_within_ten_by_LPS * ( (1/len(ten_LPS3)) + (1/len(ten)) ))

        d_thirty_by_LPS_1 = d_crit_thirty_by_LPS * math.sqrt (MSE_within_thirty_by_LPS * ( (1/len(thirty_LPS)) + (1/len(thirty)) ))
        d_thirty_by_LPS_2 = d_crit_thirty_by_LPS * math.sqrt (MSE_within_thirty_by_LPS * ( (1/len(thirty_LPS2)) + (1/len(thirty)) ))
        d_thirty_by_LPS_3 = d_crit_thirty_by_LPS * math.sqrt (MSE_within_thirty_by_LPS * ( (1/len(thirty_LPS3)) + (1/len(thirty)) ))

        d_hundred_by_LPS_1 = d_crit_hundred_by_LPS * math.sqrt (MSE_within_hundred_by_LPS * ( (1/len(hundred_LPS)) + (1/len(hundred)) ))
        d_hundred_by_LPS_2 = d_crit_hundred_by_LPS * math.sqrt (MSE_within_hundred_by_LPS * ( (1/len(hundred_LPS2)) + (1/len(hundred)) ))
        d_hundred_by_LPS_3 = d_crit_hundred_by_LPS * math.sqrt (MSE_within_hundred_by_LPS * ( (1/len(hundred_LPS3)) + (1/len(hundred)) ))

        d_threehundred_by_LPS_1 = d_crit_threehundred_by_LPS * math.sqrt (MSE_within_threehundred_by_LPS * ( (1/len(threehundred_LPS)) + (1/len(threehundred)) ))
        d_threehundred_by_LPS_2 = d_crit_threehundred_by_LPS * math.sqrt (MSE_within_threehundred_by_LPS * ( (1/len(threehundred_LPS2)) + (1/len(threehundred)) ))
        d_threehundred_by_LPS_3 = d_crit_threehundred_by_LPS * math.sqrt (MSE_within_threehundred_by_LPS * ( (1/len(threehundred_LPS3)) + (1/len(threehundred)) ))

        d_Mt_by_LPS_1 = d_crit_Mt_by_LPS * math.sqrt (MSE_within_Mt_by_LPS * ( (1/len(Mt_LPS)) + (1/len(Mt)) ))
        d_Mt_by_LPS_2 = d_crit_Mt_by_LPS * math.sqrt (MSE_within_Mt_by_LPS * ( (1/len(Mt_LPS2)) + (1/len(Mt)) ))
        d_Mt_by_LPS_3 = d_crit_Mt_by_LPS * math.sqrt (MSE_within_Mt_by_LPS * ( (1/len(Mt_LPS3)) + (1/len(Mt)) ))

        #################################################

        SED_control_1 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        SED_LPS2_1 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
        SED_LPS2_2 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
        SED_LPS2_3 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
        SED_LPS2_4 = math.sqrt(MSE_within_LPS2 * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

        SED_LPS3_1 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
        SED_LPS3_2 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
        SED_LPS3_3 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
        SED_LPS3_4 = math.sqrt(MSE_within_LPS3 * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

        SED_BBS_by_LPS_1 = math.sqrt(MSE_within_BBS_by_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
        SED_BBS_by_LPS_2 = math.sqrt(MSE_within_BBS_by_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
        SED_BBS_by_LPS_3 = math.sqrt(MSE_within_BBS_by_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))

        SED_ten_by_LPS_1 = math.sqrt(MSE_within_ten_by_LPS * ( (1/len(ten)) + (1/len(ten_LPS)) ))
        SED_ten_by_LPS_2 = math.sqrt(MSE_within_ten_by_LPS * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
        SED_ten_by_LPS_3 = math.sqrt(MSE_within_ten_by_LPS * ( (1/len(ten)) + (1/len(ten_LPS3)) ))

        SED_thirty_by_LPS_1 = math.sqrt(MSE_within_thirty_by_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
        SED_thirty_by_LPS_2 = math.sqrt(MSE_within_thirty_by_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
        SED_thirty_by_LPS_3 = math.sqrt(MSE_within_thirty_by_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))

        SED_hundred_by_LPS_1 = math.sqrt(MSE_within_hundred_by_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
        SED_hundred_by_LPS_2 = math.sqrt(MSE_within_hundred_by_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
        SED_hundred_by_LPS_3 = math.sqrt(MSE_within_hundred_by_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))

        SED_threehundred_by_LPS_1 = math.sqrt(MSE_within_threehundred_by_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))
        SED_threehundred_by_LPS_2 = math.sqrt(MSE_within_threehundred_by_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))
        SED_threehundred_by_LPS_3 = math.sqrt(MSE_within_threehundred_by_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

        SED_Mt_by_LPS_1 = math.sqrt(MSE_within_Mt_by_LPS * ( (1/len(Mt)) + (1/len(Mt_LPS)) ))
        SED_Mt_by_LPS_2 = math.sqrt(MSE_within_Mt_by_LPS * ( (1/len(Mt)) + (1/len(Mt_LPS2)) ))
        SED_Mt_by_LPS_3 = math.sqrt(MSE_within_Mt_by_LPS * ( (1/len(Mt)) + (1/len(Mt_LPS3)) ))

        #################################################

        q_control_1 = MD_control_1/SED_control_1
        q_control_2 = MD_control_2/SED_control_2
        q_control_3 = MD_control_3/SED_control_3
        q_control_4 = MD_control_4/SED_control_4

        q_LPS_1 = MD_LPS_1/SED_LPS_1
        q_LPS_2 = MD_LPS_2/SED_LPS_2
        q_LPS_3 = MD_LPS_3/SED_LPS_3
        q_LPS_4 = MD_LPS_4/SED_LPS_4

        q_LPS2_1 = MD_LPS2_1/SED_LPS2_1
        q_LPS2_2 = MD_LPS2_2/SED_LPS2_2
        q_LPS2_3 = MD_LPS2_3/SED_LPS2_3
        q_LPS2_4 = MD_LPS2_4/SED_LPS2_4

        q_LPS3_1 = MD_LPS3_1/SED_LPS3_1
        q_LPS3_2 = MD_LPS3_2/SED_LPS3_2
        q_LPS3_3 = MD_LPS3_3/SED_LPS3_3
        q_LPS3_4 = MD_LPS3_4/SED_LPS3_4

        q_BBS_by_LPS_1 = MD_BBS_by_LPS_1/SED_BBS_by_LPS_1
        q_BBS_by_LPS_2 = MD_BBS_by_LPS_2/SED_BBS_by_LPS_2
        q_BBS_by_LPS_3 = MD_BBS_by_LPS_3/SED_BBS_by_LPS_3

        q_ten_by_LPS_1 = MD_ten_by_LPS_1/SED_ten_by_LPS_1
        q_ten_by_LPS_2 = MD_ten_by_LPS_2/SED_ten_by_LPS_2
        q_ten_by_LPS_3 = MD_ten_by_LPS_3/SED_ten_by_LPS_3

        q_thirty_by_LPS_1 = MD_thirty_by_LPS_1/SED_thirty_by_LPS_1
        q_thirty_by_LPS_2 = MD_thirty_by_LPS_2/SED_thirty_by_LPS_2
        q_thirty_by_LPS_3 = MD_thirty_by_LPS_3/SED_thirty_by_LPS_3

        q_hundred_by_LPS_1 = MD_hundred_by_LPS_1/SED_hundred_by_LPS_1
        q_hundred_by_LPS_2 = MD_hundred_by_LPS_2/SED_hundred_by_LPS_2
        q_hundred_by_LPS_3 = MD_hundred_by_LPS_3/SED_hundred_by_LPS_3

        q_threehundred_by_LPS_1 = MD_threehundred_by_LPS_1/SED_threehundred_by_LPS_1
        q_threehundred_by_LPS_2 = MD_threehundred_by_LPS_2/SED_threehundred_by_LPS_2
        q_threehundred_by_LPS_3 = MD_threehundred_by_LPS_3/SED_threehundred_by_LPS_3

        q_Mt_by_LPS_1 = MD_Mt_by_LPS_1/SED_Mt_by_LPS_1
        q_Mt_by_LPS_2 = MD_Mt_by_LPS_2/SED_Mt_by_LPS_2
        q_Mt_by_LPS_3 = MD_Mt_by_LPS_3/SED_Mt_by_LPS_3

        q_list = [q_control_1, q_control_2, q_control_3, q_control_4, q_LPS_1, q_LPS_2, q_LPS_3, q_LPS_4, q_LPS2_1, q_LPS2_2, q_LPS2_3, q_LPS2_4, q_LPS3_1, q_LPS3_2, q_LPS3_3, q_LPS3_4,
                  q_BBS_by_LPS_1,q_BBS_by_LPS_2,q_BBS_by_LPS_3,q_ten_by_LPS_1,q_ten_by_LPS_2,q_ten_by_LPS_3,q_thirty_by_LPS_1,q_thirty_by_LPS_2,q_thirty_by_LPS_3,
                  q_hundred_by_LPS_1,q_hundred_by_LPS_2,q_hundred_by_LPS_3,q_threehundred_by_LPS_1,q_threehundred_by_LPS_2,q_threehundred_by_LPS_3,
                  q_Mt_by_LPS_1, q_Mt_by_LPS_2, q_Mt_by_LPS_3]

        #################################################

        #FOR TWO-WAY ANOVA WHERE 1 FAMILY ERROR IS USED

        test_parameter = input("Is this a test using 1 pooled MSE? <y/n> ")

        if test_parameter == "y" or test_parameter == "Y":
        
            print("NEW RESULTS STARTING BELOW")

            d_control_1 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(ten)) + (1/len(BBS)) ))
            d_control_2 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(thirty)) + (1/len(BBS)) ))
            d_control_3 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(hundred)) + (1/len(BBS)) ))
            d_control_4 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(threehundred)) + (1/len(BBS)) ))

            d_LPS_1 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(ten_LPS)) + (1/len(BBS_LPS)) ))
            d_LPS_2 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(thirty_LPS)) + (1/len(BBS_LPS)) ))
            d_LPS_3 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(hundred_LPS)) + (1/len(BBS_LPS)) ))
            d_LPS_4 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(threehundred_LPS)) + (1/len(BBS_LPS)) ))

            d_LPS2_1 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(ten_LPS2)) + (1/len(BBS_LPS2)) ))
            d_LPS2_2 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(thirty_LPS2)) + (1/len(BBS_LPS2)) ))
            d_LPS2_3 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(hundred_LPS2)) + (1/len(BBS_LPS2)) ))
            d_LPS2_4 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(threehundred_LPS2)) + (1/len(BBS_LPS2)) ))

            d_LPS3_1 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(ten_LPS3)) + (1/len(BBS_LPS3)) ))
            d_LPS3_2 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(thirty_LPS3)) + (1/len(BBS_LPS3)) ))
            d_LPS3_3 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(hundred_LPS3)) + (1/len(BBS_LPS3)) ))
            d_LPS3_4 = d_crit_for_vaccae * math.sqrt (Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(threehundred_LPS3)) + (1/len(BBS_LPS3)) ))

            d_BBS_by_LPS_1 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS_LPS)) + (1/len(BBS)) ))
            d_BBS_by_LPS_2 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS_LPS2)) + (1/len(BBS)) ))
            d_BBS_by_LPS_3 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS_LPS3)) + (1/len(BBS)) ))

            d_ten_by_LPS_1 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten_LPS)) + (1/len(ten)) ))
            d_ten_by_LPS_2 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten_LPS2)) + (1/len(ten)) ))
            d_ten_by_LPS_3 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten_LPS3)) + (1/len(ten)) ))

            d_thirty_by_LPS_1 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty_LPS)) + (1/len(thirty)) ))
            d_thirty_by_LPS_2 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty_LPS2)) + (1/len(thirty)) ))
            d_thirty_by_LPS_3 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty_LPS3)) + (1/len(thirty)) ))

            d_hundred_by_LPS_1 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred_LPS)) + (1/len(hundred)) ))
            d_hundred_by_LPS_2 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred_LPS2)) + (1/len(hundred)) ))
            d_hundred_by_LPS_3 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred_LPS3)) + (1/len(hundred)) ))

            d_threehundred_by_LPS_1 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred_LPS)) + (1/len(threehundred)) ))
            d_threehundred_by_LPS_2 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred_LPS2)) + (1/len(threehundred)) ))
            d_threehundred_by_LPS_3 = d_crit_for_LPS * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred_LPS3)) + (1/len(threehundred)) ))

            d_Mt_by_LPS_1 = d_crit_for_LPS_Mt * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt_LPS)) + (1/len(Mt)) ))
            d_Mt_by_LPS_2 = d_crit_for_LPS_Mt * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt_LPS2)) + (1/len(Mt)) ))
            d_Mt_by_LPS_3 = d_crit_for_LPS_Mt * math.sqrt (Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt_LPS3)) + (1/len(Mt)) ))



            SED_control_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS)) + (1/len(ten)) ))
            SED_control_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS)) + (1/len(thirty)) ))
            SED_control_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS)) + (1/len(hundred)) ))
            SED_control_4 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS)) + (1/len(threehundred)) ))

            SED_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
            SED_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
            SED_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
            SED_LPS_4 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

            SED_LPS2_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS2)) + (1/len(ten_LPS2)) ))
            SED_LPS2_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS2)) + (1/len(thirty_LPS2)) ))
            SED_LPS2_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS2)) + (1/len(hundred_LPS2)) ))
            SED_LPS2_4 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS2)) + (1/len(threehundred_LPS2)) ))

            SED_LPS3_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS3)) + (1/len(ten_LPS3)) ))
            SED_LPS3_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS3)) + (1/len(thirty_LPS3)) ))
            SED_LPS3_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS3)) + (1/len(hundred_LPS3)) ))
            SED_LPS3_4 = math.sqrt(Pooled_MSE_within_all_LPS_for_vaccae * ( (1/len(BBS_LPS3)) + (1/len(threehundred_LPS3)) ))

            SED_BBS_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS)) ))
            SED_BBS_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS2)) ))
            SED_BBS_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(BBS)) + (1/len(BBS_LPS3)) ))

            SED_ten_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten)) + (1/len(ten_LPS)) ))
            SED_ten_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten)) + (1/len(ten_LPS2)) ))
            SED_ten_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(ten)) + (1/len(ten_LPS3)) ))

            SED_thirty_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS)) ))
            SED_thirty_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS2)) ))
            SED_thirty_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(thirty)) + (1/len(thirty_LPS3)) ))

            SED_hundred_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS)) ))
            SED_hundred_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS2)) ))
            SED_hundred_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(hundred)) + (1/len(hundred_LPS3)) ))

            SED_threehundred_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS)) ))
            SED_threehundred_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS2)) ))
            SED_threehundred_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS * ( (1/len(threehundred)) + (1/len(threehundred_LPS3)) ))

            SED_Mt_by_LPS_1 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt)) + (1/len(Mt_LPS)) ))
            SED_Mt_by_LPS_2 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt)) + (1/len(Mt_LPS2)) ))
            SED_Mt_by_LPS_3 = math.sqrt(Pooled_MSE_within_all_LPS_for_LPS_Mt * ( (1/len(Mt)) + (1/len(Mt_LPS3)) ))



            q_control_1 = MD_control_1/SED_control_1
            q_control_2 = MD_control_2/SED_control_2
            q_control_3 = MD_control_3/SED_control_3
            q_control_4 = MD_control_4/SED_control_4

            q_LPS_1 = MD_LPS_1/SED_LPS_1
            q_LPS_2 = MD_LPS_2/SED_LPS_2
            q_LPS_3 = MD_LPS_3/SED_LPS_3
            q_LPS_4 = MD_LPS_4/SED_LPS_4

            q_LPS2_1 = MD_LPS2_1/SED_LPS2_1
            q_LPS2_2 = MD_LPS2_2/SED_LPS2_2
            q_LPS2_3 = MD_LPS2_3/SED_LPS2_3
            q_LPS2_4 = MD_LPS2_4/SED_LPS2_4

            q_LPS3_1 = MD_LPS3_1/SED_LPS3_1
            q_LPS3_2 = MD_LPS3_2/SED_LPS3_2
            q_LPS3_3 = MD_LPS3_3/SED_LPS3_3
            q_LPS3_4 = MD_LPS3_4/SED_LPS3_4

            q_BBS_by_LPS_1 = MD_BBS_by_LPS_1/SED_BBS_by_LPS_1
            q_BBS_by_LPS_2 = MD_BBS_by_LPS_2/SED_BBS_by_LPS_2
            q_BBS_by_LPS_3 = MD_BBS_by_LPS_3/SED_BBS_by_LPS_3

            q_ten_by_LPS_1 = MD_ten_by_LPS_1/SED_ten_by_LPS_1
            q_ten_by_LPS_2 = MD_ten_by_LPS_2/SED_ten_by_LPS_2
            q_ten_by_LPS_3 = MD_ten_by_LPS_3/SED_ten_by_LPS_3

            q_thirty_by_LPS_1 = MD_thirty_by_LPS_1/SED_thirty_by_LPS_1
            q_thirty_by_LPS_2 = MD_thirty_by_LPS_2/SED_thirty_by_LPS_2
            q_thirty_by_LPS_3 = MD_thirty_by_LPS_3/SED_thirty_by_LPS_3

            q_hundred_by_LPS_1 = MD_hundred_by_LPS_1/SED_hundred_by_LPS_1
            q_hundred_by_LPS_2 = MD_hundred_by_LPS_2/SED_hundred_by_LPS_2
            q_hundred_by_LPS_3 = MD_hundred_by_LPS_3/SED_hundred_by_LPS_3

            q_threehundred_by_LPS_1 = MD_threehundred_by_LPS_1/SED_threehundred_by_LPS_1
            q_threehundred_by_LPS_2 = MD_threehundred_by_LPS_2/SED_threehundred_by_LPS_2
            q_threehundred_by_LPS_3 = MD_threehundred_by_LPS_3/SED_threehundred_by_LPS_3

            q_Mt_by_LPS_1 = MD_Mt_by_LPS_1/SED_Mt_by_LPS_1
            q_Mt_by_LPS_2 = MD_Mt_by_LPS_2/SED_Mt_by_LPS_2
            q_Mt_by_LPS_3 = MD_Mt_by_LPS_3/SED_Mt_by_LPS_3

            q_list = [q_control_1, q_control_2, q_control_3, q_control_4, q_LPS_1, q_LPS_2, q_LPS_3, q_LPS_4, q_LPS2_1, q_LPS2_2, q_LPS2_3, q_LPS2_4, q_LPS3_1, q_LPS3_2, q_LPS3_3, q_LPS3_4,
                      q_BBS_by_LPS_1,q_BBS_by_LPS_2,q_BBS_by_LPS_3,q_ten_by_LPS_1,q_ten_by_LPS_2,q_ten_by_LPS_3,q_thirty_by_LPS_1,q_thirty_by_LPS_2,q_thirty_by_LPS_3,
                      q_hundred_by_LPS_1,q_hundred_by_LPS_2,q_hundred_by_LPS_3,q_threehundred_by_LPS_1,q_threehundred_by_LPS_2,q_threehundred_by_LPS_3,
                      q_Mt_by_LPS_1, q_Mt_by_LPS_2, q_Mt_by_LPS_3]

            print("NEW q values:")
            for a in q_list:
                print(a)

            print("NEW RESULTS END HERE (just q values above)")

        #################################################

        print ("========================================================")
        print ("")
        print ("Your Dunnett's Test Results Are:")
        print ("")

        if MD_control_1 > d_control_1:
            print ("BBS vs. 10:  *  | p=") #, p_control_1)
        else:
            print ("BBS vs. 10:  ns | p=")#, p_control_1)

        if MD_control_2 > d_control_2:
            print ("BBS vs. 30:  *  | p=")#, p_control_2)
        else:
            print ("BBS vs. 30:  ns | p=")#, p_control_2)

        if MD_control_3 > d_control_3:
            print ("BBS vs. 100: *  | p=")#, p_control_3)
        else:
            print ("BBS vs. 100: ns | p=")#, p_control_3)

        if MD_control_4 > d_control_4:
            print ("BBS vs. 300: *  | p=")#, p_control_4)
        else:
            print ("BBS vs. 300: ns | p=")#, p_control_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS_1 > d_LPS_1:
            print ("BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS_1)
        else:
            print ("BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS_1)
            
        if MD_LPS_2 > d_LPS_2:
            print ("BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS_2)
        else:
            print ("BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS_2)

        if MD_LPS_3 > d_LPS_3:
            print ("BBS LPS vs. 100 LPS: *  | p=")#, p_LPS_3)
        else:
            print ("BBS LPS vs. 100 LPS: ns | p=")#, p_LPS_3)

        if MD_LPS_4 > d_LPS_4:
            print ("BBS LPS vs. 300 LPS: *  | p=")#, p_LPS_4)
        else:
            print ("BBS LPS vs. 300 LPS: ns | p=")#, p_LPS_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS2_1 > d_LPS2_1:
            print ("2nd Dose BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS2_1)
        else:
            print ("2nd Dose BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS2_1)
            
        if MD_LPS2_2 > d_LPS2_2:
            print ("2nd Dose BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS2_2)
        else:
            print ("2nd Dose BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS2_2)

        if MD_LPS2_3 > d_LPS2_3:
            print ("2nd Dose BBS LPS vs. 100 LPS: *  | p=")#, p_LPS2_3)
        else:
            print ("2nd Dose BBS LPS vs. 100 LPS: ns | p=")#, p_LPS2_3)

        if MD_LPS2_4 > d_LPS2_4:
            print ("2nd Dose BBS LPS vs. 300 LPS: *  | p=")#, p_LPS2_4)
        else:
            print ("2nd Dose BBS LPS vs. 300 LPS: ns | p=")#, p_LPS2_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS3_1 > d_LPS3_1:
            print ("3rd Dose BBS LPS vs. 10 LPS:  *  | p=")#, p_LPS3_1)
        else:
            print ("3rd Dose BBS LPS vs. 10 LPS:  ns | p=")#, p_LPS3_1)
            
        if MD_LPS3_2 > d_LPS3_2:
            print ("3rd Dose BBS LPS vs. 30 LPS:  *  | p=")#, p_LPS3_2)
        else:
            print ("3rd Dose BBS LPS vs. 30 LPS:  ns | p=")#, p_LPS3_2)

        if MD_LPS3_3 > d_LPS3_3:
            print ("3rd Dose BBS LPS vs. 100 LPS: *  | p=")#, p_LPS3_3)
        else:
            print ("3rd Dose BBS LPS vs. 100 LPS: ns | p=")#, p_LPS3_3)

        if MD_LPS3_4 > d_LPS3_4:
            print ("3rd Dose BBS LPS vs. 300 LPS: *  | p=")#, p_LPS3_4)
        else:
            print ("3rd Dose BBS LPS vs. 300 LPS: ns | p=")#, p_LPS3_4)

        print ("")

        print ("q values: ")
        for a in q_list:
            print (a)
        
        print ("")
        print ("========================================================")
        print ("")
        
##        try:
##            shell = sys.stdout.shell
##        except AttributeError:
##            raise RuntimeError("you must run this program in IDLE")
##
##        shell.write("Thank ", "COMMENT")
##        shell.write("You ","KEYWORD")
##        shell.write("For ","STRING")
##        shell.write("Using ","DEFINITION")
##        shell.write("The ","BUILTIN")
##        shell.write("qPCR ","COMMENT")
##        shell.write("Analyzer ","KEYWORD")
##        shell.write("3000!","STRING")
##        print ("")
##        print ("")
##        shell.write("-EH","BUILTIN")
##        print ("")
##        print ("")
##
##    else:
##        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
    print ("========================================================")

    if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
        print ("")
        print ("1: LPS Dose Response")
        print ("")
        print ("2: Control vs. 1 LPS Dose")
        print ("")

        experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

        if experiment != "1" and experiment != "2":
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        print ("========================================================")

        print ("")
        print ("1: Analyze Cq Values")
        print ("")
        print ("2: Analyze Log2 Values")
        print ("")

        program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

        if program == "1" and experiment == "2":
            print ("========================================================")
            print ("")
            print ("1: Analyze a Data Set Using 1 Minimum Value")
            print ("")
            print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
            print ("")
            
            meta_program = input ("Which Method Would You Like to Use? ")
            
            if meta_program == "1":
                print ("========================================================")
                PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
            elif meta_program == "2":
                print ("========================================================")
                #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
            else:
                raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
        elif program == "1" and experiment == "1":
            print ("========================================================")
            print ("")
            print ("1: Analyze a Data Set Using 1 Minimum Value")
            print ("")
            print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
            print ("")
            
            meta_program = input ("Which Method Would You Like to Use? ")
            
            if meta_program == "1":
                print ("========================================================")
                PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
            elif meta_program == "2":
                print ("========================================================")
                #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
            else:
                raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
        elif program == "2" and experiment == "1":
            print ("")
            print ("========================================================")
            meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
            if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
            elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
            else:
                raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        elif program == "2" and experiment == "2":
            print ("")
            print ("========================================================")
            meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
            if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
            elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
            else:
                raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


    elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 

        #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
        print ("")
        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")
        print ("========================================================")

    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    return;

def PCR_posthoc_test (N, BBS, ten, thirty, hundred, threehundred, Mt, BBS_LPS, ten_LPS, thirty_LPS, hundred_LPS,
                 threehundred_LPS, Mt_LPS):

    import statistics
    import math
    import scipy.stats
    
    print ("")
    print ("1: Unprotected Fisher's LSD: doesn't correct for multiple comparisons")
    print ("")
    print ("2: Tukey's HSD: corrects for multiple comparisons and assumes every treatment has an equal n (no NaN values)")
    print ("")
    print ("3: Tukey-Cramer Method: corrects for multiple comparisons, but doesn't assume equal n for each treatment")
    print ("")
    answer = input ("Which Post-Hoc Test Would You Like to Run? (Enter The Number): ")

    if answer == "1":
    
        #Fisher LSD: any difference larger than the LSD is significant

        #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/New T Critical Values.csv", header = None)
        crit_table = [12.7065, 4.3026, 3.1824, 2.7764, 2.5706, 2.4469, 2.3646, 2.306, 2.2621, 2.2282, 2.201, 2.1788, 2.1604,	2.1448,	2.1314,	2.1199,	2.1098,	2.1009,	2.093, 2.086, 2.0796, 2.0739, 2.0686, 2.0639, 2.0596, 2.0555, 2.0518, 2.0484, 2.0452, 2.0423, 2.0395, 2.0369, 2.0345, 2.0322, 2.0301, 2.0281, 2.0262, 2.0244, 2.0227, 2.0211, 2.0196, 2.0181, 2.0167, 2.0154, 2.0141, 2.0129, 2.0117, 2.0106, 2.0096, 2.0086, 2.0076, 2.0066, 2.0057, 2.0049, 2.0041, 2.0032, 2.0025, 2.0017, 2.001, 2.0003, 1.9996, 1.999, 1.9983, 1.9977, 1.9971, 1.9966, 1.996, 1.9955, 1.995, 1.9944, 1.9939, 1.9935, 1.993,1.9925, 1.9921, 1.9917, 1.9913, 1.9909, 1.9904, 1.9901, 1.9897, 1.9893, 1.9889, 1.9886, 1.9883, 1.9879, 1.9876, 1.9873, 1.987, 1.9867, 1.9864, 1.9861, 1.9858, 1.9855, 1.9852, 1.985, 1.9847, 1.9845, 1.9842, 1.984, 1.9837, 1.9835, 1.9833, 1.983, 1.9828, 1.9826, 1.9824, 1.9822, 1.982, 1.9818, 1.9816, 1.9814, 1.9812, 1.981, 1.9808, 1.9806, 1.9805, 1.9803, 1.9801, 1.9799]
        
        #################################################
      
        df_BBS = len (BBS) - 1
        df_ten = len (ten) - 1
        df_thirty = len (thirty) - 1
        df_hundred = len (hundred) - 1
        df_threehundred = len (threehundred) - 1
        df_Mt = len (Mt) - 1

        df_BBS_LPS = len (BBS_LPS) - 1
        df_ten_LPS = len (ten_LPS) - 1
        df_thirty_LPS = len (thirty_LPS) - 1
        df_hundred_LPS = len (hundred_LPS) - 1
        df_threehundred_LPS = len (threehundred_LPS) - 1
        df_Mt_LPS = len (Mt_LPS) - 1

        #################################################
        
        if len (BBS) == 1:
            BBS_var = 0
        else:
            BBS_var = statistics.variance (BBS)

        if len(ten) == 1:
            ten_var = 0
        else:
            ten_var = statistics.variance (ten)

        if len(thirty) == 1:
            thirty_var = 0
        else:
            thirty_var = statistics.variance (thirty)

        if len(hundred) == 1:
            hundred_var = 0
        else:
            hundred_var = statistics.variance (hundred)

        if len (threehundred) == 1:
            threehundred_var = 0
        else:
            threehundred_var = statistics.variance (threehundred)

        if len (Mt) == 1:
            Mt_var = 0
        else:
            Mt_var = statistics.variance (Mt)

        ############################

        if len (BBS_LPS) == 1:
            BBS_LPS_var = 0
        else:
            BBS_LPS_var = statistics.variance (BBS_LPS)

        if len(ten_LPS) == 1:
            ten_LPS_var = 0
        else:
            ten_LPS_var = statistics.variance (ten_LPS)

        if len(thirty_LPS) == 1:
            thirty_LPS_var = 0
        else:
            thirty_LPS_var = statistics.variance (thirty_LPS)

        if len(hundred_LPS) == 1:
            hundred_LPS_var = 0
        else:
            hundred_LPS_var = statistics.variance (hundred_LPS)

        if len (threehundred_LPS) == 1:
            threehundred_LPS_var = 0
        else:
            threehundred_LPS_var = statistics.variance (threehundred_LPS)

        if len (Mt_LPS) == 1:
            Mt_LPS_var = 0
        else:
            Mt_LPS_var = statistics.variance (Mt_LPS)

        #################################################

        SS_within_control = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var)
        df_within_control = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred
        df_within_control_actual = int (float(df_within_control) - 1)
        MSE_within_control = SS_within_control / df_within_control

        SS_within_LPS = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var)
        df_within_LPS = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS
        df_within_LPS_actual = int (float (df_within_LPS) - 1)
        MSE_within_LPS = SS_within_LPS / df_within_LPS
        
        #################################################

        #t_crit_control = crit_table.iloc[df_within_control_actual, 0]
        #t_crit_LPS = crit_table.iloc[df_within_LPS_actual, 0]

        t_crit_control = crit_table[df_within_control_actual]
        t_crit_LPS = crit_table[df_within_LPS_actual]

        #################################################
     
        LSD_BBS_ten = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        LSD_BBS_thirty = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        LSD_BBS_hundred = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        LSD_BBS_threehundred = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))
        LSD_BBS_Mt = abs (((sum (BBS) / len (BBS)) - (sum (Mt) / len (Mt))))

        LSD_BBS_ten_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        LSD_BBS_thirty_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        LSD_BBS_hundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        LSD_BBS_threehundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))
        LSD_BBS_Mt_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (Mt_LPS) / len (Mt_LPS))))

        #################################################
        
        LSD_control_1 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(ten))))
        LSD_control_2 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(thirty))))
        LSD_control_3 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(hundred))))
        LSD_control_4 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(threehundred))))
        LSD_control_5 = t_crit_control * math.sqrt (MSE_within_control * ((1/len(BBS)) + (1/len(Mt))))

        LSD_LPS_1 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(ten_LPS))))
        LSD_LPS_2 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(thirty_LPS))))
        LSD_LPS_3 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(hundred_LPS))))
        LSD_LPS_4 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(threehundred_LPS))))
        LSD_LPS_5 = t_crit_LPS * math.sqrt (MSE_within_LPS * ((1/len(BBS_LPS)) + (1/len(Mt_LPS))))

        #################################################
        
        SED_control_1 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt(MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt(MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        #################################################

        tscore_control_1 = LSD_BBS_ten / SED_control_1
        tscore_control_2 = LSD_BBS_thirty / SED_control_2
        tscore_control_3 = LSD_BBS_hundred / SED_control_3
        tscore_control_4 = LSD_BBS_threehundred / SED_control_4

        tscore_LPS_1 = LSD_BBS_ten_LPS / SED_LPS_1
        tscore_LPS_2 = LSD_BBS_thirty_LPS / SED_LPS_2
        tscore_LPS_3 = LSD_BBS_hundred_LPS / SED_LPS_3
        tscore_LPS_4 = LSD_BBS_threehundred_LPS / SED_LPS_4

        #################################################

        pvalue_control_1 = 2 * scipy.stats.t.sf (tscore_control_1, df_within_control)
        pvalue_control_2 = 2 * scipy.stats.t.sf (tscore_control_2, df_within_control)
        pvalue_control_3 = 2 * scipy.stats.t.sf (tscore_control_3, df_within_control)
        pvalue_control_4 = 2 * scipy.stats.t.sf (tscore_control_4, df_within_control)

        pvalue_LPS_1 = 2 * scipy.stats.t.sf (tscore_LPS_1, df_within_LPS)
        pvalue_LPS_2 = 2 * scipy.stats.t.sf (tscore_LPS_2, df_within_LPS)
        pvalue_LPS_3 = 2 * scipy.stats.t.sf (tscore_LPS_3, df_within_LPS)
        pvalue_LPS_4 = 2 * scipy.stats.t.sf (tscore_LPS_4, df_within_LPS)

        print ("========================================================")
        print ("")
        print ("Your Fisher's LSD Results Are:")
        print ("")
        
        if LSD_BBS_ten > LSD_control_1:
            print ("BBS vs. 10:  *  | p =", pvalue_control_1)
        else:
            print ("BBS vs. 10:  ns | p =", pvalue_control_1)

        if LSD_BBS_thirty > LSD_control_2:
            print ("BBS vs. 30:  *  | p =", pvalue_control_2)
        else:
            print ("BBS vs. 30:  ns | p =", pvalue_control_2)

        if LSD_BBS_hundred > LSD_control_3:
            print ("BBS vs. 100: *  | p =", pvalue_control_3)
        else:
            print ("BBS vs. 100: ns | p =", pvalue_control_3)

        if LSD_BBS_threehundred > LSD_control_4:
            print ("BBS vs. 300: *  | p =", pvalue_control_4)
        else:
            print ("BBS vs. 300: ns | p =", pvalue_control_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        
        if LSD_BBS_ten_LPS > LSD_LPS_1:
            print ("BBS LPS vs. 10 LPS:  *  | p =", pvalue_LPS_1)
        else:
            print ("BBS LPS vs. 10 LPS:  ns | p =", pvalue_LPS_1)

        if LSD_BBS_thirty_LPS > LSD_LPS_2:
            print ("BBS LPS vs. 30 LPS:  *  | p =", pvalue_LPS_2)
        else:
            print ("BBS LPS vs. 30 LPS:  ns | p =", pvalue_LPS_2)

        if LSD_BBS_hundred_LPS > LSD_LPS_3:
            print ("BBS LPS vs. 100 LPS: *  | p =", pvalue_LPS_3)
        else:
            print ("BBS LPS vs. 100 LPS: ns | p =", pvalue_LPS_3)

        if LSD_BBS_threehundred_LPS > LSD_LPS_4:
            print ("BBS LPS vs. 300 LPS: *  | p =", pvalue_LPS_4)
        else:
            print ("BBS LPS vs. 300 LPS: ns | p =", pvalue_LPS_4)

        print ("")
        print ("========================================================")
        print ("")
        
        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")

    elif answer == "2":

        #Tukey's HSD: if the HSD is higher than Q critical value, il est significant

        from statsmodels.stats.libqsturng import psturng
        
        #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/Q critical values.csv", header = None)
        crit_table = [37.082,	10.881,	7.502,	6.287,	5.673,	5.305,	5.06,	4.886,	4.755,	4.654,	4.574,	4.508,	4.453,	4.407,	4.367,	4.333,	4.303,	4.276,	4.253,	4.232,	4.213,	4.196,	4.18,	4.166,	4.153,	4.141,	4.13,	4.12,	4.111,	4.102,	4.094,	4.086,	4.079,	4.072,	4.066,	4.06,	4.054,	4.049,	4.044,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887]

        #################################################
      
        df_BBS = len (BBS) - 1
        df_ten = len (ten) - 1
        df_thirty = len (thirty) - 1
        df_hundred = len (hundred) - 1
        df_threehundred = len (threehundred) - 1
        df_Mt = len (Mt) - 1

        df_BBS_LPS = len (BBS_LPS) - 1
        df_ten_LPS = len (ten_LPS) - 1
        df_thirty_LPS = len (thirty_LPS) - 1
        df_hundred_LPS = len (hundred_LPS) - 1
        df_threehundred_LPS = len (threehundred_LPS) - 1
        df_Mt_LPS = len (Mt_LPS) - 1

        #################################################
        
        if len (BBS) == 1:
            BBS_var = 0
        else:
            BBS_var = statistics.variance (BBS)

        if len(ten) == 1:
            ten_var = 0
        else:
            ten_var = statistics.variance (ten)

        if len(thirty) == 1:
            thirty_var = 0
        else:
            thirty_var = statistics.variance (thirty)

        if len(hundred) == 1:
            hundred_var = 0
        else:
            hundred_var = statistics.variance (hundred)

        if len (threehundred) == 1:
            threehundred_var = 0
        else:
            threehundred_var = statistics.variance (threehundred)

        if len (Mt) == 1:
            Mt_var = 0
        else:
            Mt_var = statistics.variance (Mt)

        ############################

        if len (BBS_LPS) == 1:
            BBS_LPS_var = 0
        else:
            BBS_LPS_var = statistics.variance (BBS_LPS)

        if len(ten_LPS) == 1:
            ten_LPS_var = 0
        else:
            ten_LPS_var = statistics.variance (ten_LPS)

        if len(thirty_LPS) == 1:
            thirty_LPS_var = 0
        else:
            thirty_LPS_var = statistics.variance (thirty_LPS)

        if len(hundred_LPS) == 1:
            hundred_LPS_var = 0
        else:
            hundred_LPS_var = statistics.variance (hundred_LPS)

        if len (threehundred_LPS) == 1:
            threehundred_LPS_var = 0
        else:
            threehundred_LPS_var = statistics.variance (threehundred_LPS)

        if len (Mt_LPS) == 1:
            Mt_LPS_var = 0
        else:
            Mt_LPS_var = statistics.variance (Mt_LPS)

        #################################################

        SS_within_control = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var)
        df_within_control = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred
        MSE_within_control = SS_within_control / df_within_control
        
        SS_within_LPS = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var)
        df_within_LPS = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS
        MSE_within_LPS = SS_within_LPS / df_within_LPS

        #################################################

        HSD_BBS_ten = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_thirty = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_hundred = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred)))) / (math.sqrt (MSE_within_control / float(N)))
        HSD_BBS_threehundred = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred)))) / (math.sqrt (MSE_within_control / float(N)))

        HSD_BBS_ten_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_thirty_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_hundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))
        HSD_BBS_threehundred_LPS = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS)))) / (math.sqrt (MSE_within_LPS / float(N)))

        #################################################

        #HSD_control = crit_table.iloc[int (float(df_within_control)-1), 0]
        #HSD_LPS = crit_table.iloc[int (float(df_within_LPS)-1), 0]

        HSD_control = crit_table[int (float(df_within_control)-1)]
        HSD_LPS = crit_table[int (float(df_within_LPS)-1)]

        #################################################

        MD_control_1 = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        MD_control_2 = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        MD_control_3 = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        MD_control_4 = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))

        MD_LPS_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_LPS_2 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_LPS_3 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_LPS_4 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        #################################################

        SED_control_1 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        #################################################

        q_control_1 = MD_control_1 / (SED_control_1 * 1/math.sqrt(2))
        q_control_2 = MD_control_2 / (SED_control_2 * 1/math.sqrt(2))
        q_control_3 = MD_control_3 / (SED_control_3 * 1/math.sqrt(2))
        q_control_4 = MD_control_4 / (SED_control_4 * 1/math.sqrt(2))

        q_LPS_1 = MD_LPS_1 / (SED_LPS_1 * 1/math.sqrt(2))
        q_LPS_2 = MD_LPS_2 / (SED_LPS_2 * 1/math.sqrt(2))
        q_LPS_3 = MD_LPS_3 / (SED_LPS_3 * 1/math.sqrt(2))
        q_LPS_4 = MD_LPS_4 / (SED_LPS_4 * 1/math.sqrt(2))

        #################################################

        p_control_1 = psturng (q_control_1, 5, df_within_control)
        p_control_2 = psturng (q_control_2, 5, df_within_control)
        p_control_3 = psturng (q_control_3, 5, df_within_control)
        p_control_4 = psturng (q_control_4, 5, df_within_control)

        p_LPS_1 = psturng (q_LPS_1, 5, df_within_LPS)
        p_LPS_2 = psturng (q_LPS_2, 5, df_within_LPS)
        p_LPS_3 = psturng (q_LPS_3, 5, df_within_LPS)
        p_LPS_4 = psturng (q_LPS_4, 5, df_within_LPS)

        #################################################

        print ("========================================================")
        print ("")
        print ("Your Tukey's HSD Results Are:")
        print ("")

        if HSD_BBS_ten > HSD_control:
            print ("BBS vs. 10:  *  | p =", p_control_1)
        else:
            print ("BBS vs. 10:  ns | p =", p_control_1)

        if HSD_BBS_thirty > HSD_control:
            print ("BBS vs. 30:  *  | p =", p_control_2)
        else:
            print ("BBS vs. 30:  ns | p =", p_control_2)

        if HSD_BBS_hundred > HSD_control:
            print ("BBS vs. 100: *  | p =", p_control_3)
        else:
            print ("BBS vs. 100: ns | p =", p_control_3)

        if HSD_BBS_threehundred > HSD_control:
            print ("BBS vs. 300: *  | p =", p_control_4)
        else:
            print ("BBS vs. 300: ns | p =", p_control_4)

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if HSD_BBS_ten_LPS > HSD_LPS:
            print ("BBS LPS vs. 10 LPS:  *  | p =", p_LPS_1)
        else:
            print ("BBS LPS vs. 10 LPS:  ns | p =", p_LPS_1)
            
        if HSD_BBS_thirty_LPS > HSD_LPS:
            print ("BBS LPS vs. 30 LPS:  *  | p =", p_LPS_2)
        else:
            print ("BBS LPS vs. 30 LPS:  ns | p =", p_LPS_2)

        if HSD_BBS_hundred_LPS > HSD_LPS:
            print ("BBS LPS vs. 100 LPS: *  | p =", p_LPS_3)
        else:
            print ("BBS LPS vs. 100 LPS: ns | p =", p_LPS_3)

        if HSD_BBS_threehundred_LPS > HSD_LPS:
            print ("BBS LPS vs. 300 LPS: *  | p =", p_LPS_4)
        else:
            print ("BBS LPS vs. 300 LPS: ns | p =", p_LPS_4)

        print ("")
        print ("========================================================")
        print ("")

        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")
        
    elif answer == "3":

        #Tukey-Kramer: if the MD is greater than the CD, significant

        from statsmodels.stats.libqsturng import psturng
        
        #crit_table = pandas.read_csv ("/Users/evan/Desktop/Python Tests/Q critical values.csv", header = None)
        crit_table = [37.082,	10.881,	7.502,	6.287,	5.673,	5.305,	5.06,	4.886,	4.755,	4.654,	4.574,	4.508,	4.453,	4.407,	4.367,	4.333,	4.303,	4.276,	4.253,	4.232,	4.213,	4.196,	4.18,	4.166,	4.153,	4.141,	4.13,	4.12,	4.111,	4.102,	4.094,	4.086,	4.079,	4.072,	4.066,	4.06,	4.054,	4.049,	4.044,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.039,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	4.008,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.977,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.947,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.917,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887,	3.887]

        #################################################
      
        df_BBS = len (BBS) - 1
        df_ten = len (ten) - 1
        df_thirty = len (thirty) - 1
        df_hundred = len (hundred) - 1
        df_threehundred = len (threehundred) - 1
        df_Mt = len (Mt) - 1

        df_BBS_LPS = len (BBS_LPS) - 1
        df_ten_LPS = len (ten_LPS) - 1
        df_thirty_LPS = len (thirty_LPS) - 1
        df_hundred_LPS = len (hundred_LPS) - 1
        df_threehundred_LPS = len (threehundred_LPS) - 1
        df_Mt_LPS = len (Mt_LPS) - 1

        #################################################
        
        if len (BBS) == 1:
            BBS_var = 0
        else:
            BBS_var = statistics.variance (BBS)

        if len(ten) == 1:
            ten_var = 0
        else:
            ten_var = statistics.variance (ten)

        if len(thirty) == 1:
            thirty_var = 0
        else:
            thirty_var = statistics.variance (thirty)

        if len(hundred) == 1:
            hundred_var = 0
        else:
            hundred_var = statistics.variance (hundred)

        if len (threehundred) == 1:
            threehundred_var = 0
        else:
            threehundred_var = statistics.variance (threehundred)

        if len (Mt) == 1:
            Mt_var = 0
        else:
            Mt_var = statistics.variance (Mt)

        ############################

        if len (BBS_LPS) == 1:
            BBS_LPS_var = 0
        else:
            BBS_LPS_var = statistics.variance (BBS_LPS)

        if len(ten_LPS) == 1:
            ten_LPS_var = 0
        else:
            ten_LPS_var = statistics.variance (ten_LPS)

        if len(thirty_LPS) == 1:
            thirty_LPS_var = 0
        else:
            thirty_LPS_var = statistics.variance (thirty_LPS)

        if len(hundred_LPS) == 1:
            hundred_LPS_var = 0
        else:
            hundred_LPS_var = statistics.variance (hundred_LPS)

        if len (threehundred_LPS) == 1:
            threehundred_LPS_var = 0
        else:
            threehundred_LPS_var = statistics.variance (threehundred_LPS)

        if len (Mt_LPS) == 1:
            Mt_LPS_var = 0
        else:
            Mt_LPS_var = statistics.variance (Mt_LPS)

        #################################################

        SS_within_control = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var)
        df_within_control = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred
        MSE_within_control = SS_within_control / df_within_control
        

        SS_within_LPS = (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var)
        df_within_LPS = df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS
        MSE_within_LPS = SS_within_LPS / df_within_LPS

        #################################################

        #HSD_control = crit_table.iloc[int (float(df_within_control)-1), 0]
        #HSD_LPS = crit_table.iloc[int (float(df_within_LPS)-1), 0]

        HSD_control = crit_table[int (float(df_within_control)-1)]
        HSD_LPS = crit_table[int (float(df_within_LPS)-1)]

        #################################################

        CD_control_1 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(ten)) ))
        CD_control_2 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(thirty)) ))
        CD_control_3 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(hundred)) ))
        CD_control_4 = HSD_control * math.sqrt ( (MSE_within_control / 2) * ( (1/len(BBS)) + (1/len(threehundred)) ))

        CD_LPS_1 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        CD_LPS_2 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        CD_LPS_3 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        CD_LPS_4 = HSD_LPS * math.sqrt ( (MSE_within_LPS / 2) * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))
        
        #################################################

        MD_control_1 = abs (((sum (BBS) / len (BBS)) - (sum (ten) / len (ten))))
        MD_control_2 = abs (((sum (BBS) / len (BBS)) - (sum (thirty) / len (thirty))))
        MD_control_3 = abs (((sum (BBS) / len (BBS)) - (sum (hundred) / len (hundred))))
        MD_control_4 = abs (((sum (BBS) / len (BBS)) - (sum (threehundred) / len (threehundred))))

        MD_LPS_1 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (ten_LPS) / len (ten_LPS))))
        MD_LPS_2 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (thirty_LPS) / len (thirty_LPS))))
        MD_LPS_3 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (hundred_LPS) / len (hundred_LPS))))
        MD_LPS_4 = abs (((sum (BBS_LPS) / len (BBS_LPS)) - (sum (threehundred_LPS) / len (threehundred_LPS))))

        #################################################

        SED_control_1 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(ten)) ))
        SED_control_2 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(thirty)) ))
        SED_control_3 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(hundred)) ))
        SED_control_4 = math.sqrt (MSE_within_control * ( (1/len(BBS)) + (1/len(threehundred)) ))

        SED_LPS_1 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(ten_LPS)) ))
        SED_LPS_2 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(thirty_LPS)) ))
        SED_LPS_3 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(hundred_LPS)) ))
        SED_LPS_4 = math.sqrt (MSE_within_LPS * ( (1/len(BBS_LPS)) + (1/len(threehundred_LPS)) ))

        #################################################

        q_control_1 = MD_control_1 / (SED_control_1 * 1/math.sqrt(2))
        q_control_2 = MD_control_2 / (SED_control_2 * 1/math.sqrt(2))
        q_control_3 = MD_control_3 / (SED_control_3 * 1/math.sqrt(2))
        q_control_4 = MD_control_4 / (SED_control_4 * 1/math.sqrt(2))

        q_LPS_1 = MD_LPS_1 / (SED_LPS_1 * 1/math.sqrt(2))
        q_LPS_2 = MD_LPS_2 / (SED_LPS_2 * 1/math.sqrt(2))
        q_LPS_3 = MD_LPS_3 / (SED_LPS_3 * 1/math.sqrt(2))
        q_LPS_4 = MD_LPS_4 / (SED_LPS_4 * 1/math.sqrt(2))

        #################################################

        p_control_1 = psturng (q_control_1, 5, df_within_control)
        p_control_2 = psturng (q_control_2, 5, df_within_control)
        p_control_3 = psturng (q_control_3, 5, df_within_control)
        p_control_4 = psturng (q_control_4, 5, df_within_control)

        p_LPS_1 = psturng (q_LPS_1, 5, df_within_LPS)
        p_LPS_2 = psturng (q_LPS_2, 5, df_within_LPS)
        p_LPS_3 = psturng (q_LPS_3, 5, df_within_LPS)
        p_LPS_4 = psturng (q_LPS_4, 5, df_within_LPS)

        #################################################

        print ("========================================================")
        print ("")
        print ("Your Tukey-Kramer Results Are:")
        print ("")

        if MD_control_1 > CD_control_1:
            print ("BBS vs. 10:  *  ")
        else:
            print ("BBS vs. 10:  ns ")

        if MD_control_2 > CD_control_2:
            print ("BBS vs. 30:  *  ")
        else:
            print ("BBS vs. 30:  ns ")

        if MD_control_3 > CD_control_3:
            print ("BBS vs. 100: *  ")
        else:
            print ("BBS vs. 100: ns ")

        if MD_control_4 > CD_control_4:
            print ("BBS vs. 300: *  ")
        else:
            print ("BBS vs. 300: ns ")

        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")

        if MD_LPS_1 > CD_LPS_1:
            print ("BBS LPS vs. 10 LPS:  *  ")
        else:
            print ("BBS LPS vs. 10 LPS:  ns ")
            
        if MD_LPS_2 > CD_LPS_2:
            print ("BBS LPS vs. 30 LPS:  *  ")
        else:
            print ("BBS LPS vs. 30 LPS:  ns ")

        if MD_LPS_3 > CD_LPS_3:
            print ("BBS LPS vs. 100 LPS: *  ")
        else:
            print ("BBS LPS vs. 100 LPS: ns ")

        if MD_LPS_4 > CD_LPS_4:
            print ("BBS LPS vs. 300 LPS: *  ")
        else:
            print ("BBS LPS vs. 300 LPS: ns ")

        print ("")
        print ("*** p values are coming in a later version ***")

        print ("")
        print ("========================================================")
        print ("")
        
        try:
            shell = sys.stdout.shell
        except AttributeError:
            raise RuntimeError("you must run this program in IDLE")

        shell.write("Thank ", "COMMENT")
        shell.write("You ","KEYWORD")
        shell.write("For ","STRING")
        shell.write("Using ","DEFINITION")
        shell.write("The ","BUILTIN")
        shell.write("qPCR ","COMMENT")
        shell.write("Analyzer ","KEYWORD")
        shell.write("3000!","STRING")
        print ("")
        print ("")
        shell.write("-EH","BUILTIN")
        print ("")
        print ("")  

    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
            
    return;

def One_Way_ANOVA (n, dayta, LPS_factor):

    flex = n
    import scipy.stats

    if LPS_factor == False:
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        count = 0
        while count < (len(dayta)):
            BBS.append (dayta[count])
            ten.append (dayta[count + 1])
            thirty.append (dayta[count + 2])
            hundred.append (dayta[count + 3])
            threehundred.append (dayta[count + 4])
            Mt.append (dayta[count + 5])

            BBS_LPS.append (dayta[count + 6])
            ten_LPS.append (dayta[count + 7])
            thirty_LPS.append (dayta[count + 8])
            hundred_LPS.append (dayta[count + 9])
            threehundred_LPS.append (dayta[count + 10])
            Mt_LPS.append (dayta[count + 11])

            count += 12

        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)

        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)

        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)

        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)

        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)

        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)

        ################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)

        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)

        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)

        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)

        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)

        control_results = scipy.stats.f_oneway (New_BBS, New_ten, New_thirty, New_hundred, New_threehundred)
        LPS_results = scipy.stats.f_oneway (New_BBS_LPS, New_ten_LPS, New_thirty_LPS, New_hundred_LPS, New_threehundred_LPS)

        control_t_results = scipy.stats.ttest_ind (New_BBS, New_Mt, equal_var = True)
        LPS_t_results = scipy.stats.ttest_ind (New_BBS_LPS, New_Mt_LPS, equal_var = True)

        print ("========================================================")
        print ("")
        print ("Your One-Way ANOVA Values (BBS to 300) Are: ")
        print ("")
        print ("Control F value, P value:", control_results)
        print ("")
        print ("LPS F value, P value:", LPS_results)
        print ("")
        
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("")
        print ("Your T Test Values (BBS vs. M. tub) Are: ")
        print ("")
        print ("Control t-statistic, P value:", control_t_results)
        print ("")
        print ("LPS t-statistic, P value:", LPS_t_results)
        print ("")
        print ("========================================================")
        
        Next = input ("Run Post-Hoc Test? <Y/N> ")

        if Next == "Y" or Next == "y" or Next == "Yes" or Next == "yes" or Next == "YES":
            PCR_posthoc_test (N = n, BBS = New_BBS, ten = New_ten, thirty = New_thirty, hundred = New_hundred,
                         threehundred = New_threehundred, Mt = New_Mt, BBS_LPS = New_BBS_LPS, ten_LPS = New_ten_LPS,
                         thirty_LPS = New_thirty_LPS, hundred_LPS = New_hundred_LPS, threehundred_LPS =
                         New_threehundred_LPS, Mt_LPS = New_Mt_LPS)

        elif Next == "N" or Next == "n" or Next == "No" or Next == "no" or Next == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    if LPS_factor == True:
        from statsmodels.formula.api import ols
        from statsmodels.stats.anova import anova_lm
        import statistics
        from scipy import stats
        
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        BBS_LPS2 = []
        ten_LPS2 = []
        thirty_LPS2 = []
        hundred_LPS2 = []
        threehundred_LPS2 = []
        Mt_LPS2 = []

        BBS_LPS3 = []
        ten_LPS3 = []
        thirty_LPS3 = []
        hundred_LPS3 = []
        threehundred_LPS3 = []
        Mt_LPS3 = []

        count = 0
        while count < (len(dayta)):
            BBS.append (dayta[count])
            ten.append (dayta[count + 1])
            thirty.append (dayta[count + 2])
            hundred.append (dayta[count + 3])
            threehundred.append (dayta[count + 4])
            Mt.append (dayta[count + 5])

            BBS_LPS.append (dayta[count + 6])
            ten_LPS.append (dayta[count + 7])
            thirty_LPS.append (dayta[count + 8])
            hundred_LPS.append (dayta[count + 9])
            threehundred_LPS.append (dayta[count + 10])
            Mt_LPS.append (dayta[count + 11])
            
            BBS_LPS2.append (dayta[count + 12])
            ten_LPS2.append (dayta[count + 13])
            thirty_LPS2.append (dayta[count + 14])
            hundred_LPS2.append (dayta[count + 15])
            threehundred_LPS2.append (dayta[count + 16])
            Mt_LPS2.append (dayta[count + 17])

            BBS_LPS3.append (dayta[count + 18])
            ten_LPS3.append (dayta[count + 19])
            thirty_LPS3.append (dayta[count + 20])
            hundred_LPS3.append (dayta[count + 21])
            threehundred_LPS3.append (dayta[count + 22])
            Mt_LPS3.append (dayta[count + 23])

            count += 24

##        Treatments = ( (["BBS"]*int(n)) + (["10"]*int(n)) + (["30"]*int(n)) + (["100"]*int(n)) + (["300"]*int(n)) )
##        Control = BBS + ten + thirty + hundred + threehundred
##        LPS_1 = BBS_LPS + ten_LPS + thirty_LPS + hundred_LPS + threehundred_LPS 
##        LPS_2 = BBS_LPS2 + ten_LPS2 + thirty_LPS2 + hundred_LPS2 + threehundred_LPS2 
##        LPS_3 = BBS_LPS3 + ten_LPS3 + thirty_LPS3 + hundred_LPS3 + threehundred_LPS3 
##        
##        data = {"Bacteria_Dose": Treatments,
##                "Control": Control,
##                "LPS_1": LPS_1,
##                "LPS_2": LPS_2,
##                "LPS_3": LPS_3}
##
##        df = pandas.DataFrame (data, columns = ["Bacteria_Dose", "Control","LPS_1", "LPS_2", "LPS_3"])
##        df_melt = pandas.melt(df, id_vars = ["Bacteria_Dose"], value_vars = ["Control","LPS_1", "LPS_2", "LPS_3"])
##        df_melt.columns = ["Bacteria_Dose", "LPS_Dose", "Log2_Value"]
##        print (df_melt)
##        model = ols("Log2_Value ~ (Bacteria_Dose) + (LPS_Dose) + (Bacteria_Dose):(LPS_Dose)", data = df_melt).fit()
##        aov_table = anova_lm(model, typ = 2)
        
        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)
        sBBS_mean = (sum (New_BBS) / len (New_BBS))
        
        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)
        sten_mean = (sum (New_ten) / len (New_ten))
        
        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)
        sthirty_mean = (sum (New_thirty) / len (New_thirty))
        
        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)
        shundred_mean = (sum (New_hundred) / len (New_hundred))
        
        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)
        sthreehundred_mean = (sum (New_threehundred) / len (New_threehundred))
        
        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)
        sMt_mean = (sum (New_Mt) / len (New_Mt))
        
        #################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)
        sBBS_LPS_mean = (sum (New_BBS_LPS) / len (New_BBS_LPS))
        
        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)
        sten_LPS_mean = (sum (New_ten_LPS) / len (New_ten_LPS))
        
        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)
        sthirty_LPS_mean = (sum (New_thirty_LPS) / len (New_thirty_LPS))
        
        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)
        shundred_LPS_mean = (sum (New_hundred_LPS) / len (New_hundred_LPS))
        
        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)
        sthreehundred_LPS_mean = (sum (New_threehundred_LPS) / len (New_threehundred_LPS))
        
        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)
        sMt_LPS_mean = (sum (New_Mt_LPS) / len (New_Mt_LPS))

        ###################################

        New_BBS_LPS2 = []
        for m in BBS_LPS2:
            if type(m) == str:
                continue
            New_BBS_LPS2.append (m)
        sBBS_LPS2_mean = (sum (New_BBS_LPS2) / len (New_BBS_LPS2))

        New_ten_LPS2 = []
        for n in ten_LPS2:
            if type(n) == str:
                continue
            New_ten_LPS2.append (n)
        sten_LPS2_mean = (sum (New_ten_LPS2) / len (New_ten_LPS2))

        New_thirty_LPS2 = []
        for o in thirty_LPS2:
            if type(o) == str:
                continue
            New_thirty_LPS2.append (o)
        sthirty_LPS2_mean = (sum (New_thirty_LPS2) / len (New_thirty_LPS2))

        New_hundred_LPS2 = []
        for p in hundred_LPS2:
            if type(p) == str:
                continue
            New_hundred_LPS2.append (p)
        shundred_LPS2_mean = (sum (New_hundred_LPS2) / len (New_hundred_LPS2))

        New_threehundred_LPS2 = []
        for q in threehundred_LPS2:
            if type(q) == str:
                continue
            New_threehundred_LPS2.append (q)
        sthreehundred_LPS2_mean = (sum (New_threehundred_LPS2) / len (New_threehundred_LPS2))

        New_Mt_LPS2 = []
        for r in Mt_LPS2:
            if type(r) == str:
                continue
            New_Mt_LPS2.append (r)
        sMt_LPS2_mean = (sum (New_Mt_LPS2) / len (New_Mt_LPS2))

        ###################################
        
        New_BBS_LPS3 = []
        for s in BBS_LPS3:
            if type(s) == str:
                continue
            New_BBS_LPS3.append (s)
        sBBS_LPS3_mean = (sum (New_BBS_LPS3) / len (New_BBS_LPS3))

        New_ten_LPS3 = []
        for t in ten_LPS3:
            if type(t) == str:
                continue
            New_ten_LPS3.append (t)
        sten_LPS3_mean = (sum (New_ten_LPS3) / len (New_ten_LPS3))

        New_thirty_LPS3 = []
        for u in thirty_LPS3:
            if type(u) == str:
                continue
            New_thirty_LPS3.append (u)
        sthirty_LPS3_mean = (sum (New_thirty_LPS3) / len (New_thirty_LPS3))

        New_hundred_LPS3 = []
        for v in hundred_LPS3:
            if type(v) == str:
                continue
            New_hundred_LPS3.append (v)
        shundred_LPS3_mean = (sum (New_hundred_LPS3) / len (New_hundred_LPS3))

        New_threehundred_LPS3 = []
        for w in threehundred_LPS3:
            if type(w) == str:
                continue
            New_threehundred_LPS3.append (w)
        sthreehundred_LPS3_mean = (sum (New_threehundred_LPS3) / len (New_threehundred_LPS3))

        New_Mt_LPS3 = []
        for x in Mt_LPS3:
            if type(x) == str:
                continue
            New_Mt_LPS3.append (x)
        sMt_LPS3_mean = (sum (New_Mt_LPS3) / len (New_Mt_LPS3))

        ###################################

        #Two-WAy ANOVA

        try:
            BBS_var = statistics.variance (New_BBS)
        except:
            BBS_var = 0
            print ("Control BBS variance has been assumed to be 0, since there isn't at least 2 samples")
        try:
            ten_var = statistics.variance (New_ten)
        except:
            ten_var = 0
            print ("Control 10 variance has been assumed to be 0, since there isn't at least 2 samples")

        thirty_var = statistics.variance (New_thirty)
        hundred_var = statistics.variance (New_hundred)
        threehundred_var = statistics.variance (New_threehundred)

        BBS_LPS_var = statistics.variance (New_BBS_LPS)
        ten_LPS_var = statistics.variance (New_ten_LPS)
        thirty_LPS_var = statistics.variance (New_thirty_LPS)
        hundred_LPS_var = statistics.variance (New_hundred_LPS)
        threehundred_LPS_var = statistics.variance (New_threehundred_LPS)

        BBS_LPS2_var = statistics.variance (New_BBS_LPS2)
        ten_LPS2_var = statistics.variance (New_ten_LPS2)
        thirty_LPS2_var = statistics.variance (New_thirty_LPS2)
        hundred_LPS2_var = statistics.variance (New_hundred_LPS2)
        threehundred_LPS2_var = statistics.variance (New_threehundred_LPS2)

        BBS_LPS3_var = statistics.variance (New_BBS_LPS3)
        ten_LPS3_var = statistics.variance (New_ten_LPS3)
        thirty_LPS3_var = statistics.variance (New_thirty_LPS3)
        hundred_LPS3_var = statistics.variance (New_hundred_LPS3)
        threehundred_LPS3_var = statistics.variance (New_threehundred_LPS3)

        Control_mean = statistics.mean (New_BBS + New_ten + New_thirty + New_hundred + New_threehundred)
        LPS_mean = statistics.mean (New_BBS_LPS + New_ten_LPS + New_thirty_LPS + New_hundred_LPS + New_threehundred_LPS)
        LPS2_mean = statistics.mean (New_BBS_LPS2 + New_ten_LPS2 + New_thirty_LPS2 + New_hundred_LPS2 + New_threehundred_LPS2)
        LPS3_mean = statistics.mean (New_BBS_LPS3 + New_ten_LPS3 + New_thirty_LPS3 + New_hundred_LPS3 + New_threehundred_LPS3)
        LPS_dose_mean = statistics.mean(New_BBS + New_ten + New_thirty + New_hundred + New_threehundred + New_BBS_LPS + New_ten_LPS + New_thirty_LPS + New_hundred_LPS + New_threehundred_LPS + New_BBS_LPS2 + New_ten_LPS2 + New_thirty_LPS2 + New_hundred_LPS2 + New_threehundred_LPS2 + New_BBS_LPS3 + New_ten_LPS3 + New_thirty_LPS3 + New_hundred_LPS3 + New_threehundred_LPS3) 

        BBS_mean = statistics.mean (New_BBS + New_BBS_LPS + New_BBS_LPS2 + New_BBS_LPS3)
        ten_mean = statistics.mean (New_ten + New_ten_LPS + New_ten_LPS2 + New_ten_LPS3)
        thirty_mean = statistics.mean (New_thirty + New_thirty_LPS + New_thirty_LPS2 + New_thirty_LPS3)
        hundred_mean = statistics.mean (New_hundred + New_hundred_LPS + New_hundred_LPS2 + New_hundred_LPS3)
        threehundred_mean = statistics.mean (New_threehundred + New_threehundred_LPS + New_threehundred_LPS2 + New_threehundred_LPS3)
        Bacteria_dose_mean = LPS_dose_mean
        
        df_BBS = len (BBS) - 1
        df_ten = len (ten) - 1
        df_thirty = len (thirty) - 1
        df_hundred = len (hundred) - 1
        df_threehundred = len (threehundred) - 1
        df_Mt = len (Mt) - 1

        df_BBS_LPS = len (BBS_LPS) - 1
        df_ten_LPS = len (ten_LPS) - 1
        df_thirty_LPS = len (thirty_LPS) - 1
        df_hundred_LPS = len (hundred_LPS) - 1
        df_threehundred_LPS = len (threehundred_LPS) - 1
        df_Mt_LPS = len (Mt_LPS) - 1

        df_BBS_LPS2 = len (BBS_LPS2) - 1
        df_ten_LPS2 = len (ten_LPS2) - 1
        df_thirty_LPS2 = len (thirty_LPS2) - 1
        df_hundred_LPS2 = len (hundred_LPS2) - 1
        df_threehundred_LPS2 = len (threehundred_LPS2) - 1
        df_Mt_LPS2 = len (Mt_LPS2) - 1

        df_BBS_LPS3 = len (BBS_LPS3) - 1
        df_ten_LPS3 = len (ten_LPS3) - 1
        df_thirty_LPS3 = len (thirty_LPS3) - 1
        df_hundred_LPS3 = len (hundred_LPS3) - 1
        df_threehundred_LPS3 = len (threehundred_LPS3) - 1
        df_Mt_LPS3 = len (Mt_LPS3) - 1

        SS_within = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var) + (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var) + (df_BBS_LPS2 * BBS_LPS2_var) + (df_ten_LPS2 * ten_LPS2_var) + (df_thirty_LPS2 * thirty_LPS2_var) + (df_hundred_LPS2 * hundred_LPS2_var) + (df_threehundred_LPS2 * threehundred_LPS2_var) + (df_BBS_LPS3 * BBS_LPS3_var) + (df_ten_LPS3 * ten_LPS3_var) + (df_thirty_LPS3 * thirty_LPS3_var) + (df_hundred_LPS3 * hundred_LPS3_var) + (df_threehundred_LPS3 * threehundred_LPS3_var) 
        df_within = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred + df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS + df_BBS_LPS2 + df_ten_LPS2 + df_thirty_LPS2 + df_hundred_LPS2 + df_threehundred_LPS2 + df_BBS_LPS3 + df_ten_LPS3 + df_thirty_LPS3 + df_hundred_LPS3 + df_threehundred_LPS3
        MSE_within = SS_within / df_within

        ###Flex?
        SS_LPS_dose = int(flex) * 5 * (((Control_mean - LPS_dose_mean)**2) + ((LPS_mean - LPS_dose_mean)**2) + ((LPS2_mean - LPS_dose_mean)**2) + ((LPS3_mean - LPS_dose_mean)**2))
        df_LPS_dose = 3
        MSE_LPS_dose = SS_LPS_dose/df_LPS_dose

        ###Flex?
        SS_bacteria_dose = int(flex) * 4 * (((BBS_mean - Bacteria_dose_mean)**2) + ((ten_mean - Bacteria_dose_mean)**2) + ((thirty_mean - Bacteria_dose_mean)**2) + ((hundred_mean - Bacteria_dose_mean)**2) + ((threehundred_mean - Bacteria_dose_mean)**2))
        df_bacteria_dose = 4
        MSE_bacteria_dose = SS_bacteria_dose/df_bacteria_dose

        #Yij = replicate mean; Yi = bacteria dose mean; Yj = LPS dose mean; Y = grand mean
        ###Flex?
        SS_interaction = int(flex) * (((sBBS_mean-BBS_mean-Control_mean+LPS_dose_mean)**2) + ((sten_mean-ten_mean-Control_mean+LPS_dose_mean)**2) + ((sthirty_mean-thirty_mean-Control_mean+LPS_dose_mean)**2) + ((shundred_mean-hundred_mean-Control_mean+LPS_dose_mean)**2) + ((sthreehundred_mean-threehundred_mean-Control_mean+LPS_dose_mean)**2) + ((sBBS_LPS_mean-BBS_mean-LPS_mean+LPS_dose_mean)**2) + ((sten_LPS_mean-ten_mean-LPS_mean+LPS_dose_mean)**2) + ((sthirty_LPS_mean-thirty_mean-LPS_mean+LPS_dose_mean)**2) + ((shundred_LPS_mean-hundred_mean-LPS_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS_mean-threehundred_mean-LPS_mean+LPS_dose_mean)**2) + ((sBBS_LPS2_mean-BBS_mean-LPS2_mean+LPS_dose_mean)**2) + ((sten_LPS2_mean-ten_mean-LPS2_mean+LPS_dose_mean)**2) + ((sthirty_LPS2_mean-thirty_mean-LPS2_mean+LPS_dose_mean)**2) + ((shundred_LPS2_mean-hundred_mean-LPS2_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS2_mean-threehundred_mean-LPS2_mean+LPS_dose_mean)**2) + ((sBBS_LPS3_mean-BBS_mean-LPS3_mean+LPS_dose_mean)**2) + ((sten_LPS3_mean-ten_mean-LPS3_mean+LPS_dose_mean)**2) + ((sthirty_LPS3_mean-thirty_mean-LPS3_mean+LPS_dose_mean)**2) + ((shundred_LPS3_mean-hundred_mean-LPS3_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS3_mean-threehundred_mean-LPS3_mean+LPS_dose_mean)**2))
        df_interaction = 12
        MSE_interaction = SS_interaction/df_interaction

        F_LPS_dose = MSE_LPS_dose / MSE_within
        F_bacteria_dose = MSE_bacteria_dose / MSE_within
        F_interaction = MSE_interaction / MSE_within

        p_LPS_dose = stats.f.sf(F_LPS_dose, df_LPS_dose, df_within)
        p_bacteria_dose = stats.f.sf(F_bacteria_dose, df_bacteria_dose, df_within)
        p_interaction = stats.f.sf (F_interaction, df_interaction, df_within)
        
        control_results = scipy.stats.f_oneway (New_BBS, New_ten, New_thirty, New_hundred, New_threehundred)
        LPS_results = scipy.stats.f_oneway (New_BBS_LPS, New_ten_LPS, New_thirty_LPS, New_hundred_LPS, New_threehundred_LPS)
        LPS2_results = scipy.stats.f_oneway (New_BBS_LPS2,New_ten_LPS2,New_thirty_LPS2,New_hundred_LPS2,New_threehundred_LPS2)
        LPS3_results = scipy.stats.f_oneway (New_BBS_LPS3,New_ten_LPS3,New_thirty_LPS3,New_hundred_LPS3,New_threehundred_LPS3)
        
        control_t_results = scipy.stats.ttest_ind (New_BBS, New_Mt, equal_var = True)
        LPS_t_results = scipy.stats.ttest_ind (New_BBS_LPS, New_Mt_LPS, equal_var = True)
        LPS2_t_results = scipy.stats.ttest_ind (New_BBS_LPS2, New_Mt_LPS2, equal_var = True)
        LPS3_t_results = scipy.stats.ttest_ind (New_BBS_LPS3, New_Mt_LPS3, equal_var = True)

        print ("========================================================")
        
        print ("")
        print ("Your One-Way ANOVA Values (BBS to 300) Are: ")
        print ("")
        print ("Control F value, P value:", control_results)
        print ("")
        print ("1st LPS Dose F value, P value:", LPS_results)
        print ("")
        print ("2nd LPS Dose F value, P value:", LPS2_results)
        print ("")
        print ("3rd LPS Dose F value, P value:", LPS3_results)
        print ("")
##        
##        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

##        print("")
##        print ("Your Two-Way ANOVA Values (BBS to 300) Are: ")
##        print ("")
##        print ("LPS_Dose F and P value: ", "F =", F_LPS_dose, "P =", p_LPS_dose)
##        print ("")
##        print ("Bacteria_Dose F and P value:", "F =", F_bacteria_dose, "P =", p_bacteria_dose)
##        print ("")
##        print ("LPS_Dose and Bacteria_Dose Interaction F and P value:", "F =", F_interaction, "P =", p_interaction)
##        print ("")
        
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("")
        print ("Your T Test Values (BBS vs. M. tub) Are: ")
        print ("")
        print ("Control t-statistic, P value:", control_t_results)
        print ("")
        print ("1st LPS Dose t-statistic, P value:", LPS_t_results)
        print ("")
        print ("2nd LPS Dose t-statistic, P value:", LPS2_t_results)
        print ("")
        print ("3rd LPS Dose t-statistic, P value:", LPS3_t_results)
        print ("")
        print ("========================================================")
        
        Next = input ("Run Post-Hoc Test? <Y/N> ")

        if Next == "Y" or Next == "y" or Next == "Yes" or Next == "yes" or Next == "YES":
            LPS_PCR_posthoc_test (N = flex, BBS = New_BBS, ten = New_ten, thirty = New_thirty, hundred = New_hundred,
                         threehundred = New_threehundred, Mt = New_Mt, BBS_LPS = New_BBS_LPS, ten_LPS = New_ten_LPS,
                         thirty_LPS = New_thirty_LPS, hundred_LPS = New_hundred_LPS, threehundred_LPS =
                         New_threehundred_LPS, Mt_LPS = New_Mt_LPS, BBS_LPS2 = New_BBS_LPS2, ten_LPS2 =
                         New_ten_LPS2, thirty_LPS2 = New_thirty_LPS2, hundred_LPS2 = New_hundred_LPS2,
                         threehundred_LPS2 = New_threehundred_LPS2, Mt_LPS2 = New_Mt_LPS2, BBS_LPS3 =
                         New_BBS_LPS3, ten_LPS3 = New_ten_LPS3, thirty_LPS3 = New_thirty_LPS3, hundred_LPS3 =
                         New_hundred_LPS3, threehundred_LPS3 = New_threehundred_LPS3, Mt_LPS3 = New_Mt_LPS3)

        elif Next == "N" or Next == "n" or Next == "No" or Next == "no" or Next == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
    
    return;

def PCR_grapher (N, Data, LPS_factor):
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import statistics
    
    if LPS_factor == False:
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        count = 0
        while count < (len(Data)):
            BBS.append (Data[count])
            ten.append (Data[count + 1])
            thirty.append (Data[count + 2])
            hundred.append (Data[count + 3])
            threehundred.append (Data[count + 4])
            Mt.append (Data[count + 5])

            BBS_LPS.append (Data[count + 6])
            ten_LPS.append (Data[count + 7])
            thirty_LPS.append (Data[count + 8])
            hundred_LPS.append (Data[count + 9])
            threehundred_LPS.append (Data[count + 10])
            Mt_LPS.append (Data[count + 11])

            count += 12
            
        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)
        BBS_mean = (sum (New_BBS) / len (New_BBS))
        
        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)
        ten_mean = (sum (New_ten) / len (New_ten))

        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)
        thirty_mean = (sum (New_thirty) / len (New_thirty))

        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)
        hundred_mean = (sum (New_hundred) / len (New_hundred))

        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)
        threehundred_mean = (sum (New_threehundred) / len (New_threehundred))

        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)
        Mt_mean = (sum (New_Mt) / len (New_Mt))

        #################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)
        BBS_LPS_mean = (sum (New_BBS_LPS) / len (New_BBS_LPS))

        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)
        ten_LPS_mean = (sum (New_ten_LPS) / len (New_ten_LPS))

        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)
        thirty_LPS_mean = (sum (New_thirty_LPS) / len (New_thirty_LPS))

        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)
        hundred_LPS_mean = (sum (New_hundred_LPS) / len (New_hundred_LPS))

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)
        threehundred_LPS_mean = (sum (New_threehundred_LPS) / len (New_threehundred_LPS))

        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)
        Mt_LPS_mean = (sum (New_Mt_LPS) / len (New_Mt_LPS))

        ###################################

        BBS_SEM = (statistics.stdev (New_BBS) / math.sqrt(len (New_BBS)))
        ten_SEM = (statistics.stdev (New_ten) / math.sqrt(len (New_ten)))
        thirty_SEM = (statistics.stdev (New_thirty) / math.sqrt(len (New_thirty)))
        hundred_SEM = (statistics.stdev (New_hundred) / math.sqrt(len (New_hundred)))
        threehundred_SEM = (statistics.stdev (New_threehundred) / math.sqrt(len (New_threehundred)))
        Mt_SEM = (statistics.stdev (New_Mt) / math.sqrt(len (New_Mt)))

        BBS_LPS_SEM = (statistics.stdev (New_BBS_LPS) / math.sqrt(len (New_BBS_LPS)))
        ten_LPS_SEM = (statistics.stdev (New_ten_LPS) / math.sqrt(len (New_ten_LPS)))
        thirty_LPS_SEM = (statistics.stdev (New_thirty_LPS) / math.sqrt(len (New_thirty_LPS)))
        hundred_LPS_SEM = (statistics.stdev (New_hundred_LPS) / math.sqrt(len (New_hundred_LPS)))
        threehundred_LPS_SEM = (statistics.stdev (New_threehundred_LPS) / math.sqrt(len (New_threehundred_LPS)))
        Mt_LPS_SEM = (statistics.stdev (New_Mt_LPS) / math.sqrt(len (New_Mt_LPS)))
        
        print ("========================================================")
        print ("")
        print ("Your Average Log2 Values +/- SEM Are: ")
        print ("")
        print ("Control:")
        print (BBS_mean, "+/-", BBS_SEM)
        print (ten_mean, "+/-", ten_SEM)
        print (thirty_mean, "+/-", thirty_SEM)
        print (hundred_mean, "+/-", hundred_SEM)
        print (threehundred_mean, "+/-", threehundred_SEM)
        print (Mt_mean, "+/-", Mt_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("LPS:")
        print (BBS_LPS_mean, "+/-", BBS_LPS_SEM)
        print (ten_LPS_mean, "+/-", ten_LPS_SEM)
        print (thirty_LPS_mean, "+/-", thirty_LPS_SEM)
        print (hundred_LPS_mean, "+/-", hundred_LPS_SEM)
        print (threehundred_LPS_mean, "+/-", threehundred_LPS_SEM)
        print (Mt_LPS_mean, "+/-", Mt_LPS_SEM)
        print ("")
        print ("========================================================")
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = N, dayta = Data, LPS_factor = False)

        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        #control_Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean]
        #LPS_Data = [BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #control_Labels = ("BBS", "10", "30", "100", "300", "M. tub")
        #LPS_Labels = ("LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(control_Labels))
        #x_pos2 = np.arange(len(LPS_Labels))

        #control_error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM]
        #LPS_error = [BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]


        #data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]
        
        #labels = ("BBS", "10", "30", "100", "300", "M. tub", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(labels))

        #plt.bar (x_pos, data, yerr = error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #plt.xticks (x_pos, labels)
        #plt.show()
        
        #fig, ax = plt.subplots()
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, ecolor = "black", capsize = 10)
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xlabel ("Treatment")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #plt.show()

        #fig2, ax2 = plt.subplots()
        #ax2.bar (x_pos, LPS_Data, yerr = LPS_error, align = "center", alpha = 0.75, ecolor = "black", capsize = 10)
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos)
        #ax2.set_xticklabels (LPS_Labels)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #plt.show()

        #fig, (ax, ax2) = plt.subplots (nrows = 2)
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        

        #ax2.bar (x_pos2, LPS_Data, yerr = LPS_error, align = "center", alpha = 1, color = "red", ecolor = "black", capsize = 10)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos2)
        #ax2.set_xticklabels (LPS_Labels)
        #plt.show ()
        
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        
        Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, 0, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        Labels = ("BBS", "10", "30", "100", "300", "TB", "", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS TB")
        
        x_pos = np.arange(len(Labels))
       
        Error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, 0, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]

        fig, ax = plt.subplots (figsize = (9.25, 7.25))
        barlist = plt.bar (x_pos, Data, align = "center", yerr = Error, ecolor = "black", capsize = 5, alpha = 1)
        
        barlist [0].set_color ("blue")
        barlist [1].set_color ("blue")
        barlist [2].set_color ("blue")
        barlist [3].set_color ("blue")
        barlist [4].set_color ("blue")
        barlist [5].set_color ("blue")

        barlist [7].set_color ("red")
        barlist [8].set_color ("red")
        barlist [9].set_color ("red")
        barlist [10].set_color ("red")
        barlist [11].set_color ("red")
        barlist [12].set_color ("red")

        plt.xticks (x_pos, Labels, rotation = 45)
        plt.ylabel ("Relative Gene Expression (Log2 Fold Change)")
        
        plt.show()
        
        #plt.bar (x_pos, control_Data, align = "center", alpha = 1)
        #plt.xticks (x_pos, control_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("Control Graph")
        #plt.show()
        
        #plt.bar (x_pos_2, LPS_Data, align = "center", alpha = 0.75)
        #plt.xticks (x_pos_2, LPS_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("LPS Graph")
        #plt.show ()

    if LPS_factor == True:
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        BBS_LPS2 = []
        ten_LPS2 = []
        thirty_LPS2 = []
        hundred_LPS2 = []
        threehundred_LPS2 = []
        Mt_LPS2 = []

        BBS_LPS3 = []
        ten_LPS3 = []
        thirty_LPS3 = []
        hundred_LPS3 = []
        threehundred_LPS3 = []
        Mt_LPS3 = []

        count = 0
        while count < (len(Data)):
            BBS.append (Data[count])
            ten.append (Data[count + 1])
            thirty.append (Data[count + 2])
            hundred.append (Data[count + 3])
            threehundred.append (Data[count + 4])
            Mt.append (Data[count + 5])

            BBS_LPS.append (Data[count + 6])
            ten_LPS.append (Data[count + 7])
            thirty_LPS.append (Data[count + 8])
            hundred_LPS.append (Data[count + 9])
            threehundred_LPS.append (Data[count + 10])
            Mt_LPS.append (Data[count + 11])

            BBS_LPS2.append (Data[count + 12])
            ten_LPS2.append (Data[count + 13])
            thirty_LPS2.append (Data[count + 14])
            hundred_LPS2.append (Data[count + 15])
            threehundred_LPS2.append (Data[count + 16])
            Mt_LPS2.append (Data[count + 17])

            BBS_LPS3.append (Data[count + 18])
            ten_LPS3.append (Data[count + 19])
            thirty_LPS3.append (Data[count + 20])
            hundred_LPS3.append (Data[count + 21])
            threehundred_LPS3.append (Data[count + 22])
            Mt_LPS3.append (Data[count + 23])
            
            count += 24
            
        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)
        BBS_mean = (sum (New_BBS) / len (New_BBS))
        
        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)
        ten_mean = (sum (New_ten) / len (New_ten))

        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)
        thirty_mean = (sum (New_thirty) / len (New_thirty))

        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)
        hundred_mean = (sum (New_hundred) / len (New_hundred))

        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)
        threehundred_mean = (sum (New_threehundred) / len (New_threehundred))

        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)
        Mt_mean = (sum (New_Mt) / len (New_Mt))

        #################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)
        BBS_LPS_mean = (sum (New_BBS_LPS) / len (New_BBS_LPS))

        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)
        ten_LPS_mean = (sum (New_ten_LPS) / len (New_ten_LPS))

        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)
        thirty_LPS_mean = (sum (New_thirty_LPS) / len (New_thirty_LPS))

        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)
        hundred_LPS_mean = (sum (New_hundred_LPS) / len (New_hundred_LPS))

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)
        threehundred_LPS_mean = (sum (New_threehundred_LPS) / len (New_threehundred_LPS))

        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)
        Mt_LPS_mean = (sum (New_Mt_LPS) / len (New_Mt_LPS))

        ###################################

        New_BBS_LPS2 = []
        for m in BBS_LPS2:
            if type(m) == str:
                continue
            New_BBS_LPS2.append (m)
        BBS_LPS2_mean = (sum (New_BBS_LPS2) / len (New_BBS_LPS2))

        New_ten_LPS2 = []
        for n in ten_LPS2:
            if type(n) == str:
                continue
            New_ten_LPS2.append (n)
        ten_LPS2_mean = (sum (New_ten_LPS2) / len (New_ten_LPS2))

        New_thirty_LPS2 = []
        for o in thirty_LPS2:
            if type(o) == str:
                continue
            New_thirty_LPS2.append (o)
        thirty_LPS2_mean = (sum (New_thirty_LPS2) / len (New_thirty_LPS2))

        New_hundred_LPS2 = []
        for p in hundred_LPS2:
            if type(p) == str:
                continue
            New_hundred_LPS2.append (p)
        hundred_LPS2_mean = (sum (New_hundred_LPS2) / len (New_hundred_LPS2))

        New_threehundred_LPS2 = []
        for q in threehundred_LPS2:
            if type(q) == str:
                continue
            New_threehundred_LPS2.append (q)
        threehundred_LPS2_mean = (sum (New_threehundred_LPS2) / len (New_threehundred_LPS2))

        New_Mt_LPS2 = []
        for r in Mt_LPS2:
            if type(r) == str:
                continue
            New_Mt_LPS2.append (r)
        Mt_LPS2_mean = (sum (New_Mt_LPS2) / len (New_Mt_LPS2))

        ###################################
        
        New_BBS_LPS3 = []
        for s in BBS_LPS3:
            if type(s) == str:
                continue
            New_BBS_LPS3.append (s)
        BBS_LPS3_mean = (sum (New_BBS_LPS3) / len (New_BBS_LPS3))

        New_ten_LPS3 = []
        for t in ten_LPS3:
            if type(t) == str:
                continue
            New_ten_LPS3.append (t)
        ten_LPS3_mean = (sum (New_ten_LPS3) / len (New_ten_LPS3))

        New_thirty_LPS3 = []
        for u in thirty_LPS3:
            if type(u) == str:
                continue
            New_thirty_LPS3.append (u)
        thirty_LPS3_mean = (sum (New_thirty_LPS3) / len (New_thirty_LPS3))

        New_hundred_LPS3 = []
        for v in hundred_LPS3:
            if type(v) == str:
                continue
            New_hundred_LPS3.append (v)
        hundred_LPS3_mean = (sum (New_hundred_LPS3) / len (New_hundred_LPS3))

        New_threehundred_LPS3 = []
        for w in threehundred_LPS3:
            if type(w) == str:
                continue
            New_threehundred_LPS3.append (w)
        threehundred_LPS3_mean = (sum (New_threehundred_LPS3) / len (New_threehundred_LPS3))

        New_Mt_LPS3 = []
        for x in Mt_LPS3:
            if type(x) == str:
                continue
            New_Mt_LPS3.append (x)
        Mt_LPS3_mean = (sum (New_Mt_LPS3) / len (New_Mt_LPS3))

        ###################################
        
        BBS_SEM = (statistics.stdev (New_BBS) / math.sqrt(len (New_BBS)))
        ten_SEM = (statistics.stdev (New_ten) / math.sqrt(len (New_ten)))
        thirty_SEM = (statistics.stdev (New_thirty) / math.sqrt(len (New_thirty)))
        hundred_SEM = (statistics.stdev (New_hundred) / math.sqrt(len (New_hundred)))
        threehundred_SEM = (statistics.stdev (New_threehundred) / math.sqrt(len (New_threehundred)))
        Mt_SEM = (statistics.stdev (New_Mt) / math.sqrt(len (New_Mt)))

        BBS_LPS_SEM = (statistics.stdev (New_BBS_LPS) / math.sqrt(len (New_BBS_LPS)))
        ten_LPS_SEM = (statistics.stdev (New_ten_LPS) / math.sqrt(len (New_ten_LPS)))
        thirty_LPS_SEM = (statistics.stdev (New_thirty_LPS) / math.sqrt(len (New_thirty_LPS)))
        hundred_LPS_SEM = (statistics.stdev (New_hundred_LPS) / math.sqrt(len (New_hundred_LPS)))
        threehundred_LPS_SEM = (statistics.stdev (New_threehundred_LPS) / math.sqrt(len (New_threehundred_LPS)))
        Mt_LPS_SEM = (statistics.stdev (New_Mt_LPS) / math.sqrt(len (New_Mt_LPS)))

        BBS_LPS2_SEM = (statistics.stdev (New_BBS_LPS2) / math.sqrt(len(New_BBS_LPS2)))
        ten_LPS2_SEM = (statistics.stdev (New_ten_LPS2) / math.sqrt(len(New_ten_LPS2)))
        thirty_LPS2_SEM = (statistics.stdev (New_thirty_LPS2) / math.sqrt(len(New_thirty_LPS2)))
        hundred_LPS2_SEM = (statistics.stdev (New_hundred_LPS2) / math.sqrt(len(New_hundred_LPS2)))
        threehundred_LPS2_SEM = (statistics.stdev (New_threehundred_LPS2) / math.sqrt(len(New_threehundred_LPS2)))
        Mt_LPS2_SEM = (statistics.stdev (New_Mt_LPS2) / math.sqrt(len(New_Mt_LPS2)))

        BBS_LPS3_SEM = (statistics.stdev (New_BBS_LPS3) / math.sqrt(len(New_BBS_LPS3)))
        ten_LPS3_SEM = (statistics.stdev (New_ten_LPS3) / math.sqrt(len(New_ten_LPS3)))
        thirty_LPS3_SEM = (statistics.stdev (New_thirty_LPS3) / math.sqrt(len(New_thirty_LPS3)))
        hundred_LPS3_SEM = (statistics.stdev (New_hundred_LPS3) / math.sqrt(len(New_hundred_LPS3)))
        threehundred_LPS3_SEM = (statistics.stdev (New_threehundred_LPS3) / math.sqrt(len(New_threehundred_LPS3)))
        Mt_LPS3_SEM = (statistics.stdev (New_Mt_LPS3) / math.sqrt(len(New_Mt_LPS3)))
        
        print ("========================================================")
        print ("")
        print ("Your Average Log2 Values +/- SEM Are: ")
        print ("")
        print ("Control:")
        print (BBS_mean, "+/-", BBS_SEM)
        print (ten_mean, "+/-", ten_SEM)
        print (thirty_mean, "+/-", thirty_SEM)
        print (hundred_mean, "+/-", hundred_SEM)
        print (threehundred_mean, "+/-", threehundred_SEM)
        print (Mt_mean, "+/-", Mt_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("1st LPS Dose:")
        print (BBS_LPS_mean, "+/-", BBS_LPS_SEM)
        print (ten_LPS_mean, "+/-", ten_LPS_SEM)
        print (thirty_LPS_mean, "+/-", thirty_LPS_SEM)
        print (hundred_LPS_mean, "+/-", hundred_LPS_SEM)
        print (threehundred_LPS_mean, "+/-", threehundred_LPS_SEM)
        print (Mt_LPS_mean, "+/-", Mt_LPS_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("2nd LPS Dose:")
        print (BBS_LPS2_mean, "+/-", BBS_LPS2_SEM)
        print (ten_LPS2_mean, "+/-", ten_LPS2_SEM)
        print (thirty_LPS2_mean, "+/-", thirty_LPS2_SEM)
        print (hundred_LPS2_mean, "+/-", hundred_LPS2_SEM)
        print (threehundred_LPS2_mean, "+/-", threehundred_LPS2_SEM)
        print (Mt_LPS2_mean, "+/-", Mt_LPS2_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("3rd LPS Dose:")
        print (BBS_LPS3_mean, "+/-", BBS_LPS3_SEM)
        print (ten_LPS3_mean, "+/-", ten_LPS3_SEM)
        print (thirty_LPS3_mean, "+/-", thirty_LPS3_SEM)
        print (hundred_LPS3_mean, "+/-", hundred_LPS3_SEM)
        print (threehundred_LPS3_mean, "+/-", threehundred_LPS3_SEM)
        print (Mt_LPS3_mean, "+/-", Mt_LPS3_SEM)
        print ("")
        print ("========================================================")
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = N, dayta = Data, LPS_factor = True)

        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        #control_Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean]
        #LPS_Data = [BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #control_Labels = ("BBS", "10", "30", "100", "300", "M. tub")
        #LPS_Labels = ("LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(control_Labels))
        #x_pos2 = np.arange(len(LPS_Labels))

        #control_error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM]
        #LPS_error = [BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]

        #data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]
        
        #labels = ("BBS", "10", "30", "100", "300", "M. tub", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(labels))

        #plt.bar (x_pos, data, yerr = error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #plt.xticks (x_pos, labels)
        #plt.show()
        
        #fig, ax = plt.subplots()
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, ecolor = "black", capsize = 10)
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xlabel ("Treatment")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #plt.show()

        #fig2, ax2 = plt.subplots()
        #ax2.bar (x_pos, LPS_Data, yerr = LPS_error, align = "center", alpha = 0.75, ecolor = "black", capsize = 10)
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos)
        #ax2.set_xticklabels (LPS_Labels)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #plt.show()

        #fig, (ax, ax2) = plt.subplots (nrows = 2)
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax2.bar (x_pos2, LPS_Data, yerr = LPS_error, align = "center", alpha = 1, color = "red", ecolor = "black", capsize = 10)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos2)
        #ax2.set_xticklabels (LPS_Labels)
        #plt.show ()

        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        
        Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, 0, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]
        Labels = ("BBS", "10", "30", "100", "300", "TB", "", "BBS", "10", "30", "100", "300", "TB")
        x_pos = np.arange(len(Labels))
        Error = [(0,0,0,0,0,0,0,0,0,0,0,0,0), (BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, 0, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM)]

        Data2 = [BBS_LPS2_mean, ten_LPS2_mean, thirty_LPS2_mean, hundred_LPS2_mean, threehundred_LPS2_mean, Mt_LPS2_mean, 0, BBS_LPS3_mean, ten_LPS3_mean, thirty_LPS3_mean, hundred_LPS3_mean, threehundred_LPS3_mean, Mt_LPS3_mean]
        Labels2 = ("BBS", "10", "30", "100", "300", "TB", "", "BBS", "10", "30", "100", "300", "TB")
        x_pos2 = np.arange(len(Labels2))
        Error2 = [(0,0,0,0,0,0,0,0,0,0,0,0,0), (BBS_LPS2_SEM,ten_LPS2_SEM,thirty_LPS2_SEM,hundred_LPS2_SEM,threehundred_LPS2_SEM,Mt_LPS2_SEM,0,BBS_LPS3_SEM,ten_LPS3_SEM,thirty_LPS3_SEM,hundred_LPS3_SEM,threehundred_LPS3_SEM,Mt_LPS3_SEM)]
        
        fig, (ax,ax2) = plt.subplots (nrows = 2)
        #barlist = ax.bar (x_pos, Data, yerr = Error, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#B7D7EE", "#FFB5B1", "#FFB5B1", "#FFB5B1", "#FFB5B1", "#FACFBB", "b", "#7DB3D1", "#CF7E81", "#CF7E81", "#CF7E81", "#CF7E81", "#e6a27d"])
        barlist = ax.bar (x_pos, Data, yerr = Error, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#ADD4ED", "#FDB5AF", "#FDB5AF", "#FDB5AF", "#FDB5AF", "#FACFBB", "b", "#73B0D0", "#E08A81", "#E08A81", "#E08A81", "#E08A81", "#e6a27d"])
        ax.set_title ("Mean Log2 Values +/- 1 SEM")
        ax.set_ylabel ("Relative Gene Expression")
        ax.set_xticks(x_pos)
        ax.set_xticklabels (Labels, rotation = 25)
##        barlist[0].set_color ("#FFB5B5")
##        barlist[1].set_hatch ("/")
##        barlist[2].set_hatch ("/")
##        barlist[3].set_hatch ("/")
##        barlist[4].set_hatch ("/")
##        barlist[5].set_hatch ("/")
##        barlist[7].set_hatch ("//")
##        barlist[8].set_hatch ("//")
##        barlist[9].set_hatch ("//")
##        barlist[10].set_hatch ("//")
##        barlist[11].set_hatch ("//")
##        barlist[12].set_hatch ("//")
        
##        barlist[0].set_color ("#70C2FF")
##        barlist[1].set_color ("#FF95C4")
##        barlist[2].set_color ("#FF70AF")
##        barlist[3].set_color ("#FF4998")
##        barlist[4].set_color ("#FF006F")
##        barlist[5].set_color ("#000000")
##        barlist[7].set_color ("#70C2FF")
##        barlist[8].set_color ("#FF95C4")
##        barlist[9].set_color ("#FF70AF")
##        barlist[10].set_color ("#FF4998")
##        barlist[11].set_color ("#FF006F")
##        barlist[12].set_color ("#000000")

        #barlist2 = ax2.bar(x_pos2, Data2, yerr = Error2, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#4E8FB3", "#B5495B", "#B5495B", "#B5495B", "#B5495B", "#e8814c", "b", "#006B96", "#9B002E", "#9B002E", "#9B002E", "#9B002E", "#dc5c00"])
        barlist2 = ax2.bar(x_pos2, Data2, yerr = Error2, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#3A8BB3", "#C45E54", "#C45E54", "#C45E54", "#C45E54", "#e8814c", "b", "#006796", "#A73326", "#A73326", "#A73326", "#A73326", "#dc5c00"])
        ax2.set_ylabel ("Realtive Gene Expression")
        ax2.set_xticks(x_pos2)
        ax2.set_xticklabels(Labels2, rotation = 25)
        
##        plt.rcParams["font.weight"] = "bold"
##        plt.rcParams["axes.labelweight"] = "bold"
##        
##        Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, 0, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]
##        Labels = ("BBS", "10", "30", "100", "300", "TB", "", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS TB")
##        x_pos = np.arange(len(Labels))
##        Error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, 0, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]
##
##        Data2 = [BBS_LPS2_mean, ten_LPS2_mean, thirty_LPS2_mean, hundred_LPS2_mean, threehundred_LPS2_mean, Mt_LPS2_mean, 0, BBS_LPS3_mean, ten_LPS3_mean, thirty_LPS3_mean, hundred_LPS3_mean, threehundred_LPS3_mean, Mt_LPS3_mean]
##        Labels2 = ("LPS_2 BBS","LPS_2 10","LPS_2 30","LPS_2 100","LPS_2 300","LPS_2 TB", "", "LPS_3 BBS","LPS_3 10","LPS_3 30","LPS_3 100","LPS_3 300","LPS_3 TB")
##        x_pos2 = np.arange(len(Labels2))
##        Error2 = [BBS_LPS2_SEM,ten_LPS2_SEM,thirty_LPS2_SEM,hundred_LPS2_SEM,threehundred_LPS2_SEM,Mt_LPS2_SEM,0,BBS_LPS3_SEM,ten_LPS3_SEM,thirty_LPS3_SEM,hundred_LPS3_SEM,threehundred_LPS3_SEM,Mt_LPS3_SEM]
##        
##        fig, (ax,ax2) = plt.subplots (nrows = 2)
##        barlist = ax.bar (x_pos, Data, yerr = Error, align = "center", alpha = 1, capsize = 5)
##        ax.set_title ("Mean Log2 Values +/- 1 SEM")
##        ax.set_ylabel ("Relative Gene Expression")
##        ax.set_xticks(x_pos)
##        ax.set_xticklabels (Labels, rotation = 25)
##        barlist[0].set_color ("blue")
##        barlist[1].set_color ("blue")
##        barlist[2].set_color ("blue")
##        barlist[3].set_color ("blue")
##        barlist[4].set_color ("blue")
##        barlist[5].set_color ("blue")
##        barlist[7].set_color ("red")
##        barlist[8].set_color ("red")
##        barlist[9].set_color ("red")
##        barlist[10].set_color ("red")
##        barlist[11].set_color ("red")
##        barlist[12].set_color ("red")
##
##        barlist2 = ax2.bar(x_pos2, Data2, yerr = Error2, align = "center", alpha = 1, capsize = 5)
##        ax2.set_ylabel ("Realtive Gene Expression")
##        ax2.set_xticks(x_pos2)
##        ax2.set_xticklabels(Labels2, rotation = 25)
##        barlist2[0].set_color ("green")
##        barlist2[1].set_color ("green")
##        barlist2[2].set_color ("green")
##        barlist2[3].set_color ("green")
##        barlist2[4].set_color ("green")
##        barlist2[5].set_color ("green")
##        barlist2[7].set_color ("yellow")
##        barlist2[8].set_color ("yellow")
##        barlist2[9].set_color ("yellow")
##        barlist2[10].set_color ("yellow")
##        barlist2[11].set_color ("yellow")
##        barlist2[12].set_color ("yellow")

        #plt.xticks (x_pos, Labels, rotation = 45)
        #plt.ylabel ("Relative Gene Expression (Log2 Fold Change)")
        
        plt.show()
        
        #plt.bar (x_pos, control_Data, align = "center", alpha = 1)
        #plt.xticks (x_pos, control_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("Control Graph")
        #plt.show()
        
        #plt.bar (x_pos_2, LPS_Data, align = "center", alpha = 0.75)
        #plt.xticks (x_pos_2, LPS_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("LPS Graph")
        #plt.show ()
     
    return;

def PCR_grapher_2 (N, data, LPS_factor):
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import statistics

    csv = pandas.read_csv(data, header = None, na_filter = False)

    if LPS_factor == False:
        List_50 = []
        count = 0
        while count < (float(N) * 12):
            List_50.append (csv.iloc[count, 0])
            
            count += 1
            
        List = []
        for a in List_50:
            try:
                y = int (float(a))
                List.append (float(a))
            except:
                List.append (a)

        Data = []
        for a in List:
            if type(a) == float:
                Data.append(a)
            else:
                Data.append("NaN")

    if LPS_factor == True:
        List_50 = []
        count = 0
        while count < ((float(N) * 2) * 12):
            List_50.append (csv.iloc[count, 0])
            
            count += 1
            
        List = []
        for a in List_50:
            try:
                y = int (float(a))
                List.append (float(a))
            except:
                List.append (a)

        Data = []
        for a in List:
            if type(a) == float:
                Data.append(a)
            else:
                Data.append("NaN")
    
    if LPS_factor == False:
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        count = 0
        while count < (len(Data)):
            BBS.append (Data[count])
            ten.append (Data[count + 1])
            thirty.append (Data[count + 2])
            hundred.append (Data[count + 3])
            threehundred.append (Data[count + 4])
            Mt.append (Data[count + 5])

            BBS_LPS.append (Data[count + 6])
            ten_LPS.append (Data[count + 7])
            thirty_LPS.append (Data[count + 8])
            hundred_LPS.append (Data[count + 9])
            threehundred_LPS.append (Data[count + 10])
            Mt_LPS.append (Data[count + 11])

            count += 12
            
        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)
        BBS_mean = (sum (New_BBS) / len (New_BBS))
        
        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)
        ten_mean = (sum (New_ten) / len (New_ten))

        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)
        thirty_mean = (sum (New_thirty) / len (New_thirty))

        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)
        hundred_mean = (sum (New_hundred) / len (New_hundred))

        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)
        threehundred_mean = (sum (New_threehundred) / len (New_threehundred))

        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)
        Mt_mean = (sum (New_Mt) / len (New_Mt))

        #################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)
        BBS_LPS_mean = (sum (New_BBS_LPS) / len (New_BBS_LPS))

        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)
        ten_LPS_mean = (sum (New_ten_LPS) / len (New_ten_LPS))

        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)
        thirty_LPS_mean = (sum (New_thirty_LPS) / len (New_thirty_LPS))

        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)
        hundred_LPS_mean = (sum (New_hundred_LPS) / len (New_hundred_LPS))

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)
        threehundred_LPS_mean = (sum (New_threehundred_LPS) / len (New_threehundred_LPS))

        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)
        Mt_LPS_mean = (sum (New_Mt_LPS) / len (New_Mt_LPS))

        ###################################

        BBS_SEM = (statistics.stdev (New_BBS) / math.sqrt(len (New_BBS)))
        ten_SEM = (statistics.stdev (New_ten) / math.sqrt(len (New_ten)))
        thirty_SEM = (statistics.stdev (New_thirty) / math.sqrt(len (New_thirty)))
        hundred_SEM = (statistics.stdev (New_hundred) / math.sqrt(len (New_hundred)))
        threehundred_SEM = (statistics.stdev (New_threehundred) / math.sqrt(len (New_threehundred)))
        Mt_SEM = (statistics.stdev (New_Mt) / math.sqrt(len (New_Mt)))

        BBS_LPS_SEM = (statistics.stdev (New_BBS_LPS) / math.sqrt(len (New_BBS_LPS)))
        ten_LPS_SEM = (statistics.stdev (New_ten_LPS) / math.sqrt(len (New_ten_LPS)))
        thirty_LPS_SEM = (statistics.stdev (New_thirty_LPS) / math.sqrt(len (New_thirty_LPS)))
        hundred_LPS_SEM = (statistics.stdev (New_hundred_LPS) / math.sqrt(len (New_hundred_LPS)))
        threehundred_LPS_SEM = (statistics.stdev (New_threehundred_LPS) / math.sqrt(len (New_threehundred_LPS)))
        Mt_LPS_SEM = (statistics.stdev (New_Mt_LPS) / math.sqrt(len (New_Mt_LPS)))
        
        print ("========================================================")
        print ("")
        print ("Your Average Log2 Values +/- SEM Are: ")
        print ("")
        print ("Control:")
        print (BBS_mean, "+/-", BBS_SEM)
        print (ten_mean, "+/-", ten_SEM)
        print (thirty_mean, "+/-", thirty_SEM)
        print (hundred_mean, "+/-", hundred_SEM)
        print (threehundred_mean, "+/-", threehundred_SEM)
        print (Mt_mean, "+/-", Mt_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("LPS:")
        print (BBS_LPS_mean, "+/-", BBS_LPS_SEM)
        print (ten_LPS_mean, "+/-", ten_LPS_SEM)
        print (thirty_LPS_mean, "+/-", thirty_LPS_SEM)
        print (hundred_LPS_mean, "+/-", hundred_LPS_SEM)
        print (threehundred_LPS_mean, "+/-", threehundred_LPS_SEM)
        print (Mt_LPS_mean, "+/-", Mt_LPS_SEM)
        print ("")
        print ("========================================================")
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = N, dayta = Data, LPS_factor = False)

        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        #control_Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean]
        #LPS_Data = [BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #control_Labels = ("BBS", "10", "30", "100", "300", "M. tub")
        #LPS_Labels = ("LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(control_Labels))
        #x_pos2 = np.arange(len(LPS_Labels))

        #control_error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM]
        #LPS_error = [BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]


        #data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]
        
        #labels = ("BBS", "10", "30", "100", "300", "M. tub", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(labels))

        #plt.bar (x_pos, data, yerr = error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #plt.xticks (x_pos, labels)
        #plt.show()
        
        #fig, ax = plt.subplots()
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, ecolor = "black", capsize = 10)
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xlabel ("Treatment")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #plt.show()

        #fig2, ax2 = plt.subplots()
        #ax2.bar (x_pos, LPS_Data, yerr = LPS_error, align = "center", alpha = 0.75, ecolor = "black", capsize = 10)
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos)
        #ax2.set_xticklabels (LPS_Labels)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #plt.show()

        #fig, (ax, ax2) = plt.subplots (nrows = 2)
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        

        #ax2.bar (x_pos2, LPS_Data, yerr = LPS_error, align = "center", alpha = 1, color = "red", ecolor = "black", capsize = 10)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos2)
        #ax2.set_xticklabels (LPS_Labels)
        #plt.show ()
        
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        
        Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, 0, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        Labels = ("BBS", "10", "30", "100", "300", "TB", "", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS TB")
        
        x_pos = np.arange(len(Labels))
       
        Error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, 0, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]

        fig, ax = plt.subplots (figsize = (9.25, 7.25))
        barlist = plt.bar (x_pos, Data, align = "center", yerr = Error, ecolor = "black", capsize = 5, alpha = 1)
        
        barlist [0].set_color ("blue")
        barlist [1].set_color ("blue")
        barlist [2].set_color ("blue")
        barlist [3].set_color ("blue")
        barlist [4].set_color ("blue")
        barlist [5].set_color ("blue")

        barlist [7].set_color ("red")
        barlist [8].set_color ("red")
        barlist [9].set_color ("red")
        barlist [10].set_color ("red")
        barlist [11].set_color ("red")
        barlist [12].set_color ("red")

        plt.xticks (x_pos, Labels, rotation = 45)
        plt.ylabel ("Relative Gene Expression (Log2 Fold Change)")
        
        plt.show()
        
        #plt.bar (x_pos, control_Data, align = "center", alpha = 1)
        #plt.xticks (x_pos, control_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("Control Graph")
        #plt.show()
        
        #plt.bar (x_pos_2, LPS_Data, align = "center", alpha = 0.75)
        #plt.xticks (x_pos_2, LPS_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("LPS Graph")
        #plt.show ()

    if LPS_factor == True:
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        BBS_LPS2 = []
        ten_LPS2 = []
        thirty_LPS2 = []
        hundred_LPS2 = []
        threehundred_LPS2 = []
        Mt_LPS2 = []

        BBS_LPS3 = []
        ten_LPS3 = []
        thirty_LPS3 = []
        hundred_LPS3 = []
        threehundred_LPS3 = []
        Mt_LPS3 = []

        count = 0
        while count < (len(Data)):
            BBS.append (Data[count])
            ten.append (Data[count + 1])
            thirty.append (Data[count + 2])
            hundred.append (Data[count + 3])
            threehundred.append (Data[count + 4])
            Mt.append (Data[count + 5])

            BBS_LPS.append (Data[count + 6])
            ten_LPS.append (Data[count + 7])
            thirty_LPS.append (Data[count + 8])
            hundred_LPS.append (Data[count + 9])
            threehundred_LPS.append (Data[count + 10])
            Mt_LPS.append (Data[count + 11])

            BBS_LPS2.append (Data[count + 12])
            ten_LPS2.append (Data[count + 13])
            thirty_LPS2.append (Data[count + 14])
            hundred_LPS2.append (Data[count + 15])
            threehundred_LPS2.append (Data[count + 16])
            Mt_LPS2.append (Data[count + 17])

            BBS_LPS3.append (Data[count + 18])
            ten_LPS3.append (Data[count + 19])
            thirty_LPS3.append (Data[count + 20])
            hundred_LPS3.append (Data[count + 21])
            threehundred_LPS3.append (Data[count + 22])
            Mt_LPS3.append (Data[count + 23])
            
            count += 24
            
        New_BBS = []
        for a in BBS:
            if type(a) == str:
                continue
            New_BBS.append (a)
        try:
            BBS_mean = (sum (New_BBS) / len (New_BBS))
        except:
            BBS_mean = 0
        
        New_ten = []
        for b in ten:
            if type(b) == str:
                continue
            New_ten.append (b)
        try:
            ten_mean = (sum (New_ten) / len (New_ten))
        except:
            ten_mean = 0

        New_thirty = []
        for c in thirty:
            if type(c) == str:
                continue
            New_thirty.append (c)
        try:
            thirty_mean = (sum (New_thirty) / len (New_thirty))
        except:
            thirty_mean = 0

        New_hundred = []
        for d in hundred:
            if type(d) == str:
                continue
            New_hundred.append (d)
        try:
            hundred_mean = (sum (New_hundred) / len (New_hundred))
        except:
            hundred_mean = 0

        New_threehundred = []
        for e in threehundred:
            if type(e) == str:
                continue
            New_threehundred.append (e)
        try:
            threehundred_mean = (sum (New_threehundred) / len (New_threehundred))
        except:
            threehundred_mean = 0

        New_Mt = []
        for f in Mt:
            if type(f) == str:
                continue
            New_Mt.append (f)
        try:
            Mt_mean = (sum (New_Mt) / len (New_Mt))
        except:
            Mt_mean = 0

        #################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            if type(g) == str:
                continue
            New_BBS_LPS.append (g)
        BBS_LPS_mean = (sum (New_BBS_LPS) / len (New_BBS_LPS))

        New_ten_LPS = []
        for h in ten_LPS:
            if type(h) == str:
                continue
            New_ten_LPS.append (h)
        ten_LPS_mean = (sum (New_ten_LPS) / len (New_ten_LPS))

        New_thirty_LPS = []
        for i in thirty_LPS:
            if type(i) == str:
                continue
            New_thirty_LPS.append (i)
        thirty_LPS_mean = (sum (New_thirty_LPS) / len (New_thirty_LPS))

        New_hundred_LPS = []
        for j in hundred_LPS:
            if type(j) == str:
                continue
            New_hundred_LPS.append (j)
        hundred_LPS_mean = (sum (New_hundred_LPS) / len (New_hundred_LPS))

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            if type(k) == str:
                continue
            New_threehundred_LPS.append (k)
        threehundred_LPS_mean = (sum (New_threehundred_LPS) / len (New_threehundred_LPS))

        New_Mt_LPS = []
        for l in Mt_LPS:
            if type(l) == str:
                continue
            New_Mt_LPS.append (l)
        Mt_LPS_mean = (sum (New_Mt_LPS) / len (New_Mt_LPS))

        ###################################

        New_BBS_LPS2 = []
        for m in BBS_LPS2:
            if type(m) == str:
                continue
            New_BBS_LPS2.append (m)
        BBS_LPS2_mean = (sum (New_BBS_LPS2) / len (New_BBS_LPS2))

        New_ten_LPS2 = []
        for n in ten_LPS2:
            if type(n) == str:
                continue
            New_ten_LPS2.append (n)
        ten_LPS2_mean = (sum (New_ten_LPS2) / len (New_ten_LPS2))

        New_thirty_LPS2 = []
        for o in thirty_LPS2:
            if type(o) == str:
                continue
            New_thirty_LPS2.append (o)
        thirty_LPS2_mean = (sum (New_thirty_LPS2) / len (New_thirty_LPS2))

        New_hundred_LPS2 = []
        for p in hundred_LPS2:
            if type(p) == str:
                continue
            New_hundred_LPS2.append (p)
        hundred_LPS2_mean = (sum (New_hundred_LPS2) / len (New_hundred_LPS2))

        New_threehundred_LPS2 = []
        for q in threehundred_LPS2:
            if type(q) == str:
                continue
            New_threehundred_LPS2.append (q)
        threehundred_LPS2_mean = (sum (New_threehundred_LPS2) / len (New_threehundred_LPS2))

        New_Mt_LPS2 = []
        for r in Mt_LPS2:
            if type(r) == str:
                continue
            New_Mt_LPS2.append (r)
        Mt_LPS2_mean = (sum (New_Mt_LPS2) / len (New_Mt_LPS2))

        ###################################
        
        New_BBS_LPS3 = []
        for s in BBS_LPS3:
            if type(s) == str:
                continue
            New_BBS_LPS3.append (s)
        BBS_LPS3_mean = (sum (New_BBS_LPS3) / len (New_BBS_LPS3))

        New_ten_LPS3 = []
        for t in ten_LPS3:
            if type(t) == str:
                continue
            New_ten_LPS3.append (t)
        ten_LPS3_mean = (sum (New_ten_LPS3) / len (New_ten_LPS3))

        New_thirty_LPS3 = []
        for u in thirty_LPS3:
            if type(u) == str:
                continue
            New_thirty_LPS3.append (u)
        thirty_LPS3_mean = (sum (New_thirty_LPS3) / len (New_thirty_LPS3))

        New_hundred_LPS3 = []
        for v in hundred_LPS3:
            if type(v) == str:
                continue
            New_hundred_LPS3.append (v)
        hundred_LPS3_mean = (sum (New_hundred_LPS3) / len (New_hundred_LPS3))

        New_threehundred_LPS3 = []
        for w in threehundred_LPS3:
            if type(w) == str:
                continue
            New_threehundred_LPS3.append (w)
        threehundred_LPS3_mean = (sum (New_threehundred_LPS3) / len (New_threehundred_LPS3))

        New_Mt_LPS3 = []
        for x in Mt_LPS3:
            if type(x) == str:
                continue
            New_Mt_LPS3.append (x)
        Mt_LPS3_mean = (sum (New_Mt_LPS3) / len (New_Mt_LPS3))

        ###################################

        try:
            BBS_SEM = (statistics.stdev (New_BBS) / math.sqrt(len (New_BBS)))
        except:
            BBS_SEM = 0

        try:
            ten_SEM = (statistics.stdev (New_ten) / math.sqrt(len (New_ten)))
        except:
            ten_SEM = 0

        try:
            thirty_SEM = (statistics.stdev (New_thirty) / math.sqrt(len (New_thirty)))
        except:
            thirty_SEM = 0

        try:
            hundred_SEM = (statistics.stdev (New_hundred) / math.sqrt(len (New_hundred)))
        except:
            hundred_SEM = 0

        try:
            threehundred_SEM = (statistics.stdev (New_threehundred) / math.sqrt(len (New_threehundred)))
        except:
            threehundred_SEM = 0

        try:
            Mt_SEM = (statistics.stdev (New_Mt) / math.sqrt(len (New_Mt)))
        except:
            Mt_SEM = 0
        
        try:
            BBS_LPS_SEM = (statistics.stdev (New_BBS_LPS) / math.sqrt(len (New_BBS_LPS)))
        except:
            BBS_LPS_SEM = 0
            
        try:
            ten_LPS_SEM = (statistics.stdev (New_ten_LPS) / math.sqrt(len (New_ten_LPS)))
        except:
            ten_LPS_SEM = 0

        try:
            thirty_LPS_SEM = (statistics.stdev (New_thirty_LPS) / math.sqrt(len (New_thirty_LPS)))
        except:
            thirty_LPS_SEM = 0

        try:
            hundred_LPS_SEM = (statistics.stdev (New_hundred_LPS) / math.sqrt(len (New_hundred_LPS)))
        except:
            hundred_LPS_SEM = 0

        try:
            threehundred_LPS_SEM = (statistics.stdev (New_threehundred_LPS) / math.sqrt(len (New_threehundred_LPS)))
        except:
            threehundred_LPS_SEM = 0

        try:
            Mt_LPS_SEM = (statistics.stdev (New_Mt_LPS) / math.sqrt(len (New_Mt_LPS)))
        except:
            Mt_LPS_SEM = 0
        
        try:
            BBS_LPS2_SEM = (statistics.stdev (New_BBS_LPS2) / math.sqrt(len(New_BBS_LPS2)))
        except:
            BBS_LPS2_SEM = 0

        try:
            ten_LPS2_SEM = (statistics.stdev (New_ten_LPS2) / math.sqrt(len(New_ten_LPS2)))
        except:
            ten_LPS2_SEM = 0

        try:
            thirty_LPS2_SEM = (statistics.stdev (New_thirty_LPS2) / math.sqrt(len(New_thirty_LPS2)))
        except:
            thirty_LPS2_SEM = 0

        try:
            hundred_LPS2_SEM = (statistics.stdev (New_hundred_LPS2) / math.sqrt(len(New_hundred_LPS2)))
        except:
            hundred_LPS2_SEM = 0

        try:
            threehundred_LPS2_SEM = (statistics.stdev (New_threehundred_LPS2) / math.sqrt(len(New_threehundred_LPS2)))
        except:
            threehundred_LPS2_SEM = 0

        try:
            Mt_LPS2_SEM = (statistics.stdev (New_Mt_LPS2) / math.sqrt(len(New_Mt_LPS2)))
        except:
            Mt_LPS2_SEM = 0

        try:
            BBS_LPS3_SEM = (statistics.stdev (New_BBS_LPS3) / math.sqrt(len(New_BBS_LPS3)))
        except:
            BBS_LPS3_SEM = 0

        try:
            ten_LPS3_SEM = (statistics.stdev (New_ten_LPS3) / math.sqrt(len(New_ten_LPS3)))
        except:
            ten_LPS3_SEM = 0

        try:
            thirty_LPS3_SEM = (statistics.stdev (New_thirty_LPS3) / math.sqrt(len(New_thirty_LPS3)))
        except:
            thirty_LPS3_SEM = 0

        try:
            hundred_LPS3_SEM = (statistics.stdev (New_hundred_LPS3) / math.sqrt(len(New_hundred_LPS3)))
        except:
            hundred_LPS3_SEM = 0

        try:
            threehundred_LPS3_SEM = (statistics.stdev (New_threehundred_LPS3) / math.sqrt(len(New_threehundred_LPS3)))
        except:
            threehundred_LPS3_SEM = 0

        try:
            Mt_LPS3_SEM = (statistics.stdev (New_Mt_LPS3) / math.sqrt(len(New_Mt_LPS3)))
        except:
            Mt_LPS3_SEM = 0
        
        
        print ("========================================================")
        print ("")
        print ("Your Average Log2 Values +/- SEM Are: ")
        print ("")
        print ("Control:")
        print (BBS_mean, "+/-", BBS_SEM)
        print (ten_mean, "+/-", ten_SEM)
        print (thirty_mean, "+/-", thirty_SEM)
        print (hundred_mean, "+/-", hundred_SEM)
        print (threehundred_mean, "+/-", threehundred_SEM)
        print (Mt_mean, "+/-", Mt_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("1st LPS Dose:")
        print (BBS_LPS_mean, "+/-", BBS_LPS_SEM)
        print (ten_LPS_mean, "+/-", ten_LPS_SEM)
        print (thirty_LPS_mean, "+/-", thirty_LPS_SEM)
        print (hundred_LPS_mean, "+/-", hundred_LPS_SEM)
        print (threehundred_LPS_mean, "+/-", threehundred_LPS_SEM)
        print (Mt_LPS_mean, "+/-", Mt_LPS_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("2nd LPS Dose:")
        print (BBS_LPS2_mean, "+/-", BBS_LPS2_SEM)
        print (ten_LPS2_mean, "+/-", ten_LPS2_SEM)
        print (thirty_LPS2_mean, "+/-", thirty_LPS2_SEM)
        print (hundred_LPS2_mean, "+/-", hundred_LPS2_SEM)
        print (threehundred_LPS2_mean, "+/-", threehundred_LPS2_SEM)
        print (Mt_LPS2_mean, "+/-", Mt_LPS2_SEM)
        print ("")
        print ("////////////////////////////////////////////////////////")
        print ("")
        print ("3rd LPS Dose:")
        print (BBS_LPS3_mean, "+/-", BBS_LPS3_SEM)
        print (ten_LPS3_mean, "+/-", ten_LPS3_SEM)
        print (thirty_LPS3_mean, "+/-", thirty_LPS3_SEM)
        print (hundred_LPS3_mean, "+/-", hundred_LPS3_SEM)
        print (threehundred_LPS3_mean, "+/-", threehundred_LPS3_SEM)
        print (Mt_LPS3_mean, "+/-", Mt_LPS3_SEM)
        print ("")
        print ("========================================================")
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = N, dayta = Data, LPS_factor = True)

        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

        #control_Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean]
        #LPS_Data = [BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #control_Labels = ("BBS", "10", "30", "100", "300", "M. tub")
        #LPS_Labels = ("LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(control_Labels))
        #x_pos2 = np.arange(len(LPS_Labels))

        #control_error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM]
        #LPS_error = [BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]

        #data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]

        #error = [BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM]
        
        #labels = ("BBS", "10", "30", "100", "300", "M. tub", "LPS BBS", "LPS 10", "LPS 30", "LPS 100", "LPS 300", "LPS M. tub")

        #x_pos = np.arange(len(labels))

        #plt.bar (x_pos, data, yerr = error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #plt.xticks (x_pos, labels)
        #plt.show()
        
        #fig, ax = plt.subplots()
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, ecolor = "black", capsize = 10)
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xlabel ("Treatment")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #plt.show()

        #fig2, ax2 = plt.subplots()
        #ax2.bar (x_pos, LPS_Data, yerr = LPS_error, align = "center", alpha = 0.75, ecolor = "black", capsize = 10)
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos)
        #ax2.set_xticklabels (LPS_Labels)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #plt.show()

        #fig, (ax, ax2) = plt.subplots (nrows = 2)
        #ax.bar (x_pos, control_Data, yerr = control_error, align = "center", alpha = 1, color = "blue", ecolor = "black", capsize = 10)
        #ax.set_title ("Control Log2 Graph +/- SEM")
        #ax.set_ylabel ("Relative Gene Expression")
        #ax.set_xticks (x_pos)
        #ax.set_xticklabels (control_Labels)
        #ax2.bar (x_pos2, LPS_Data, yerr = LPS_error, align = "center", alpha = 1, color = "red", ecolor = "black", capsize = 10)
        #ax2.set_title ("LPS Log2 Graph +/- SEM")
        #ax2.set_ylabel ("Relative Gene Expression")
        #ax2.set_xlabel ("Treatment")
        #ax2.set_xticks (x_pos2)
        #ax2.set_xticklabels (LPS_Labels)
        #plt.show ()
        
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        
        Data = [BBS_mean, ten_mean, thirty_mean, hundred_mean, threehundred_mean, Mt_mean, 0, BBS_LPS_mean, ten_LPS_mean, thirty_LPS_mean, hundred_LPS_mean, threehundred_LPS_mean, Mt_LPS_mean]
        Labels = ("BBS", "10", "30", "100", "300", "TB", "", "BBS", "10", "30", "100", "300", "TB")
        x_pos = np.arange(len(Labels))
        Error = [(0,0,0,0,0,0,0,0,0,0,0,0,0), (BBS_SEM, ten_SEM, thirty_SEM, hundred_SEM, threehundred_SEM, Mt_SEM, 0, BBS_LPS_SEM, ten_LPS_SEM, thirty_LPS_SEM, hundred_LPS_SEM, threehundred_LPS_SEM, Mt_LPS_SEM)]

        Data2 = [BBS_LPS2_mean, ten_LPS2_mean, thirty_LPS2_mean, hundred_LPS2_mean, threehundred_LPS2_mean, Mt_LPS2_mean, 0, BBS_LPS3_mean, ten_LPS3_mean, thirty_LPS3_mean, hundred_LPS3_mean, threehundred_LPS3_mean, Mt_LPS3_mean]
        Labels2 = ("BBS", "10", "30", "100", "300", "TB", "", "BBS", "10", "30", "100", "300", "TB")
        x_pos2 = np.arange(len(Labels2))
        Error2 = [(0,0,0,0,0,0,0,0,0,0,0,0,0), (BBS_LPS2_SEM,ten_LPS2_SEM,thirty_LPS2_SEM,hundred_LPS2_SEM,threehundred_LPS2_SEM,Mt_LPS2_SEM,0,BBS_LPS3_SEM,ten_LPS3_SEM,thirty_LPS3_SEM,hundred_LPS3_SEM,threehundred_LPS3_SEM,Mt_LPS3_SEM)]
        
        fig, (ax,ax2) = plt.subplots (nrows = 2)
        #barlist = ax.bar (x_pos, Data, yerr = Error, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#B7D7EE", "#FFB5B1", "#FFB5B1", "#FFB5B1", "#FFB5B1", "#FACFBB", "b", "#7DB3D1", "#CF7E81", "#CF7E81", "#CF7E81", "#CF7E81", "#e6a27d"])
        barlist = ax.bar (x_pos, Data, yerr = Error, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#ADD4ED", "#FDB5AF", "#FDB5AF", "#FDB5AF", "#FDB5AF", "#FACFBB", "b", "#73B0D0", "#E08A81", "#E08A81", "#E08A81", "#E08A81", "#e6a27d"])
        ax.set_title ("Mean Log2 Values +/- 1 SEM")
        ax.set_ylabel ("Relative Gene Expression")
        ax.set_xticks(x_pos)
        ax.set_xticklabels (Labels, rotation = 25)
##        barlist[0].set_color ("#FFB5B5")
##        barlist[1].set_hatch ("/")
##        barlist[2].set_hatch ("/")
##        barlist[3].set_hatch ("/")
##        barlist[4].set_hatch ("/")
##        barlist[5].set_hatch ("/")
##        barlist[7].set_hatch ("//")
##        barlist[8].set_hatch ("//")
##        barlist[9].set_hatch ("//")
##        barlist[10].set_hatch ("//")
##        barlist[11].set_hatch ("//")
##        barlist[12].set_hatch ("//")
        
##        barlist[0].set_color ("#70C2FF")
##        barlist[1].set_color ("#FF95C4")
##        barlist[2].set_color ("#FF70AF")
##        barlist[3].set_color ("#FF4998")
##        barlist[4].set_color ("#FF006F")
##        barlist[5].set_color ("#000000")
##        barlist[7].set_color ("#70C2FF")
##        barlist[8].set_color ("#FF95C4")
##        barlist[9].set_color ("#FF70AF")
##        barlist[10].set_color ("#FF4998")
##        barlist[11].set_color ("#FF006F")
##        barlist[12].set_color ("#000000")

        #barlist2 = ax2.bar(x_pos2, Data2, yerr = Error2, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#4E8FB3", "#B5495B", "#B5495B", "#B5495B", "#B5495B", "#e8814c", "b", "#006B96", "#9B002E", "#9B002E", "#9B002E", "#9B002E", "#dc5c00"])
        barlist2 = ax2.bar(x_pos2, Data2, yerr = Error2, align = "center", alpha = 1, capsize = 5, edgecolor = "k", linewidth = 1.5, color = ["#3A8BB3", "#C45E54", "#C45E54", "#C45E54", "#C45E54", "#e8814c", "b", "#006796", "#A73326", "#A73326", "#A73326", "#A73326", "#dc5c00"])
        ax2.set_ylabel ("Realtive Gene Expression")
        ax2.set_xticks(x_pos2)
        ax2.set_xticklabels(Labels2, rotation = 25)
##        barlist2[0].set_hatch ("x")
##        barlist2[1].set_hatch ("x")
##        barlist2[2].set_hatch ("x")
##        barlist2[3].set_hatch ("x")
##        barlist2[4].set_hatch ("x")
##        barlist2[5].set_hatch ("x")
##        barlist2[7].set_hatch (".")
##        barlist2[8].set_hatch (".")
##        barlist2[9].set_hatch (".")
##        barlist2[10].set_hatch (".")
##        barlist2[11].set_hatch (".")
##        barlist2[12].set_hatch (".")

##        barlist2[0].set_color ("#70C2FF")
##        barlist2[1].set_color ("#FF95C4")
##        barlist2[2].set_color ("#FF70AF")
##        barlist2[3].set_color ("#FF4998")
##        barlist2[4].set_color ("#FF006F")
##        barlist2[5].set_color ("#000000")
##        barlist2[7].set_color ("#70C2FF")
##        barlist2[8].set_color ("#FF95C4")
##        barlist2[9].set_color ("#FF70AF")
##        barlist2[10].set_color ("#FF4998")
##        barlist2[11].set_color ("#FF006F")
##        barlist2[12].set_color ("#000000")

        #plt.xticks (x_pos, Labels, rotation = 45)
        #plt.ylabel ("Relative Gene Expression (Log2 Fold Change)")
        
        plt.show()
        
        #plt.bar (x_pos, control_Data, align = "center", alpha = 1)
        #plt.xticks (x_pos, control_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("Control Graph")
        #plt.show()
        
        #plt.bar (x_pos_2, LPS_Data, align = "center", alpha = 0.75)
        #plt.xticks (x_pos_2, LPS_Labels)
        #plt.ylabel ("Relative Gene Expression")
        #plt.title ("LPS Graph")
        #plt.show ()
     
    return;

def PCR_analyzer (n, data, LPS_factor):

    print ("========================================================")
    if LPS_factor == True:
        n *= 2
    
    csv = pandas.read_csv(data, header = None, na_filter = False)

    List_50 = []
    List_51 = []
    count = 0
    while count < ((float(n) * 2) * 12):
        List_50.append (csv.iloc[count, 0])
        List_51.append (csv.iloc[count, 1])
        
        count += 1

    List = []
    for a in List_50:
        try:
            y = int (float(a))
            List.append (float(a))
        except:
            List.append (a)
            
    List_2 = []
    for b in List_51:
        try:
            z = int (float(b))
            List_2.append (float(b))
        except:
            List_2.append (b)

    List_3 = []
    count_2 = 0
    count_3 = 1
    a = 0
    b = 0
    while count_2 < (len(List) - 12):
        
        if (count_3 % 13) == 0:
            count_2 += 12
            count_3 = 1
            a = 1
            b = 1
        elif b == 1:
            count_2 += 1
            count_3 += 1
            b = 0
        
        Replicate_1 = List[count_2]
        Replicate_2 = List[12 + count_2]

        if type(List[count_2]) == str: 
            Replicate_Average = List[12 + count_2]
        elif type(List[12 + count_2]) == str:
            Replicate_Average = List[count_2]
        else:
            Replicate_Average = (Replicate_1 + Replicate_2) / 2
            
        List_3.append (Replicate_Average)

        if a == 1:
            a = 0
        else:
            count_2 += 1
            count_3 += 1

    # List 3 is the cytokine average Cq values
    print ("Cytokine average Cq values:")
    for a in List_3:
        print(a)

    
    List_7 = []
    count_8 = 0
    count_9 = 1
    c = 0
    d = 0
    while count_8 < (len(List_2) - 12):
        
        if (count_9 % 13) == 0:
            count_8 += 12
            count_9 = 1
            c = 1
            d = 1
        elif d == 1:
            count_8 += 1
            count_9 += 1
            d = 0
        
        Replicate_3 = List_2[count_8]
        Replicate_4 = List_2[12 + count_8]

        if type(List_2[count_8]) == str: 
            Replicate_Average_2 = List_2[12 + count_8]
        elif type(List_2[12 + count_8]) == str:
            Replicate_Average_2 = List_2[count_8] 
        else:
            Replicate_Average_2 = (Replicate_3 + Replicate_4) / 2
            
        List_7.append (Replicate_Average_2)

        if c == 1:
            c = 0
        else:
            count_8 += 1
            count_9 += 1

    List_4 = []
    count_4 = 0
    while count_4 < (12 * float(n)):
        C_avg = List_3[count_4]
        B_actin = List_7[count_4]

        if type(List_3[count_4]) == str:
            Normalized = "NaN"
        elif type(List_7[count_4]) == str:
            Normalized = "NaN"
        else:
            Normalized = C_avg - B_actin

        List_4.append (Normalized)

        count_4 += 1

    
    List_Floaty = []
    count_9 = 0
    while count_9 < len(List_4):
        if type (List_4[count_9]) == str:
            Floaty = -1
        else:
            Floaty = float(List_4[count_9])
        List_Floaty.append (Floaty)
        count_9 += 1


    Max_Normal = max (List_Floaty)


    List_5 = []
    count_5 = 0
    while count_5 < (12 * float(n)):
        if type (List_4[count_5]) == str:
            Delta_Delta = "NaN"
        else:
            Delta_Delta = List_4[count_5] - Max_Normal
            
        List_5.append (Delta_Delta)
        
        count_5 += 1



    List_6 = []
    count_6 = 0
    while count_6 < (12 * float(n)):
        if type (List_5[count_6]) == str:
            Log2 = "NaN"
        else:
            Log2 = 2**(List_5[count_6] * -1)
            
        List_6.append (Log2)
        
        count_6 += 1


    print ("")
    print ("Your Log2 Values Are:")
    print ("")

    count_7 = 0
    while count_7 < len(List_6):
        print (List_6[count_7])
        
        count_7 += 1

    if LPS_factor == True:
        n = n/2

    print ("")
    print ("========================================================")
    
    step2 = input ("Would You Like to Graph Your Log2 Values? <Y/N> ")

    if step2 == "Y" or step2 == "y" or step2 == "Yes" or step2 == "yes" or step2 == "YES":
        PCR_grapher (N = n, Data = List_6, LPS_factor = LPS_factor) 
    elif step2 == "N" or step2 == "n" or step2 == "No" or step2 == "no" or step2 == "NO":
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = n, dayta = List_6, LPS_factor = LPS_factor)
        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    return;

def PCR_stats (n, data, LPS_factor):

    flex = n
    
    import scipy.stats

    if LPS_factor == False:
        
        csv = pandas.read_csv (data, header = None, na_filter = False)

        dayta = []
        count = 0
        while count < (float(n) * 12):
            dayta.append (csv.iloc[count, 0])
            
            count += 1
      
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        count = 0
        while count < (len(dayta)):
            BBS.append (dayta[count])
            ten.append (dayta[count + 1])
            thirty.append (dayta[count + 2])
            hundred.append (dayta[count + 3])
            threehundred.append (dayta[count + 4])
            Mt.append (dayta[count + 5])

            BBS_LPS.append (dayta[count + 6])
            ten_LPS.append (dayta[count + 7])
            thirty_LPS.append (dayta[count + 8])
            hundred_LPS.append (dayta[count + 9])
            threehundred_LPS.append (dayta[count + 10])
            Mt_LPS.append (dayta[count + 11])

            count += 12

        New_BBS = []
        for a in BBS:
            try:
                z = int (float(a))
                New_BBS.append (float(a))
            except:
                continue

        New_ten = []
        for b in ten:
            try:
                z = int (float(b))
                New_ten.append (float(b))
            except:
                continue

        New_thirty = []
        for c in thirty:
            try:
                z = int (float(c))
                New_thirty.append (float(c))
            except:
                continue

        New_hundred = []
        for d in hundred:
            try:
                z = int (float(d))
                New_hundred.append (float(d))
            except:
                continue

        New_threehundred = []
        for e in threehundred:
            try:
                z = int (float(e))
                New_threehundred.append (float(e))
            except:
                continue

        New_Mt = []
        for f in Mt:
            try:
                z = int (float(f))
                New_Mt.append (float(f))
            except:
                continue

        ################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            try:
                z = int (float(g))
                New_BBS_LPS.append (float(g))
            except:
                continue

        New_ten_LPS = []
        for h in ten_LPS:
            try:
                z = int (float(h))
                New_ten_LPS.append (float(h))
            except:
                continue

        New_thirty_LPS = []
        for i in thirty_LPS:
            try:
                z = int (float(i))
                New_thirty_LPS.append (float(i))
            except:
                continue

        New_hundred_LPS = []
        for j in hundred_LPS:
            try:
                z = int (float(j))
                New_hundred_LPS.append (float(j))
            except:
                continue

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            try:
                z = int (float(k))
                New_threehundred_LPS.append (float(k))
            except:
                continue

        New_Mt_LPS = []
        for l in Mt_LPS:
            try:
                z = int (float(l))
                New_Mt_LPS.append (float(l))
            except:
                continue

        ################################

        control_results = scipy.stats.f_oneway (New_BBS, New_ten, New_thirty, New_hundred, New_threehundred)
        LPS_results = scipy.stats.f_oneway (New_BBS_LPS, New_ten_LPS, New_thirty_LPS, New_hundred_LPS, New_threehundred_LPS)

        control_t_results = scipy.stats.ttest_ind (New_BBS, New_Mt, equal_var = True)
        LPS_t_results = scipy.stats.ttest_ind (New_BBS_LPS, New_Mt_LPS, equal_var = True)

        ################################
        
        print ("========================================================")
        print ("")
        print ("Your One-Way ANOVA Values (BBS to 300) Are: ")
        print ("")
        print ("Control F value, P value:", control_results)
        print ("")
        print ("LPS F value, P value:", LPS_results)
        print ("")
        
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("")
        print ("Your T Test Values (BBS vs. M. tub) Are: ")
        print ("")
        print ("Control t-statistic, P value:", control_t_results)
        print ("")
        print ("LPS t-statistic, P value:", LPS_t_results)
        print ("")
        print ("========================================================")

        Next = input ("Run Post-Hoc Test? <Y/N> ")

        if Next == "Y" or Next == "y" or Next == "Yes" or Next == "yes" or Next == "YES":
            PCR_posthoc_test (N = flex, BBS = New_BBS, ten = New_ten, thirty = New_thirty, hundred = New_hundred,
                         threehundred = New_threehundred, Mt = New_Mt, BBS_LPS = New_BBS_LPS, ten_LPS = New_ten_LPS,
                         thirty_LPS = New_thirty_LPS, hundred_LPS = New_hundred_LPS, threehundred_LPS =
                         New_threehundred_LPS, Mt_LPS = New_Mt_LPS)

        elif Next == "N" or Next == "n" or Next == "No" or Next == "no" or Next == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 

                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    elif LPS_factor == True:
        from statsmodels.formula.api import ols
        from statsmodels.stats.anova import anova_lm
        import numpy
        import statistics
        from scipy import stats
        
        csv = pandas.read_csv (data, header = None, na_filter = False)

        Dayta,count = [],0
        while count < (float(n) * 24):
            Dayta.append (csv.iloc[count, 0])
            count += 1
                
        dayta,count = [],0
        while count < len(Dayta):
            try:
                z = int(float(Dayta[count]))
                dayta.append(float(Dayta[count]))
            except:
                dayta.append ("NaN")
            count += 1
        
        BBS = []
        ten = []
        thirty = []
        hundred = []
        threehundred = []
        Mt = []

        BBS_LPS = []
        ten_LPS = []
        thirty_LPS = []
        hundred_LPS = []
        threehundred_LPS = []
        Mt_LPS = []

        BBS_LPS2 = []
        ten_LPS2 = []
        thirty_LPS2 = []
        hundred_LPS2 = []
        threehundred_LPS2 = []
        Mt_LPS2 = []

        BBS_LPS3 = []
        ten_LPS3 = []
        thirty_LPS3 = []
        hundred_LPS3 = []
        threehundred_LPS3 = []
        Mt_LPS3 = []

        count = 0

        while count < (len(dayta)):
            BBS.append (dayta[count])
            ten.append (dayta[count + 1])
            thirty.append (dayta[count + 2])
            hundred.append (dayta[count + 3])
            threehundred.append (dayta[count + 4])
            Mt.append (dayta[count + 5])

            BBS_LPS.append (dayta[count + 6])
            ten_LPS.append (dayta[count + 7])
            thirty_LPS.append (dayta[count + 8])
            hundred_LPS.append (dayta[count + 9])
            threehundred_LPS.append (dayta[count + 10])
            Mt_LPS.append (dayta[count + 11])

            BBS_LPS2.append (dayta[count + 12])
            ten_LPS2.append (dayta[count + 13])
            thirty_LPS2.append (dayta[count + 14])
            hundred_LPS2.append (dayta[count + 15])
            threehundred_LPS2.append (dayta[count + 16])
            Mt_LPS2.append (dayta[count + 17])

            BBS_LPS3.append (dayta[count + 18])
            ten_LPS3.append (dayta[count + 19])
            thirty_LPS3.append (dayta[count + 20])
            hundred_LPS3.append (dayta[count + 21])
            threehundred_LPS3.append (dayta[count + 22])
            Mt_LPS3.append (dayta[count + 23])
            
            count += 24
            
        New_BBS = []
        for a in BBS:
            try:
                z = int (float(a))
                New_BBS.append (float(a))
            except:
                continue

        New_ten = []
        for b in ten:
            try:
                z = int (float(b))
                New_ten.append (float(b))
            except:
                continue

        New_thirty = []
        for c in thirty:
            try:
                z = int (float(c))
                New_thirty.append (float(c))
            except:
                continue

        New_hundred = []
        for d in hundred:
            try:
                z = int (float(d))
                New_hundred.append (float(d))
            except:
                continue

        New_threehundred = []
        for e in threehundred:
            try:
                z = int (float(e))
                New_threehundred.append (float(e))
            except:
                continue

        New_Mt = []
        for f in Mt:
            try:
                z = int (float(f))
                New_Mt.append (float(f))
            except:
                continue

        ################################

        New_BBS_LPS = []
        for g in BBS_LPS:
            try:
                z = int (float(g))
                New_BBS_LPS.append (float(g))
            except:
                continue

        New_ten_LPS = []
        for h in ten_LPS:
            try:
                z = int (float(h))
                New_ten_LPS.append (float(h))
            except:
                continue

        New_thirty_LPS = []
        for i in thirty_LPS:
            try:
                z = int (float(i))
                New_thirty_LPS.append (float(i))
            except:
                continue

        New_hundred_LPS = []
        for j in hundred_LPS:
            try:
                z = int (float(j))
                New_hundred_LPS.append (float(j))
            except:
                continue

        New_threehundred_LPS = []
        for k in threehundred_LPS:
            try:
                z = int (float(k))
                New_threehundred_LPS.append (float(k))
            except:
                continue

        New_Mt_LPS = []
        for l in Mt_LPS:
            try:
                z = int (float(l))
                New_Mt_LPS.append (float(l))
            except:
                continue

        ###################################

        New_BBS_LPS2 = []
        for m in BBS_LPS2:
            try:
                z = int (float(m))
                New_BBS_LPS2.append (float(m))
            except:
                continue

        New_ten_LPS2 = []
        for n in ten_LPS2:
            try:
                z = int (float(n))
                New_ten_LPS2.append (float(n))
            except:
                continue

        New_thirty_LPS2 = []
        for o in thirty_LPS2:
            try:
                z = int (float(o))
                New_thirty_LPS2.append (float(o))
            except:
                continue

        New_hundred_LPS2 = []
        for p in hundred_LPS2:
            try:
                z = int (float(p))
                New_hundred_LPS2.append (float(p))
            except:
                continue

        New_threehundred_LPS2 = []
        for q in threehundred_LPS2:
            try:
                z = int (float(q))
                New_threehundred_LPS2.append (float(q))
            except:
                continue

        New_Mt_LPS2 = []
        for r in Mt_LPS2:
            try:
                z = int (float(r))
                New_Mt_LPS2.append (float(r))
            except:
                continue

        ###################################
        
        New_BBS_LPS3 = []
        for s in BBS_LPS3:
            try:
                z = int (float(s))
                New_BBS_LPS3.append (float(s))
            except:
                continue

        New_ten_LPS3 = []
        for t in ten_LPS3:
            try:
                z = int (float(t))
                New_ten_LPS3.append (float(t))
            except:
                continue

        New_thirty_LPS3 = []
        for u in thirty_LPS3:
            try:
                z = int (float(u))
                New_thirty_LPS3.append (float(u))
            except:
                continue

        New_hundred_LPS3 = []
        for v in hundred_LPS3:
            try:
                z = int (float(v))
                New_hundred_LPS3.append (float(v))
            except:
                continue

        New_threehundred_LPS3 = []
        for w in threehundred_LPS3:
            try:
                z = int (float(w))
                New_threehundred_LPS3.append (float(w))
            except:
                continue

        New_Mt_LPS3 = []
        for x in Mt_LPS3:
            try:
                z = int (float(x))
                New_Mt_LPS3.append (float(x))
            except:
                continue

        ###################################

##        #Two-WAy ANOVA
##
##        sBBS_mean = sum(New_BBS) / len(New_BBS)
##        sten_mean = sum(New_ten) / len(New_ten)
##        sthirty_mean = sum(New_thirty) / len(New_thirty)
##        shundred_mean = sum(New_hundred) / len(New_hundred)
##        sthreehundred_mean = sum(New_threehundred) / len(New_threehundred)
##
##        sBBS_LPS_mean = sum(New_BBS_LPS) / len(New_BBS_LPS)
##        sten_LPS_mean = sum(New_ten_LPS) / len(New_ten_LPS)
##        sthirty_LPS_mean = sum(New_thirty_LPS) / len(New_thirty_LPS)
##        shundred_LPS_mean = sum(New_hundred_LPS) / len(New_hundred_LPS)
##        sthreehundred_LPS_mean = sum(New_threehundred_LPS) / len(New_threehundred_LPS)
##
##        sBBS_LPS2_mean = sum(New_BBS_LPS2) / len(New_BBS_LPS2)
##        sten_LPS2_mean = sum(New_ten_LPS2) / len(New_ten_LPS2)
##        sthirty_LPS2_mean = sum(New_thirty_LPS2) / len(New_thirty_LPS2)
##        shundred_LPS2_mean = sum(New_hundred_LPS2) / len(New_hundred_LPS2)
##        sthreehundred_LPS2_mean = sum(New_threehundred_LPS2) / len(New_threehundred_LPS2)
##
##        sBBS_LPS3_mean = sum(New_BBS_LPS3) / len(New_BBS_LPS3)
##        sten_LPS3_mean = sum(New_ten_LPS3) / len(New_ten_LPS3)
##        sthirty_LPS3_mean = sum(New_thirty_LPS3) / len(New_thirty_LPS3)
##        shundred_LPS3_mean = sum(New_hundred_LPS3) / len(New_hundred_LPS3)
##        sthreehundred_LPS3_mean = sum(New_threehundred_LPS3) / len(New_threehundred_LPS3)
##
##        BBS_var = statistics.variance (New_BBS)
##        ten_var = statistics.variance (New_ten)
##        thirty_var = statistics.variance (New_thirty)
##        hundred_var = statistics.variance (New_hundred)
##        threehundred_var = statistics.variance (New_threehundred)
##
##        BBS_LPS_var = statistics.variance (New_BBS_LPS)
##        ten_LPS_var = statistics.variance (New_ten_LPS)
##        thirty_LPS_var = statistics.variance (New_thirty_LPS)
##        hundred_LPS_var = statistics.variance (New_hundred_LPS)
##        threehundred_LPS_var = statistics.variance (New_threehundred_LPS)
##
##        BBS_LPS2_var = statistics.variance (New_BBS_LPS2)
##        ten_LPS2_var = statistics.variance (New_ten_LPS2)
##        thirty_LPS2_var = statistics.variance (New_thirty_LPS2)
##        hundred_LPS2_var = statistics.variance (New_hundred_LPS2)
##        threehundred_LPS2_var = statistics.variance (New_threehundred_LPS2)
##
##        BBS_LPS3_var = statistics.variance (New_BBS_LPS3)
##        ten_LPS3_var = statistics.variance (New_ten_LPS3)
##        thirty_LPS3_var = statistics.variance (New_thirty_LPS3)
##        hundred_LPS3_var = statistics.variance (New_hundred_LPS3)
##        threehundred_LPS3_var = statistics.variance (New_threehundred_LPS3)
##
##        Control_mean = statistics.mean (New_BBS + New_ten + New_thirty + New_hundred + New_threehundred)
##        LPS_mean = statistics.mean (New_BBS_LPS + New_ten_LPS + New_thirty_LPS + New_hundred_LPS + New_threehundred_LPS)
##        LPS2_mean = statistics.mean (New_BBS_LPS2 + New_ten_LPS2 + New_thirty_LPS2 + New_hundred_LPS2 + New_threehundred_LPS2)
##        LPS3_mean = statistics.mean (New_BBS_LPS3 + New_ten_LPS3 + New_thirty_LPS3 + New_hundred_LPS3 + New_threehundred_LPS3)
##        LPS_dose_mean = statistics.mean(New_BBS + New_ten + New_thirty + New_hundred + New_threehundred + New_BBS_LPS + New_ten_LPS + New_thirty_LPS + New_hundred_LPS + New_threehundred_LPS + New_BBS_LPS2 + New_ten_LPS2 + New_thirty_LPS2 + New_hundred_LPS2 + New_threehundred_LPS2 + New_BBS_LPS3 + New_ten_LPS3 + New_thirty_LPS3 + New_hundred_LPS3 + New_threehundred_LPS3) 
##
##        BBS_mean = statistics.mean (New_BBS + New_BBS_LPS + New_BBS_LPS2 + New_BBS_LPS3)
##        ten_mean = statistics.mean (New_ten + New_ten_LPS + New_ten_LPS2 + New_ten_LPS3)
##        thirty_mean = statistics.mean (New_thirty + New_thirty_LPS + New_thirty_LPS2 + New_thirty_LPS3)
##        hundred_mean = statistics.mean (New_hundred + New_hundred_LPS + New_hundred_LPS2 + New_hundred_LPS3)
##        threehundred_mean = statistics.mean (New_threehundred + New_threehundred_LPS + New_threehundred_LPS2 + New_threehundred_LPS3)
##        Bacteria_dose_mean = LPS_dose_mean
##        
##        df_BBS = len (BBS) - 1
##        df_ten = len (ten) - 1
##        df_thirty = len (thirty) - 1
##        df_hundred = len (hundred) - 1
##        df_threehundred = len (threehundred) - 1
##        df_Mt = len (Mt) - 1
##
##        df_BBS_LPS = len (BBS_LPS) - 1
##        df_ten_LPS = len (ten_LPS) - 1
##        df_thirty_LPS = len (thirty_LPS) - 1
##        df_hundred_LPS = len (hundred_LPS) - 1
##        df_threehundred_LPS = len (threehundred_LPS) - 1
##        df_Mt_LPS = len (Mt_LPS) - 1
##
##        df_BBS_LPS2 = len (BBS_LPS2) - 1
##        df_ten_LPS2 = len (ten_LPS2) - 1
##        df_thirty_LPS2 = len (thirty_LPS2) - 1
##        df_hundred_LPS2 = len (hundred_LPS2) - 1
##        df_threehundred_LPS2 = len (threehundred_LPS2) - 1
##        df_Mt_LPS2 = len (Mt_LPS2) - 1
##
##        df_BBS_LPS3 = len (BBS_LPS3) - 1
##        df_ten_LPS3 = len (ten_LPS3) - 1
##        df_thirty_LPS3 = len (thirty_LPS3) - 1
##        df_hundred_LPS3 = len (hundred_LPS3) - 1
##        df_threehundred_LPS3 = len (threehundred_LPS3) - 1
##        df_Mt_LPS3 = len (Mt_LPS3) - 1
##
##        SS_within = (df_BBS * BBS_var) + (df_ten * ten_var) + (df_thirty * thirty_var) + (df_hundred * hundred_var) + (df_threehundred * threehundred_var) + (df_BBS_LPS * BBS_LPS_var) + (df_ten_LPS * ten_LPS_var) + (df_thirty_LPS * thirty_LPS_var) + (df_hundred_LPS * hundred_LPS_var) + (df_threehundred_LPS * threehundred_LPS_var) + (df_BBS_LPS2 * BBS_LPS2_var) + (df_ten_LPS2 * ten_LPS2_var) + (df_thirty_LPS2 * thirty_LPS2_var) + (df_hundred_LPS2 * hundred_LPS2_var) + (df_threehundred_LPS2 * threehundred_LPS2_var) + (df_BBS_LPS3 * BBS_LPS3_var) + (df_ten_LPS3 * ten_LPS3_var) + (df_thirty_LPS3 * thirty_LPS3_var) + (df_hundred_LPS3 * hundred_LPS3_var) + (df_threehundred_LPS3 * threehundred_LPS3_var) 
##        df_within = df_BBS + df_ten + df_thirty + df_hundred + df_threehundred + df_BBS_LPS + df_ten_LPS + df_thirty_LPS + df_hundred_LPS + df_threehundred_LPS + df_BBS_LPS2 + df_ten_LPS2 + df_thirty_LPS2 + df_hundred_LPS2 + df_threehundred_LPS2 + df_BBS_LPS3 + df_ten_LPS3 + df_thirty_LPS3 + df_hundred_LPS3 + df_threehundred_LPS3
##        MSE_within = SS_within / df_within
##
##        ###Flex?
##        SS_LPS_dose = int(flex) * 5 * (((Control_mean - LPS_dose_mean)**2) + ((LPS_mean - LPS_dose_mean)**2) + ((LPS2_mean - LPS_dose_mean)**2) + ((LPS3_mean - LPS_dose_mean)**2))
##        df_LPS_dose = 3
##        MSE_LPS_dose = SS_LPS_dose/df_LPS_dose
##
##        ###Flex?
##        SS_bacteria_dose = int(flex) * 4 * (((BBS_mean - Bacteria_dose_mean)**2) + ((ten_mean - Bacteria_dose_mean)**2) + ((thirty_mean - Bacteria_dose_mean)**2) + ((hundred_mean - Bacteria_dose_mean)**2) + ((threehundred_mean - Bacteria_dose_mean)**2))
##        df_bacteria_dose = 4
##        MSE_bacteria_dose = SS_bacteria_dose/df_bacteria_dose
##
##        #Yij = replicate mean; Yi = bacteria dose mean; Yj = LPS dose mean; Y = grand mean
##        ###Flex?
##        SS_interaction = int(flex) * (((sBBS_mean-BBS_mean-Control_mean+LPS_dose_mean)**2) + ((sten_mean-ten_mean-Control_mean+LPS_dose_mean)**2) + ((sthirty_mean-thirty_mean-Control_mean+LPS_dose_mean)**2) + ((shundred_mean-hundred_mean-Control_mean+LPS_dose_mean)**2) + ((sthreehundred_mean-threehundred_mean-Control_mean+LPS_dose_mean)**2) + ((sBBS_LPS_mean-BBS_mean-LPS_mean+LPS_dose_mean)**2) + ((sten_LPS_mean-ten_mean-LPS_mean+LPS_dose_mean)**2) + ((sthirty_LPS_mean-thirty_mean-LPS_mean+LPS_dose_mean)**2) + ((shundred_LPS_mean-hundred_mean-LPS_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS_mean-threehundred_mean-LPS_mean+LPS_dose_mean)**2) + ((sBBS_LPS2_mean-BBS_mean-LPS2_mean+LPS_dose_mean)**2) + ((sten_LPS2_mean-ten_mean-LPS2_mean+LPS_dose_mean)**2) + ((sthirty_LPS2_mean-thirty_mean-LPS2_mean+LPS_dose_mean)**2) + ((shundred_LPS2_mean-hundred_mean-LPS2_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS2_mean-threehundred_mean-LPS2_mean+LPS_dose_mean)**2) + ((sBBS_LPS3_mean-BBS_mean-LPS3_mean+LPS_dose_mean)**2) + ((sten_LPS3_mean-ten_mean-LPS3_mean+LPS_dose_mean)**2) + ((sthirty_LPS3_mean-thirty_mean-LPS3_mean+LPS_dose_mean)**2) + ((shundred_LPS3_mean-hundred_mean-LPS3_mean+LPS_dose_mean)**2) + ((sthreehundred_LPS3_mean-threehundred_mean-LPS3_mean+LPS_dose_mean)**2))
##        df_interaction = 12
##        MSE_interaction = SS_interaction/df_interaction
##
##        F_LPS_dose = MSE_LPS_dose / MSE_within
##        F_bacteria_dose = MSE_bacteria_dose / MSE_within
##        F_interaction = MSE_interaction / MSE_within
##
##        p_LPS_dose = stats.f.sf(F_LPS_dose, df_LPS_dose, df_within)
##        p_bacteria_dose = stats.f.sf(F_bacteria_dose, df_bacteria_dose, df_within)
##        p_interaction = stats.f.sf (F_interaction, df_interaction, df_within)
        
        control_results = scipy.stats.f_oneway (New_BBS, New_ten, New_thirty, New_hundred, New_threehundred)
        LPS_results = scipy.stats.f_oneway (New_BBS_LPS, New_ten_LPS, New_thirty_LPS, New_hundred_LPS, New_threehundred_LPS)
        LPS2_results = scipy.stats.f_oneway (New_BBS_LPS2,New_ten_LPS2,New_thirty_LPS2,New_hundred_LPS2,New_threehundred_LPS2)
        LPS3_results = scipy.stats.f_oneway (New_BBS_LPS3,New_ten_LPS3,New_thirty_LPS3,New_hundred_LPS3,New_threehundred_LPS3)
        
        control_t_results = scipy.stats.ttest_ind (New_BBS, New_Mt, equal_var = True)
        LPS_t_results = scipy.stats.ttest_ind (New_BBS_LPS, New_Mt_LPS, equal_var = True)
        LPS2_t_results = scipy.stats.ttest_ind (New_BBS_LPS2, New_Mt_LPS2, equal_var = True)
        LPS3_t_results = scipy.stats.ttest_ind (New_BBS_LPS3, New_Mt_LPS3, equal_var = True)

        print ("========================================================")
        
        print ("")
        print ("Your One-Way ANOVA Values (BBS to 300) Are: ")
        print ("")
        print ("Control F value, P value:", control_results)
        print ("")
        print ("1st LPS Dose F value, P value:", LPS_results)
        print ("")
        print ("2nd LPS Dose F value, P value:", LPS2_results)
        print ("")
        print ("3rd LPS Dose F value, P value:", LPS3_results)
        print ("")
##        
##        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        ###################################
##        control_results = scipy.stats.f_oneway (New_BBS, New_ten, New_thirty, New_hundred, New_threehundred)
##        LPS_results = scipy.stats.f_oneway (New_BBS_LPS, New_ten_LPS, New_thirty_LPS, New_hundred_LPS, New_threehundred_LPS)
##        LPS2_results = scipy.stats.f_oneway (New_BBS_LPS2,New_ten_LPS2,New_thirty_LPS2,New_hundred_LPS2,New_threehundred_LPS2)
##        LPS3_results = scipy.stats.f_oneway (New_BBS_LPS3,New_ten_LPS3,New_thirty_LPS3,New_hundred_LPS3,New_threehundred_LPS3)
        
##        control_t_results = scipy.stats.ttest_ind (New_BBS, New_Mt, equal_var = False)
##        LPS_t_results = scipy.stats.ttest_ind (New_BBS_LPS, New_Mt_LPS, equal_var = False)
##        LPS2_t_results = scipy.stats.ttest_ind (New_BBS_LPS2, New_Mt_LPS2, equal_var = False)
##        LPS3_t_results = scipy.stats.ttest_ind (New_BBS_LPS3, New_Mt_LPS3, equal_var = False)

        ###################################
        
##        print ("========================================================")
##        print ("")
##        print ("Your One-Way ANOVA Values (BBS to 300) Are: ")
##        print ("")
##        print ("Control F value, P value:", control_results)
##        print ("")
##        print ("1st LPS Dose F value, P value:", LPS_results)
##        print ("")
##        print ("2nd LPS Dose F value, P value:", LPS2_results)
##        print ("")
##        print ("3rd LPS Dose F value, P value:", LPS3_results)
##        print ("")

##        print("")
##        print ("Your Two-Way ANOVA Values (BBS to 300) Are: ")
##        print ("")
##        print ("LPS_Dose F and P value: ", "F =", F_LPS_dose, "P =", p_LPS_dose)
##        print ("")
##        print ("Bacteria_Dose F and P value:", "F =", F_bacteria_dose, "P =", p_bacteria_dose)
##        print ("")
##        print ("LPS_Dose and Bacteria_Dose Interaction F and P value:", "F =", F_interaction, "P =", p_interaction)
##        print ("")
        
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("")
        print ("Your T Test Values (BBS vs. M. tub) Are: ")
        print ("")
        print ("Control t-statistic, P value:", control_t_results)
        print ("")
        print ("1st LPS Dose t-statistic, P value:", LPS_t_results)
        print ("")
        print ("2nd LPS Dose t-statistic, P value:", LPS2_t_results)
        print ("")
        print ("3rd LPS Dose t-statistic, P value:", LPS3_t_results)
        print ("")
        print ("========================================================")
        
        Next = input ("Run Post-Hoc Test? <Y/N> ")

        if Next == "Y" or Next == "y" or Next == "Yes" or Next == "yes" or Next == "YES":
            LPS_PCR_posthoc_test (N = flex, BBS = New_BBS, ten = New_ten, thirty = New_thirty, hundred = New_hundred,
                         threehundred = New_threehundred, Mt = New_Mt, BBS_LPS = New_BBS_LPS, ten_LPS = New_ten_LPS,
                         thirty_LPS = New_thirty_LPS, hundred_LPS = New_hundred_LPS, threehundred_LPS =
                         New_threehundred_LPS, Mt_LPS = New_Mt_LPS, BBS_LPS2 = New_BBS_LPS2, ten_LPS2 =
                         New_ten_LPS2, thirty_LPS2 = New_thirty_LPS2, hundred_LPS2 = New_hundred_LPS2,
                         threehundred_LPS2 = New_threehundred_LPS2, Mt_LPS2 = New_Mt_LPS2, BBS_LPS3 =
                         New_BBS_LPS3, ten_LPS3 = New_ten_LPS3, thirty_LPS3 = New_thirty_LPS3, hundred_LPS3 =
                         New_hundred_LPS3, threehundred_LPS3 = New_threehundred_LPS3, Mt_LPS3 = New_Mt_LPS3)

        elif Next == "N" or Next == "n" or Next == "No" or Next == "no" or Next == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")


            elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 

                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")

        else:
            raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

    return;

def PCR_analyzer2 (n, data, LPS_factor):

    print ("========================================================")
    
    Listt = []
    for a in n:
        if a == ",":
            continue
        Listt.append (float(a))

    if LPS_factor == True:
        count = 0
        while count < len (Listt):
            Listt[count] = Listt[count] * 2
            count += 1
        n = sum (Listt)
    
    csv = pandas.read_csv(data, header = None, na_filter = False)

    List_50 = []
    List_51 = []
    count = 0
    while count < ((sum(Listt) * 2) * 12):
        List_50.append (csv.iloc[count, 0])
        List_51.append (csv.iloc[count, 1])
        
        count += 1
        
    List = []
    for a in List_50:
        try:
            y = int (float(a))
            List.append (float(a))
        except:
            List.append (a)
            
    List_2 = []
    for b in List_51:
        try:
            z = int (float(b))
            List_2.append (float(b))
        except:
            List_2.append (b)
            
    List_3 = []
    count_2 = 0
    count_3 = 1
    a = 0
    b = 0
    while count_2 < (len(List) - 12):
        if (count_3 % 13) == 0:
            count_2 += 12
            count_3 = 1
            a = 1
            b = 1
        elif b == 1:
            count_2 += 1
            count_3 += 1
            b = 0
        
        Replicate_1 = List[count_2]
        Replicate_2 = List[12 + count_2]

        if type (List[count_2]) == str: 
            Replicate_Average = List[12 + count_2]  
        elif type (List[12 + count_2]) == str:
            Replicate_Average = List[count_2] 
        else:
            Replicate_Average = (Replicate_1 + Replicate_2) / 2
            
        List_3.append (Replicate_Average)

        if a == 1:
            a = 0
        else:
            count_2 += 1
            count_3 += 1
    
    List_7 = []
    count_8 = 0
    count_9 = 1
    c = 0
    d = 0
    while count_8 < (len(List_2) - 12):
        if (count_9 % 13) == 0:
            count_8 += 12
            count_9 = 1
            c = 1
            d = 1
        elif d == 1:
            count_8 += 1
            count_9 += 1
            d = 0
        
        Replicate_3 = List_2[count_8]
        Replicate_4 = List_2[12 + count_8]

        if type (List_2[count_8]) == str: 
            Replicate_Average_2 = List_2[12 + count_8]  
        elif type (List_2[12 + count_8]) == str:
            Replicate_Average_2 = List_2[count_8]  
        else:
            Replicate_Average_2 = (Replicate_3 + Replicate_4) / 2
            
        List_7.append (Replicate_Average_2)

        if c == 1:
            c = 0
        else:
            count_8 += 1
            count_9 += 1
    
    List_4 = []
    count_4 = 0
    while count_4 < (12 * sum(Listt)):
        C_avg = List_3[count_4]
        B_actin = List_7[count_4]

        if type (List_3[count_4]) == str:
            Normalized = "NaN"
        elif type (List_7[count_4]) == str:
            Normalized = "NaN"
        else:
            Normalized = C_avg - B_actin

        List_4.append (Normalized)

        count_4 += 1

    count_5 = 0
    Normal_List = []
    while count_5 < len (List_4):
        Normal_List.append (float(List_4[count_5]))
        count_5 += 1

    count_set = 0
    count_set2 = 1
    lower = 0
    upper = 0
    Lower_List = [0]   ##################### maybe here? ####################
    Upper_List = [(Listt[0] * 12)]
    Max_List = []
    while count_set < (len (Listt) - 1):
        lower = Lower_List[count_set]
        upper = Upper_List[count_set]

        Lower_List.append (upper)
        Upper_List.append (Upper_List[count_set] + (12 * Listt[count_set2]))
        
        count_set += 1
        count_set2 += 1

    count_max = 0
    while count_max < len (Listt):
        Max = max (Normal_List[int(Lower_List[count_max]):int(Upper_List[count_max])])
        Max_List.append (Max)

        count_max += 1

    count_6 = 0
    count_7 = 1
    count_set3 = 0
    List_5 = []
    while count_6 < len (List_4): 
        if (count_7 % ((12 * Listt[count_set3]) + 1)) == 0:
            count_7 = 1
            count_set3 += 1
        if type (List_4[count_6]) == str:
            Delta = "NaN"
        else:
            Delta = List_4[count_6] - Max_List[count_set3]

        List_5.append (Delta)

        count_6 += 1
        count_7 += 1

    count_10 = 0
    List_6 = []
    while count_10 < len (List_5):
        if type (List_5[count_10]) == str:
            Log2 = "NaN"
        else:
            Log2 = 2**(List_5[count_10] * -1)

        List_6.append (Log2)

        count_10 += 1

    print ("")
    print ("Your Log2 Values Are:")
    print ("")

    count_11 = 0
    while count_11 < len(List_6):
        print (List_6[count_11])
        
        count_11 += 1

    n = sum(Listt)
    
    if LPS_factor == True:
        n = sum(Listt) / 2
    
    print ("")
    print ("========================================================")
    
    step2 = input ("Would You Like to Graph Your Log2 Values? <Y/N> ")

    if step2 == "Y" or step2 == "y" or step2 == "Yes" or step2 == "yes" or step2 == "YES":
        PCR_grapher (N = n, Data = List_6, LPS_factor = LPS_factor)
        
    elif step2 == "N" or step2 == "n" or step2 == "No" or step2 == "no" or step2 == "NO":
        
        answer = input ("Would You Like To Run Stats? <Y/N> ")

        if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
            One_Way_ANOVA (n = n, dayta = List_6, LPS_factor = LPS_factor)

        elif answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO":
            answer = input ("Would You Like to Analyze Another Data Set? <Y/N> ")
            print ("========================================================")

            if answer == "Y" or answer == "y" or answer == "Yes" or answer == "yes" or answer == "YES":
                print ("")
                print ("1: LPS Dose Response")
                print ("")
                print ("2: Control vs. 1 LPS Dose")
                print ("")

                experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

                if experiment != "1" and experiment != "2":
                    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                print ("========================================================")

                print ("")
                print ("1: Analyze Cq Values")
                print ("")
                print ("2: Analyze Log2 Values")
                print ("")

                program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

                if program == "1" and experiment == "2":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "1" and experiment == "1":
                    print ("========================================================")
                    print ("")
                    print ("1: Analyze a Data Set Using 1 Minimum Value")
                    print ("")
                    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
                    print ("")
                    
                    meta_program = input ("Which Method Would You Like to Use? ")
                    
                    if meta_program == "1":
                        print ("========================================================")
                        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
                    elif meta_program == "2":
                        print ("========================================================")
                        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
                        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
                elif program == "2" and experiment == "1":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

                elif program == "2" and experiment == "2":
                    print ("")
                    print ("========================================================")
                    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
                    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
                        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
                        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
                    else:
                        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

            if answer == "N" or answer == "n" or answer == "No" or answer == "no" or answer == "NO": 
                #print ("Thank You For Using The qPCR Analyzer 3000! -EH")
                print ("")
                try:
                    shell = sys.stdout.shell
                except AttributeError:
                    raise RuntimeError("you must run this program in IDLE")

                shell.write("Thank ", "COMMENT")
                shell.write("You ","KEYWORD")
                shell.write("For ","STRING")
                shell.write("Using ","DEFINITION")
                shell.write("The ","BUILTIN")
                shell.write("qPCR ","COMMENT")
                shell.write("Analyzer ","KEYWORD")
                shell.write("3000!","STRING")
                print ("")
                print ("")
                shell.write("-EH","BUILTIN")
                print ("")
                print ("")
                print ("========================================================")
                
    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again") 

    return;

#from tkinter.filedialog import askopenfilename

print ("")
print ("=====================================================================================")
print ("")
print ("INSTRUCTIONS:")
print ("")
print ("1: Data must be saved in a .csv format with cytokine data in the first row&column and beta actin data in the first row and second column")
print ("")
print ("2: Missing data values cen be left blank or filled with any letters you want (e.g. NaN, nan, askdhi, puppies)")
print ("")
print ("3: File location input must be in a file directory format. If the file directory is unknown, further instructions for MacOS can be found here (https://www.switchingtomac.com/tutorials/osx/5-ways-to-reveal-the-path-of-a-file-on-macos/) or for Windows here (https://www.wikihow.com/Find-a-File%27s-Path-on-Windows)")
print ("")
print ("4: For statistical analysis, data must be in an alternating control/LPS/control/LPS/etc. format")
print ("")
print ("=====================================================================================")
print ("")

print ("1: LPS Dose Response")
print ("")
print ("2: Control vs. 1 LPS Dose")
print ("")

experiment = input ("Which Type of Experiment Did You Run? (Enter the Number): ")

if experiment != "1" and experiment != "2":
    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

print ("========================================================")

print ("")
print ("1: Analyze Cq Values")
print ("")
print ("2: Analyze Log2 Values")
print ("")

program = input ("Which Data Type Would You Like to Analyze? (Enter the Number): ")

if program == "1" and experiment == "2":
    print ("========================================================")
    print ("")
    print ("1: Analyze a Data Set Using 1 Minimum Value")
    print ("")
    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
    print ("")
    
    meta_program = input ("Which Method Would You Like to Use? ")
    
    if meta_program == "1":
        print ("========================================================")
        PCR_analyzer (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
    elif meta_program == "2":
        print ("========================================================")
        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = False)
    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
elif program == "1" and experiment == "1":
    print ("========================================================")
    print ("")
    print ("1: Analyze a Data Set Using 1 Minimum Value")
    print ("")
    print ("2: Analyze a Combined Data Set Using A Minimum Value From Each Part")
    print ("")
    
    meta_program = input ("Which Method Would You Like to Use? ")
    
    if meta_program == "1":
        print ("========================================================")
        PCR_analyzer (n = float(input ("Number of Replicates? ")), data = input ("File Location? "), LPS_factor = True)
    elif meta_program == "2":
        print ("========================================================")
        #PCR_analyzer2 (data = askopenfilename(), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES): "))
        PCR_analyzer2 (data = input ("File Location? "), n = input ("Number of Replicates in Each Data Set In Order of Appearance in .csv File? (Separate With Commas, NO SPACES; e.g. 3,4,2): "), LPS_factor = True)
    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
elif program == "2" and experiment == "1":
    print ("")
    print ("========================================================")
    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = True)
    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

elif program == "2" and experiment == "2":
    print ("")
    print ("========================================================")
    meta_meta_program = input ("Would You Like to Graph Your Data? <Y/N> ")
    if meta_meta_program == "Y" or meta_meta_program == "y" or meta_meta_program == "Yes" or meta_meta_program == "yes" or meta_meta_program == "YES":
        PCR_grapher_2 (N = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
    elif meta_meta_program == "N" or meta_meta_program == "n" or meta_meta_program == "No" or meta_meta_program == "no" or meta_meta_program == "NO":
        PCR_stats (n = input ("Number of Replicates? "), data = input ("File Location? "), LPS_factor = False)
    else:
        raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")
    
else:
    raise SystemExit ("Input Not Recognized; Re-Run Program and Try Again")

