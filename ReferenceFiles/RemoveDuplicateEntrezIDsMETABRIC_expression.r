# Gene expression data was manually curated by searching for the duplicated Entrez ID at https://www.ncbi.nlm.nih.gov/gene/, and keeping the row with Hugo_symbol corresponding to the "Official symbol" in the search result. If the official symbol does not match any of the listed Hugo_symbols, matches found in the "Also known as:" field are used (if several of the possible gene names are found in the "Also known as:" field, whichever comes first in the list is kept). If no matches are found at this point, all rows for the duplicated gene are removed.
# This curation is based on METABRIC-data downloaded on 28.03.2019, files last modified 19.02.2019. If different files are used, this manual curation may need to be repeated, or at least verified.


rowsToRemove <- c(1907, #Entrez 84268
                  291,  #Entrez 172
                  1620, #Entrez 389634; all incorrect
                  1967, #Entrez 389634; all incorrect
                  17384, #Entrez 389634; all incorrect
                  2104, #Entrez 51312
                  1360, #Entrez 25873
                  2338, #Entrez 51237
                  2408, #Entrez 57460
                  2086, #Entrez 387590
                  3340, #Entrez 728656
                  3494, #Entrez 22897
                  2068, #Entrez 23521
                  3605, #Entrez 23521
                  8943, #Entrez 23521
                  9191, #Entrez 23521
                  9382, #Entrez 23521
                  419, #Entrez 5959
                  3115, #Entrez 374666
                  4017, #Entrez 374666
                  18259, #Entrez 374666
                  3444, #Entrez 5379
                  4041, #Entrez 5379
                  252, #Entrez 730094
                  3035, #Entrez 6638
                  4294, #Entrez 6638
                  11415, #Entrez 6638
                  3772, #Entrez 548593
                  981, #Entrez 552900
                  1264, #Entrez 81557
                  4437, #Entrez 653203; all incorrect
                  4767, #Entrez 653203; all incorrect
                  3938, #Entrez 55073
                  4828, #Entrez 80895
                  533, #Entrez 4745
                  4984, #Entrez 64853
                  5125, #Entrez 7258
                  141, #Entrez 4584
                  5401, #Entrez 26013
                  2543, #Entrez 284861; all incorrect
                  5432, #Entrez 284861; all incorrect
                  4240, #Entrez 407738
                  5647, #Entrez 57035
                  85, #Entrez 283970; all incorrect
                  5648, #Entrez 283970; all incorrect
                  5816, #Entrez 115352
                  1072, #Entrez 727956
                  6584, #Entrez 729096
                  2781, #Entrez 60312
                  40, #Entrez 6130
                  6971, #Entrez 6130
                  7024, #Entrez 23530
                  524, #Entrez 28517; all incorrect
                  7031, #Entrez 28517; all incorrect
                  11962, #Entrez 28517; all incorrect
                  12776, #Entrez 28517; all incorrect
                  13901, #Entrez 28517; all incorrect
                  16694, #Entrez 28517; all incorrect
                  16866, #Entrez 28517; all incorrect
                  7048, #Entrez 332
                  2355, #Entrez 92105
                  2136, #Entrez 161725
                  1035, #Entrez 285097; all incorrect
                  7127, #Entrez 285097; all incorrect
                  600, #Entrez 402483
                  7197, #Entrez 402483
                  7211, #Entrez 63895
                  7239, #Entrez 114817
                  10153, #Entrez 114817
                  3551, #Entrez 5026
                  7406, #Entrez 51119
                  7413, #Entrez 260294
                  1989, #Entrez 338799; all incorrect
                  7426,  #Entrez 338799; all incorrect
                  7508, #Entrez 4830
                  6074, #Entrez 57492
                  7545, #Entrez 57492
                  5288, #Entrez 54540
                  3526, #Entrez 117584
                  396, #Entrez 729264
                  7462, #Entrez 390538
                  8019, #Entrez 3206
                  5914, #Entrez 284800
                  8090, #Entrez 23380
                  8349, #Entrez 80832
                  8360, #Entrez 23266
                  1280, #Entrez 349196
                  8403, #Entrez 349196
                  11960, #Entrez 349196
                  8487, #Entrez 56895
                  5139, #Entrez 374650; all incorrect
                  8610, #Entrez 374650; all incorrect
                  16229, #Entrez 374650; all incorrect
                  5830, #Entrez 284942; all incorrect
                  8639, #Entrez 284942; all incorrect
                  18029, #Entrez 284942; all incorrect
                  8713, #Entrez 65065
                  7585, #Entrez 51561
                  8760, #Entrez 51561
                  9542, #Entrez 51561
                  17697, #Entrez 51561
                  8890, #Entrez 27185
                  8605, #Entrez 100288805; all incorrect
                  9036, #Entrez 100288805; all incorrect
                  9172, #Entrez 84818
                  4029, #Entrez 8386
                  7838, #Entrez 837
                  4966, #Entrez 23150
                  9341, #Entrez 159163
                  9397, #Entrez 6189
                  9366, #Entrez 51750
                  9469, #Entrez 728588
                  9580, #Entrez 6757
                  7894, #Entrez 653361
                  9841, #Entrez 2790
                  1462, #Entrez 115557
                  9864, #Entrez 55672
                  9880, #Entrez 152519
                  9888, #Entrez 7341
                  9911, #Entrez 9859
                  2064, #Entrez 57234; all incorrect
                  9942, #Entrez 57234; all incorrect
                  5100, #Entrez 11039
                  10065, #Entrez 8123
                  10096, #Entrez 645426
                  11524, #Entrez 645426
                  10222, #Entrez 63915
                  1638, #Entrez 100288974; all incorrect
                  10228, #Entrez 100288974; all incorrect
                  1875, #Entrez 4046
                  2976, #Entrez 1196
                  10362, #Entrez 84629
                  8444, #Entrez 3838
                  3096, #Entrez 6124
                  10740, #Entrez 91409
                  10829, #Entrez 6167
                  9344, #Entrez 6625
                  5206, #Entrez 10734
                  798, #Entrez 56117
                  10973, #Entrez 5325
                  10016, #Entrez 55668; all incorrect
                  10981, #Entrez 55668; all incorrect
                  11015, #Entrez 4493
                  10293, #Entrez 55125
                  11223, #Entrez 728441
                  8264, #Entrez 401357; all incorect
                  11228, #Entrez 401357; all incorect
                  4129, #Entrez 642477; all incorrect
                  11229, #Entrez 642477; all incorrect
                  2657, #Entrez 5046
                  11375, #Entrez 653784
                  11420, #Entrez 643008
                  7182, #Entrez 1198
                  10692, #Entrez 347127; all incorrect
                  11534, #Entrez 347127; all incorrect
                  11557, #Entrez 55695
                  11580, #Entrez 1270
                  1099, #Entrez 3543
                  11695, #Entrez 441425
                  7852, #Entrez 23042
                  4701, #Entrez 92002
                  1596, #Entrez 83658
                  4495, #Entrez 9374
                  13054, #Entrez 9374
                  12046, #Entrez 11130
                  6459, #Entrez 284001
                  1985, #Entrez 6150
                  7617, #Entrez 81554
                  11033, #Entrez 6147
                  5627, #Entrez 26873
                  12492, #Entrez 96626
                  12500, #Entrez 55027
                  12525, #Entrez 56981
                  12536, #Entrez 100133941
                  12644, #Entrez 10521
                  12658, #Entrez 23085
                  7918, #Entrez 9503
                  8382, #Entrez 143162
                  7283, #Entrez 5036
                  12880, #Entrez 56001
                  12935, #Entrez 22907
                  836, #Entrez 595101; all incorrect
                  13008, #Entrez 595101; all incorrect
                  9363, #Entrez 727849; all incorrect
                  13015, #Entrez 727849; all incorrect
                  17286, #Entrez 727849; all incorrect
                  13045, #Entrez 6125
                  13065, #Entrez 221223
                  13070, #Entrez 9790
                  12249, #Entrez 266553
                  13096, #Entrez 642968
                  13129, #Entrez 55917
                  11152, #Entrez 100132565
                  13258, #Entrez 113177
                  11492, #Entrez 92270
                  13338, #Entrez 54681
                  5952, #Entrez 10528
                  13489, #Entrez 84908
                  13508, #Entrez 10016
                  13533, #Entrez 3811
                  1247, #Entrez 2209
                  5742, #Entrez 5143
                  13667, #Entrez 728096
                  13715, #Entrez 55871
                  13760, #Entrez 116447
                  7101, #Entrez 7366
                  2538, #Entrez 7499
                  6457, #Entrez 908
                  14076, #Entrez 440279
                  7970, #Entrez 10015
                  12785, #Entrez 400464
                  14159, #Entrez 23507
                  14168, #Entrez 8531
                  14230, #Entrez 7381
                  2450, #Entrez 55180
                  14368, #Entrez 440248
                  16278, #Entrez 440248
                  13444, #Entrez 53917
                  8253, #Entrez 55100
                  12278, #Entrez 85415
                  2670, #Entrez 1325
                  10410, #Entrez 8505
                  10492, #Entrez 9946
                  4809, #Entrez 6137
                  14710, #Entrez 642265
                  14724, #Entrez 119032
                  14878, #Entrez 6170
                  14990, #Entrez 148823
                  15026, #Entrez 9701
                  14739, #Entrez 56997
                  13594, #Entrez 8227
                  7073, #Entrez 23336
                  5428, #Entrez 283788
                  15138, #Entrez 10450
                  15239, #Entrez 283755
                  15414, #Entrez 114800
                  13449, #Entrez 541471; all incorrect
                  15469, #Entrez 541471; all incorrect
                  15522, #Entrez 55559
                  15554, #Entrez 221262
                  15645, #Entrez 2969
                  14591, #Entrez 390705
                  15878, #Entrez 7818
                  15937, #Entrez 8510
                  11620, #Entrez 391627
                  16133, #Entrez 667
                  16170, #Entrez 445347
                  12091, #Entrez 100506084
                  5479, #Entrez 253725
                  16299, #Entrez 220064
                  11319, #Entrez 23284
                  6362, #Entrez 3800
                  16520, #Entrez 5303
                  3281, #Entrez 643837; all incorrect
                  16611, #Entrez 643837; all incorrect
                  7740, #Entrez 2966
                  16692, #Entrez 10625
                  16703, #Entrez 100287704
                  16884, #Entrez 7367
                  16949, #Entrez 514
                  9111, #Entrez 65996; all incorrect
                  16966, #Entrez 65996; all incorrect
                  5817, #Entrez 9639
                  14127, #Entrez 2678
                  17031, #Entrez 653479
                  17047, #Entrez 285596
                  17060, #Entrez 7531
                  13095, #Entrez 145781
                  15857, #Entrez 3514
                  17419, #Entrez 8780
                  17503, #Entrez 7335
                  17620, #Entrez 58505
                  17651, #Entrez 63826
                  15813, #Entrez 5826
                  17732, #Entrez 114780
                  17803, #Entrez 4760
                  1582, #Entrez 3126
                  17969, #Entrez 1565
                  17982, #Entrez 373863
                  18023, #Entrez 10659
                  18059, #Entrez 120863
                  18086, #Entrez 9659
                  13321, #Entrez 10301
                  18143, #Entrez 26206
                  6538, #Entrez 54442
                  6233, #Entrez 57794
                  10050, #Entrez 100506144; all incorrect
                  18205, #Entrez 100506144; all incorrect
                  12460, #Entrez 23049
                  16752, #Entrez 4090
                  18280, #Entrez 100128385
                  7014, #Entrez 7568
                  18116, #Entrez 644662; this record has been withdrawn
                  18379) #Entrez 644662; this record has been withdrawn
