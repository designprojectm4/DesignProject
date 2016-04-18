function [molar_density,mass_density] = density(T,P)
%Finds molar and mass densities of components as functions of T[KELVIN]
%   Inputs:
%       Temperature [K]
%       Pressure P=[PH;P_L][Pa]
PH=P(1);
%   Outputs: 
%       molar density[rho_H;rho_G;rho_A;rho_P;rho_W;rho_POH] kmol/m3
%       mass density[rho_H;rho_G;rho_A;rho_P;rho_W;rho_POH] kg/m3
%--------------------------------------------------------------------------
%Gas Molar Density: For Hydrogen
Z=solvePR(T,PH);
R=8.314459848;%J K-1 mol-1
rho_H=PH/(Z*R*T);%mol/m3
rho_H=rho_H*1e-3;%kmol/m3
molar_density=[rho_H;rho_G(T);rho_A(T);rho_P(T);rho_W(T);rho_POH(T)];
RFM=[2;92.09382;74.08;76.09;18.01528;60.09502];%[H2;G;A;P;POH]
mass_density=bsxfun(@times,molar_density,RFM);
%--------------------------------------------------------------------------
%Liquid Molar Density
%ASSUMED TO BE INDEPENDENT OF PRESSURE
%UNITS: kmol/m3
%GLYCEROL
function rho_G=rho_G(T)
       T_table=493:1:800;
        rho_table=[10.39067628
10.38349995
10.3763089
10.36910308
10.36188246
10.35464699
10.34739663
10.34013134
10.33285106
10.32555577
10.3182454
10.31091992
10.30357928
10.29622342
10.28885229
10.28146585
10.27406404
10.26664681
10.25921411
10.25176588
10.24430206
10.2368226
10.22932743
10.22181651
10.21428976
10.20674713
10.19918855
10.19161396
10.18402329
10.17641648
10.16879347
10.16115417
10.15349853
10.14582646
10.13813791
10.13043279
10.12271103
10.11497255
10.10721728
10.09944514
10.09165605
10.08384992
10.07602668
10.06818624
10.06032851
10.05245341
10.04456085
10.03665074
10.02872299
10.02077751
10.0128142
10.00483296
9.996833711
9.988816341
9.980780755
9.972726849
9.964654521
9.956563666
9.948454179
9.940325951
9.932178873
9.924012834
9.915827721
9.90762342
9.899399814
9.891156785
9.882894214
9.87461198
9.866309959
9.857988025
9.849646053
9.841283914
9.832901476
9.824498608
9.816075174
9.807631039
9.799166065
9.79068011
9.782173033
9.773644689
9.765094931
9.756523612
9.74793058
9.739315684
9.730678767
9.722019672
9.713338241
9.704634312
9.69590772
9.687158301
9.678385884
9.669590299
9.660771374
9.651928932
9.643062795
9.634172783
9.625258712
9.616320397
9.607357649
9.598370278
9.589358091
9.58032089
9.571258478
9.562170651
9.553057207
9.543917937
9.534752631
9.525561077
9.516343059
9.507098357
9.497826749
9.488528011
9.479201914
9.469848227
9.460466714
9.45105714
9.441619261
9.432152835
9.422657613
9.413133344
9.403579773
9.393996642
9.384383689
9.374740649
9.365067253
9.355363227
9.345628295
9.335862177
9.326064587
9.316235237
9.306373836
9.296480086
9.286553686
9.276594332
9.266601714
9.256575518
9.246515427
9.236421118
9.226292263
9.216128531
9.205929585
9.195695083
9.18542468
9.175118025
9.16477476
9.154394524
9.143976952
9.13352167
9.123028302
9.112496465
9.101925769
9.091315822
9.080666223
9.069976565
9.059246438
9.048475424
9.037663098
9.026809029
9.015912782
9.004973912
8.993991969
8.982966497
8.971897031
8.9607831
8.949624226
8.938419924
8.9271697
8.915873054
8.904529477
8.893138453
8.881699458
8.870211958
8.858675413
8.847089273
8.835452978
8.823765963
8.81202765
8.800237452
8.788394776
8.776499014
8.764549553
8.752545767
8.740487019
8.728372665
8.716202046
8.703974495
8.691689331
8.679345865
8.666943392
8.654481197
8.641958554
8.629374721
8.616728946
8.604020462
8.591248489
8.578412233
8.565510885
8.552543623
8.539509607
8.526407986
8.51323789
8.499998433
8.486688714
8.473307814
8.459854796
8.446328706
8.432728572
8.419053402
8.405302185
8.39147389
8.377567466
8.36358184
8.34951592
8.335368588
8.321138706
8.306825111
8.292426616
8.277942011
8.263370057
8.248709491
8.233959022
8.219117331
8.204183071
8.189154863
8.174031299
8.15881094
8.143492312
8.12807391
8.11255419
8.096931576
8.081204453
8.065371167
8.049430025
8.033379293
8.017217194
8.000941907
7.984551564
7.968044251
7.951418002
7.934670806
7.917800595
7.900805246
7.88368258
7.866430358
7.849046281
7.831527985
7.813873039
7.796078944
7.778143129
7.760062944
7.741835667
7.72345849
7.704928521
7.68624278
7.667398192
7.648391586
7.629219689
7.609879122
7.590366395
7.5706779
7.550809907
7.530758558
7.510519859
7.490089675
7.469463722
7.448637557
7.427606574
7.406365992
7.384910843
7.363235969
7.341335998
7.319205354
7.296838219
7.274228535
7.251369983
7.228255969
7.204879607
7.181233695
7.157310702
7.133102738
7.108601535
7.083798416
7.05868427
7.033249514
7.007484064
6.981377293
6.954917987
6.928094302
6.900893706
6.873302927
6.845307884
6.816893618
6.788044209
6.758742687
6.728970932
6.698709556
6.66793777
6.636633256
6.604771971
6.572327971
6.539273187
6.505577164
6.471206768
6.436125836
6.40029477
6.363670052
6.326203673
6.28784245];rho_table=transpose(rho_table); 
rho_G=interp1(T_table,rho_table,T,'linear','extrap');
    end
%ACETOL
function rho_A=rho_A(T)
       T_table=493:1:800;
        rho_table=[9.573069792
9.552707761
9.532159608
9.511421872
9.490491001
9.469363345
9.448035157
9.426502587
9.40476168
9.382808369
9.360638471
9.338247686
9.315631587
9.292785617
9.269705087
9.246385165
9.222820872
9.199007076
9.174938484
9.150609638
9.126014902
9.101148459
9.0760043
9.050576215
9.024857784
8.998842365
8.972523086
8.945892831
8.918944227
8.891669631
8.864061119
8.836110466
8.807809131
8.779148241
8.750118569
8.720710518
8.690914094
8.660718881
8.630114029
8.599088208
8.567629586
8.535725799
8.503363912
8.470530381
8.437211005
8.403390896
8.369054412
8.334185107
8.298765674
8.262777872
8.226202457
8.189019097
8.151206282
8.112741225
8.073599747
8.033756151
7.993183081
7.951851371
7.909729855
7.866785171
7.822981529
7.778280452
7.73264047
7.686016783
7.638360855
7.589619953
7.53973661
7.488647991
7.436285139
7.382572094
7.327424826
7.270749947
7.212443171
7.152387382
7.090450278
7.026481403
6.960308393
6.891732148
6.820520546
6.746400099
6.66904468
6.588059954
6.502961298
6.413141556
6.317822277
6.215976582
6.106200638
5.986484186
5.853762119
5.702919892
5.524118826
5.292640684
4.847063857
4.336814336
4.021503643
3.76405898
3.550569702
3.3712811
3.219142425
3.088909053
2.976539995
2.87869212
2.791828079
2.694660205
2.69039704
2.687070344
2.684398476
2.68224943
2.680518555
2.679122691
2.677995617
2.677084504
2.676347134
2.676644012
2.677796337
2.679040246
2.680348526
2.681698306
2.683070933
2.684451528
2.685828426
2.687192635
2.688537335
2.689857455
2.69114932
2.692410355
2.693638854
2.694833789
2.695994658
2.697121367
2.698214131
2.6992734
2.700299801
2.701294087
2.702257103
2.703189755
2.704092987
2.704967764
2.70581506
2.706635843
2.70743107
2.70820168
2.708948592
2.709672696
2.710374858
2.711055912
2.711716665
2.712357891
2.712980335
2.713584714
2.714171712
2.714741988
2.715296171
2.715834863
2.716358643
2.716868061
2.717363645
2.717845901
2.718315311
2.718772337
2.719217421
2.719650985
2.720073434
2.720485153
2.720886515
2.721277871
2.721659563
2.722031914
2.722395236
2.722749827
2.723095972
2.723433946
2.723764011
2.724086419
2.724401412
2.724709221
2.72501007
2.725304172
2.725591733
2.72587295
2.726148013
2.726417105
2.7266804
2.726938068
2.727190271
2.727437166
2.727678903
2.727915628
2.728147479
2.728374593
2.728597098
2.72881512
2.72902878
2.729238194
2.729443475
2.729644732
2.729842069
2.730035588
2.730225387
2.73041156
2.730594198
2.730773391
2.730949223
2.731121778
2.731291136
2.731457373
2.731620566
2.731780786
2.731938104
2.732092589
2.732244306
2.73239332
2.732539692
2.732683483
2.732824752
2.732963554
2.733099945
2.733233978
2.733365706
2.733495178
2.733622443
2.733747549
2.733870541
2.733991466
2.734110365
2.734227282
2.734342257
2.734455332
2.734566544
2.734675931
2.734783531
2.734889379
2.734993511
2.735095959
2.735196758
2.735295939
2.735393535
2.735489575
2.735584089
2.735677106
2.735768656
2.735858764
2.735947459
2.736034767
2.736120712
2.736205321
2.736288617
2.736370625
2.736451366
2.736530865
2.736609143
2.736686222
2.736762123
2.736836866
2.736910472
2.73698296
2.737054349
2.737124659
2.737193907
2.737262112
2.73732929
2.737395459
2.737460637
2.737524838
2.73758808
2.737650377
2.737711746
2.737772201
2.737831757
2.737890428
2.737948228
2.738005171
2.738061271
2.738116541
2.738170993
2.738224641
2.738277496
2.738329571
2.738380877
2.738431427
2.738481232
2.738530303
2.738578651
2.738626286
2.73867322
2.738719462
2.738765022
2.738809912
2.73885414
2.738897715
2.738940648
2.738982948
2.739024623
2.739065682
2.739106134
2.739145988
2.739185252
2.739223934
2.739262043
2.739299585
2.739336569
2.739373003
2.739408894
2.739444249
2.739479076
2.739513382
2.739547173
2.739580457];rho_table=transpose(rho_table); 
rho_A=interp1(T_table,rho_table,T,'linear','extrap');
    end
%PG
function rho_P=rho_P(T)
       T_table=493:1:800;
        rho_table=[10.42594582
10.40958355
10.39309476
10.37647762
10.35973027
10.34285081
10.32583729
10.30868773
10.29140011
10.27397237
10.25640239
10.23868802
10.22082706
10.20281726
10.18465632
10.16634189
10.14787156
10.12924289
10.11045336
10.0915004
10.07238137
10.0530936
10.03363432
10.01400072
9.994189902
9.974198919
9.954024738
9.933664254
9.913114288
9.892371576
9.871432775
9.850294452
9.828953089
9.807405075
9.7856467
9.763674159
9.741483539
9.719070824
9.696431885
9.673562476
9.650458232
9.627114662
9.603527145
9.579690924
9.555601099
9.531252624
9.506640298
9.481758759
9.456602476
9.431165745
9.405442671
9.379427177
9.353112977
9.326493575
9.299562255
9.272312067
9.244735818
9.216826054
9.188575058
9.159974824
9.131017046
9.101693104
9.071994043
9.041910555
9.01143296
8.980551184
8.949254732
8.917532668
8.885373582
8.852765564
8.819696172
8.78615239
8.7521206
8.717586534
8.682535228
8.646950974
8.610817267
8.574116743
8.536831107
8.498941076
8.460426285
8.4212652
8.381435023
8.340911581
8.299669199
8.257680569
8.214916588
8.17134619
8.126936148
8.08165084
8.03545202
7.988298502
7.940145835
7.890945926
7.840646583
7.789191005
7.736517173
7.682557133
7.627236154
7.570471713
7.512172286
7.45223588
7.390548231
7.326980591
7.261386989
7.193600722
7.12342994
7.050651881
6.975005315
6.896180396
6.813804793
6.72742425
6.636474573
6.540239914
6.437788096
6.327865504
6.208715022
6.077735113
5.930767393
5.760351756
5.550175614
5.243781056
4.605037617
4.23245006
3.934804323
3.692511371
3.492327385
3.324986209
3.183818065
3.06390269
2.961517009
2.873712137
2.797686771
2.718401538
2.715561423
2.713370126
2.711627
2.710238484
2.709130995
2.70824655
2.707539373
2.70697327
2.706618815
2.70750973
2.708529311
2.709650558
2.71084837
2.712100707
2.713388912
2.71469758
2.71601423
2.717328906
2.71863378
2.71992278
2.721191278
2.722435812
2.72365386
2.724843657
2.726004039
2.727134316
2.728234177
2.729303604
2.73034281
2.731352181
2.732332238
2.733283604
2.734206971
2.735103084
2.73597272
2.736816679
2.737635766
2.738430789
2.739202549
2.739951835
2.740679422
2.741386068
2.742072507
2.742739457
2.743387608
2.744017631
2.744630171
2.745225851
2.745805269
2.746369004
2.746917608
2.747451616
2.747971537
2.748477863
2.748971065
2.749451595
2.749919888
2.750376358
2.750821405
2.751255412
2.751678746
2.75209176
2.752494791
2.752888165
2.753272193
2.753647174
2.754013395
2.754371131
2.754720648
2.755062199
2.755396031
2.755722376
2.756041463
2.756353507
2.756658717
2.756957296
2.757249437
2.757535325
2.757815141
2.758089057
2.75835724
2.75861985
2.758877041
2.759128964
2.759375762
2.759617574
2.759854533
2.760086769
2.760314406
2.760537566
2.760756364
2.760970914
2.761181322
2.761387696
2.761590136
2.76178874
2.761983605
2.76217482
2.762362477
2.76254666
2.762727453
2.762904938
2.763079192
2.763250291
2.763418308
2.763583316
2.763745383
2.763904576
2.76406096
2.764214598
2.764365552
2.764513881
2.764659642
2.764802892
2.764943684
2.765082073
2.765218108
2.765351841
2.765483319
2.76561259
2.7657397
2.765864692
2.765987612
2.7661085
2.766227397
2.766344345
2.766459381
2.766572544
2.766683871
2.766793396
2.766901156
2.767007185
2.767111515
2.767214179
2.767315209
2.767414636
2.767512488
2.767608797
2.767703589
2.767796893
2.767888737
2.767979147
2.768068147
2.768155765
2.768242024
2.768326949
2.768410562
2.768492887
2.768573947
2.768653763
2.768732357
2.768809749
2.76888596
2.76896101
2.769034919
2.769107705
2.769179387
2.769249984
2.769319513
2.769387992
2.769455438
2.769521867
2.769587296
2.769651742
2.769715219
2.769777743
2.769839328
2.769899991
2.769959745
2.770018604
2.770076582
2.770133693
2.770189949
2.770245364];rho_table=transpose(rho_table); 
rho_P=interp1(T_table,rho_table,T,'linear','extrap');
    end
%WATER
function rho_W=rho_W(T)
       T_table=493:1:800;
        rho_table=[33.66905941
33.6064397
33.5434581
33.48011002
33.4163908
33.35229566
33.2878197
33.22295792
33.1577052
33.09205631
33.0260059
32.95954845
32.89267837
32.82538989
32.75767711
32.68953398
32.6209543
32.55193173
32.48245973
32.41253163
32.34214055
32.27127944
32.19994107
32.12811801
32.05580261
31.98298702
31.90966317
31.83582276
31.76145723
31.68655779
31.61111539
31.53512069
31.4585641
31.38143569
31.30372525
31.22542224
31.14651578
31.06699464
30.9868472
30.90606148
30.82462507
30.74252516
30.65974845
30.57628121
30.49210918
30.40721761
30.32159118
30.235214
30.14806956
30.06014072
29.97140965
29.8818578
29.79146586
29.70021372
29.60808041
29.51504405
29.42108179
29.32616976
29.230283
29.13339537
29.0354795
28.93650669
28.8364468
28.73526817
28.63293751
28.52941975
28.42467794
28.31867309
28.21136398
28.10270701
27.99265599
27.8811619
27.76817268
27.65363291
27.53748352
27.41966143
27.30009917
27.17872443
27.05545953
26.93022088
26.80291829
26.67345426
26.54172303
26.40760964
26.27098874
26.13172318
25.98966242
25.84464061
25.69647429
25.54495962
25.38986911
25.23094749
25.06790673
24.90041982
24.72811288
24.55055514
24.36724602
24.17759808
23.98091423
23.7763567
23.56290352
23.33928623
23.10389722
22.8546471
22.58873387
22.30224687
21.9894325
21.64118299
21.24141916
20.75604212
20.07708646
18.57998601
17.60954956
16.77030433
16.03825073
15.39492345
14.82588369
14.31967887
13.86710868
13.46069886
13.09431374
12.76287049
12.46212282
12.1884961
11.93895994
11.71092879
11.50218355
11.31080918
11.13514471
10.97374277
10.82533661
10.68881299
10.56319026
10.44759928
10.34126833
10.24350998
10.15371016
10.07131885
9.995842156
9.926835497
9.863897767
9.806666207
9.754812164
9.708037007
9.666068833
9.62865939
9.595581327
9.566625684
9.541599516
9.520323803
9.502631245
9.488364421
9.477374255
9.469518875
9.464603167
9.464088855
9.464065292
9.464044312
9.464025625
9.464008973
9.46399413
9.463980893
9.463969085
9.463958546
9.463949137
9.463940733
9.463933224
9.463926512
9.463920511
9.463915142
9.463910338
9.463906037
9.463914909
9.463941801
9.463978931
9.464028358
9.464092132
9.464172237
9.46427054
9.464388759
9.464528429
9.464690892
9.464877278
9.46508851
9.465325306
9.465588186
9.46587748
9.466193347
9.466535786
9.466904652
9.46729967
9.467720454
9.468166516
9.468637287
9.469132124
9.469650324
9.470191136
9.47075377
9.471337409
9.47194121
9.472564319
9.473205873
9.473865007
9.474540858
9.475232569
9.475939294
9.476660199
9.477394466
9.478141297
9.478899909
9.479669542
9.480449458
9.481238942
9.4820373
9.482843863
9.483657985
9.484479045
9.485306443
9.486139604
9.486977975
9.487821029
9.488668256
9.489519172
9.490373313
9.491230234
9.492089512
9.492950743
9.493813542
9.494677542
9.495542393
9.496407764
9.497273339
9.498138817
9.499003915
9.49986836
9.500731898
9.501594286
9.502455294
9.503314704
9.50417231
9.50502792
9.505881349
9.506732425
9.507580984
9.508426873
9.509269948
9.510110072
9.51094712
9.511780971
9.512611513
9.513438641
9.514262259
9.515082274
9.515898601
9.516711161
9.517519881
9.518324694
9.519125535
9.519922347
9.520715076
9.521503674
9.522288097
9.523068302
9.523844254
9.524615919
9.525383268
9.526146274
9.526904913
9.527659165
9.528409013
9.529154441
9.529895437
9.530631991
9.531364095
9.532091743
9.532814932
9.533533661
9.534247928
9.534957737
9.535663091
9.536363994
9.537060453
9.537752477
9.538440073
9.539123254
9.53980203
9.540476414
9.541146419
9.541812061
9.542473355
9.543130318
9.543782965
9.544431317
9.545075391
9.545715206
9.546350783
9.546982142
9.547609304
9.548232291
9.548851125
9.549465828
9.550076422
9.550682932
9.551285381
9.551883793
9.552478191
9.5530686
9.553655045];rho_table=transpose(rho_table); 
rho_W=interp1(T_table,rho_table,T,'linear','extrap');
    end
%PROPANOL
function rho_POH=rho_POH(T)
       T_table=493:1:800;
        rho_table=[7.497769373
7.452509218
7.406482624
7.359653365
7.311982113
7.263426053
7.213938422
7.163468018
7.111958532
7.059347839
7.005567113
6.950539764
6.894180161
6.836392056
6.777066646
6.716080162
6.653290806
6.588534859
6.521621621
6.452326749
6.380383303
6.305469554
6.227191765
6.145059429
6.058448417
5.966543664
5.868245994
5.762010763
5.64554535
5.515177207
5.364306648
5.17849475
4.907991059
4.382175576
4.11165037
3.882805319
3.686856749
3.517283385
3.369114339
3.23844744
3.12208988
3.017219051
2.920846218
2.827898153
2.721989008
2.714629184
2.70840709
2.703139232
2.698673356
2.694882743
2.691661641
2.688921584
2.686588406
2.684599813
2.682903391
2.681454969
2.68021727
2.679158793
2.68003553
2.681236017
2.682448879
2.683663523
2.684871618
2.686066711
2.687243882
2.688399452
2.689530739
2.690635855
2.691713545
2.692763053
2.693784012
2.694776362
2.695740278
2.696676116
2.697584369
2.698465629
2.699320564
2.70014989
2.700954359
2.701734741
2.702491814
2.703226357
2.703939141
2.704630927
2.705302459
2.705954465
2.70658765
2.707202698
2.707800273
2.708381012
2.708945533
2.709494427
2.710028266
2.710547596
2.711052944
2.711544814
2.712023691
2.712490037
2.712944298
2.713386901
2.713818253
2.714238746
2.714648755
2.71504864
2.715438745
2.715819401
2.716190923
2.716553616
2.71690777
2.717253664
2.717591567
2.717921734
2.718244413
2.718559841
2.718868243
2.71916984
2.71946484
2.719753444
2.720035847
2.720312235
2.720582787
2.720847675
2.721107065
2.721361116
2.721609982
2.721853812
2.722092746
2.722326924
2.722556477
2.722781533
2.723002214
2.72321864
2.723430926
2.723639181
2.723843512
2.724044024
2.724240814
2.724433979
2.724623613
2.724809805
2.724992641
2.725172207
2.725348583
2.725521848
2.725692078
2.725859347
2.726023725
2.726185283
2.726344087
2.726500203
2.726653691
2.726804615
2.726953033
2.727099002
2.727242577
2.727383814
2.727522763
2.727659476
2.727794003
2.727926391
2.728056686
2.728184934
2.728311179
2.728435463
2.728557828
2.728678314
2.728796961
2.728913806
2.729028886
2.729142238
2.729253897
2.729363897
2.729472272
2.729579052
2.729684271
2.729787959
2.729890145
2.729990859
2.73009013
2.730187984
2.730284449
2.730379552
2.730473317
2.73056577
2.730656935
2.730746836
2.730835497
2.730922939
2.731009185
2.731094256
2.731178174
2.731260959
2.731342631
2.731423209
2.731502714
2.731581163
2.731658574
2.731734966
2.731810357
2.731884762
2.731958199
2.732030684
2.732102234
2.732172862
2.732242585
2.732311418
2.732379375
2.732446471
2.732512719
2.732578133
2.732642726
2.732706512
2.732769503
2.732831711
2.73289315
2.732953831
2.733013765
2.733072965
2.733131441
2.733189205
2.733246267
2.733302638
2.733358329
2.733413349
2.733467709
2.733521418
2.733574487
2.733626923
2.733678737
2.733729938
2.733780535
2.733830535
2.733879949
2.733928783
2.733977047
2.734024748
2.734071894
2.734118494
2.734164554
2.734210082
2.734255086
2.734299572
2.734343548
2.73438702
2.734429997
2.734472483
2.734514486
2.734556012
2.734597068
2.73463766
2.734677793
2.734717474
2.734756709
2.734795504
2.734833864
2.734871795
2.734909303
2.734946392
2.734983068
2.735019337
2.735055203
2.735090672
2.735125749
2.735160438
2.735194744
2.735228672
2.735262227
2.735295413
2.735328235
2.735360697
2.735392804
2.735424559
2.735455968
2.735487034
2.735517761
2.735548153
2.735578214
2.735607949
2.735637361
2.735666453
2.73569523
2.735723695
2.735751852
2.735779704
2.735807255
2.735834509
2.735861467
2.735888135
2.735914515
2.735940611
2.735966425
2.735991961
2.736017222
2.736042211
2.736066931
2.736091385
2.736115576
2.736139506
2.736163179
2.736186597
2.736209764
2.736232681
2.736255352
2.736277779
2.736299965
2.736321912
2.736343623
2.7363651];rho_table=transpose(rho_table); 
rho_POH=interp1(T_table,rho_table,T,'linear','extrap');
    end
end
