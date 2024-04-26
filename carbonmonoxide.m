function [rho_carbonmonoxide, mu_carbonmonoxide, k_carbonmonoxide, C_carbonmonoxide] = carbonmonoxide(T_carbonmonoxide)
%CoolProp CO properties at 1.15 bara
%Note viscoscity and thermal conductivity are estiamted at CO2 for now
properties = [...
290	1.35431493491080	1.45246295443268e-05	0.0160195357698005	1042.33260281792
300	1.30901053943365	1.50046454675511e-05	0.0167828229431323	1042.44851466957
310	1.26664869879518	1.54819717752559e-05	0.0175538041984166	1042.64301976749
320	1.22695082457926	1.59564579189455e-05	0.0183316830089891	1042.91999713537
330	1.18967253361790	1.64279699894106e-05	0.0191157238027321	1043.28334270395
340	1.15459853763502	1.68963904354996e-05	0.0199052460647189	1043.73676420873
350	1.12153842549161	1.73616174838262e-05	0.0206996193302685	1044.28362024553
360	1.09032316069170	1.78235643668165e-05	0.0214982588679492	1044.92679873023
370	1.06080215609655	1.82821584350116e-05	0.0223006219085077	1045.66863023360
380	1.03284081752255	1.87373402073839e-05	0.0231062043136862	1046.51083177516
390	1.00631847059096	1.91890623977264e-05	0.0239145376056330	1047.45447677288
400	0.981126602663564	1.96372889439170e-05	0.0247251862966103	1048.49998700770
410	0.957167365245478	2.00819940586912e-05	0.0255377454720240	1049.64714268921
420	0.934352292820636	2.05231613146082e-05	0.0263518385885612	1050.89510699501
430	0.912601202411189	2.09607827715295e-05	0.0271671154523714	1052.24246178882
440	0.891841244744168	2.13948581517349e-05	0.0279832503284335	1053.68725158248
450	0.872006083159663	2.18253940654540e-05	0.0287999399795532	1055.22703317850
460	0.853035180601811	2.22524032878806e-05	0.0296169021068507	1056.85892879829
470	0.834873178423105	2.26759040874996e-05	0.0304338825384292	1058.57968085426
480	0.817469353477061	2.30959196046679e-05	0.0312506322699541	1060.38570685402
490	0.800777142207927	2.35124772787767e-05	0.0320669255450244	1062.27315322583
500	0.784753722272432	2.39256083219003e-05	0.0328825527360505	1064.23794712365
510	0.769359643728745	2.43353472365729e-05	0.0336973192788735	1066.27584550558
520	0.754558503065223	2.47417313751789e-05	0.0345110446817408	1068.38248098305
530	0.740316654366555	2.51448005383768e-05	0.0353235616027427	1070.55340410956
540	0.726602952767263	2.55445966099648e-05	0.0361347149904077	1072.78412192110
550	0.713388526054003	2.59411632256443e-05	0.0369443612826663	1075.07013265640
560	0.700646570874159	2.63345454732022e-05	0.0377523676598323	1077.40695667870
570	0.688352170509270	2.67247896217334e-05	0.0385586113476409	1079.79016369288
580	0.676482131594490	2.71119428776346e-05	0.0393629789667180	1082.21539640676
590	0.665014837522885	2.74960531652199e-05	0.0401653659251632	1084.67839082519
600	0.653930116576930	2.78771689299351e-05	0.0409656758511914	1087.17499339246
610	0.643209123087979	2.82553389622736e-05	0.0417638200630260	1089.70117521578
620	0.632834230145134	2.86306122406233e-05	0.0425597170734522	1092.25304361071
630	0.622788932563787	2.90030377913972e-05	0.0433532921266367	1094.82685121122
640	0.613057758986302	2.93726645649185e-05	0.0441444767650005	1097.41900288382
650	0.603626192126810	2.97395413256470e-05	0.0449332084240967	1100.02606067740
660	0.594480596292507	3.01037165554406e-05	0.0457194300535942	1102.64474703062
670	0.585608151417994	3.04652383686491e-05	0.0465030897626051	1105.27194644585
680	0.576996792939403	3.08241544379344e-05	0.0472841404877249	1107.90470582542
690	0.568635156913559	3.11805119298021e-05	0.0480625396822630	1110.54023365128
700	0.560512529855642	3.15343574489116e-05	0.0488382490252577	1113.17589817469
710	0.552618802828460	3.18857369903140e-05	0.0496112341489621	1115.80922476781
720	0.544944429368495	3.22346958988382e-05	0.0503814643835824	1118.43789257491
730	0.537480386879574	3.25812788349111e-05	0.0511489125181331	1121.05973058704
740	0.530218141165053	3.29255297461661e-05	0.0519135545763548	1123.67271325120
750	0.523149613804656	3.32674918442435e-05	0.0526753696067090	1126.27495571242
760	0.516267152113165	3.36072075862476e-05	0.0534343394855356	1128.86470877608
770	0.509563501445529	3.39447186603674e-05	0.0541904487325188	1131.44035366694
780	0.503031779637243	3.42800659752160e-05	0.0549436843376660	1134.00039665184
790	0.496665453390245	3.46132896524843e-05	0.0556940355990585	1136.54346358379
800	0.490458316433650	3.49444290225411e-05	0.0564414939706819	1139.06829441770
810	0.484404469305550	3.52735226226463e-05	0.0571860529196917	1141.57373774007
820	0.478498300617139	3.56006081974767e-05	0.0579277077925115	1144.05874534888
830	0.472734469673850	3.59257227016901e-05	0.0586664556892029	1146.52236691392
840	0.467107890340148	3.62489023042812e-05	0.0594022953455821	1148.96374474257
850	0.461613716045337	3.65701823945065e-05	0.0601352270225943	1151.38210867153
860	0.456247325837288	3.68895975891765e-05	0.0608652524024876	1153.77677110080
870	0.451004311399613	3.72071817411342e-05	0.0615923744913607	1156.14712218290
880	0.445880464955462	3.75229679487549e-05	0.0623165975276836	1158.49262517693
890	0.440871767988099	3.78369885663215e-05	0.0630379268964204	1160.81281197453
900	0.435974380714576	3.81492752151428e-05	0.0637563690484032	1163.10727880256
910	0.431184632254488	3.84598587952942e-05	0.0644719314246330	1165.37568210504
920	0.426499011440813	3.87687694978771e-05	0.0651846223852021	1167.61773460549
930	0.421914158224429	3.90760368176988e-05	0.0658944511425513	1169.83320154916
940	0.417426855628012	3.93816895662891e-05	0.0666014276987986	1172.02189712357
950	0.413034022208782	3.96857558851769e-05	0.0673055627868853	1174.18368105464
960	0.408732704992916	3.99882632593592e-05	0.0680068678153092	1176.31845537495
970	0.404520072847546	4.02892385309017e-05	0.0687053548162245	1178.42616136000
980	0.400393410259044	4.05887079126172e-05	0.0694010363967029	1180.50677662779
990	0.396350111488802	4.08866970017756e-05	0.0700939256929658	1182.56031239661
1000	0.392387675080075	4.11832307938016e-05	0.0707840363274055	1184.58681089562
1010	0.388503698691489	4.14783336959249e-05	0.0714713823682290	1186.58634292267
1020	0.384695874234793	4.17720295407499e-05	0.0721559782915641	1188.55900554348
1030	0.380961983296130	4.20643415997168e-05	0.0728378389458820	1190.50491992651
1040	0.377299892821737	4.23552925964287e-05	0.0735169795185961	1192.42422930741
1050	0.373707551050404	4.26449047198245e-05	0.0741934155047069	1194.31709707748
1060	0.370182983676396	4.29331996371791e-05	0.0748671626773711	1196.18370499016
1070	0.366724290227727	4.32201985069146e-05	0.0755382370602807	1198.02425147987
1080	0.363329640645824	4.35059219912099e-05	0.0762066549017431	1199.83895008781
1090	0.359997272053624	4.37903902683974e-05	0.0768724326503630	1201.62802798906
1100	0.356725485700104	4.40736230451380e-05	0.0775355869322298	1203.39172461583
1110	0.353512644070109	4.43556395683662e-05	0.0781961345295225	1205.13029037172
1120	0.350357168149123	4.46364586370005e-05	0.0788540923604475	1206.84398543214
1130	0.347257534833401	4.49160986134142e-05	0.0795094774604323	1208.53307862593
1140	0.344212274476488	4.51945774346626e-05	0.0801623069644992	1210.19784639393
1150	0.341219968563861	4.54719126234661e-05	0.0808125980907516	1211.83857181982
1160	0.338279247507914	4.57481212989460e-05	0.0814603681249067	1213.45554372932
1170	0.335388788556113	4.60232201871143e-05	0.0821056344058129	1215.04905585353
1180	0.332547313805569	4.62972256311173e-05	0.0827484143118962	1216.61940605279
1190	0.329753588317802	4.65701536012342e-05	0.0833887252484785	1218.16689559734
1200	0.327006418327803	4.68420197046313e-05	0.0840265846359200	1219.69182850135
1210	0.324304649541978	4.71128391948764e-05	0.0846620098985349	1221.19451090708
1220	0.321647165519831	4.73826269812130e-05	0.0852950184542374	1222.67525051605
1230	0.319032886134648	4.76513976375999e-05	0.0859256277048744	1224.13435606430
1240	0.316460766108699	4.79191654115176e-05	0.0865538550272050	1225.57213683892
1250	0.313929793618788	4.81859442325466e-05	0.0871797177644904	1226.98890223330
1260	0.311438988968243	4.84517477207194e-05	0.0878032332186565	1228.38496133854
1270	0.308987403321678	4.87165891946530e-05	0.0884244186429980	1229.76062256873
1280	0.306574117499093	4.89804816794635e-05	0.0890432912353918	1231.11619331788
1290	0.304198240826086	4.92434379144680e-05	0.0896598681319897	1232.45197964637
1300	0.301858910037164	4.95054703606790e-05	0.0902741664013627	1233.76828599505
1310	0.299555288229298	4.97665912080946e-05	0.0908862030390716	1235.06541492505
1320	0.297286563863074	5.00268123827892e-05	0.0914959949626377	1236.34366688166
1330	0.295051949808919	5.02861455538099e-05	0.0921035590068912	1237.60333998054
1340	0.292850682436050	5.05446021398829e-05	0.0927089119196748	1238.84472981479
1350	0.290682020741935	5.08021933159334e-05	0.0933120703578824	1240.06812928146
1360	0.288545245520172	5.10589300194257e-05	0.0939130508838130	1241.27382842607
1370	0.286439658564829	5.13148229565255e-05	0.0945118699618206	1242.46211430398
1380	0.284364581909386	5.15698826080914e-05	0.0951085439552445	1243.63327085726
1390	0.282319357098542	5.18241192354976e-05	0.0957030891236021	1244.78757880613
1400	0.280303344491247	5.20775428862944e-05	0.0962955216200290	1245.92531555378
];

tempinterp = properties(:,1);
interp = interp1q(tempinterp, properties, T_carbonmonoxide);

rho_carbonmonoxide = interp(2); %density, kg/m^3
mu_carbonmonoxide = interp(3); %viscocity, Pa*s
k_carbonmonoxide = interp(4); %conductive heat trans coef (W/mK)
C_carbonmonoxide = interp(5); %conductive heat trans coef (J/kgK)
