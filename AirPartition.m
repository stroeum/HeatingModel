function eta = AirPartition(E, h, i)
%
% eta = AirPartition(E, h, i)
%
% AirPartition.m  returns the fractions of energy deposited into different
% excitations and the "fast heating" process in air (N2: 78.11%, O2:20.91%,
% Ar: 0.98%).
%
%	eta - fraction of energy deposited into the process given by i
%	E   - electric field (V/m), which can be a scalar or column vector
%	h   - altitude in km (h>150km is intepreted as neutral density in cm-3)
%	i   - integer parameter specifying fractions to be calculated:
%		1   elastic collisions of O2
%		2   elastic collisions of N2
%		3   elastic collisions of Ar
%		4   excitation of rotational energy levels of O2
%		5   excitation of rotational energy levels of N2
%		6   excitation of vibrational energy levels of O2
%		7   excitation of vibrational energy levels of N2 (eta_V)
%		8   excitation of electronic energy levels of O2
%		9   excitation of electronic energy levels of N2
%		10  excitation of electronic energy levels of Ar
%		11  ionization of O2
%		12  ionization of N2
%		13  ionization of Ar
%		14  attachment to O2
%		15  fast heating process (eta_T)
%			eta_T = eta_ela_O2 + eta_ela_N2 + eta_ela_Ar +
%             			eta_rot_O2 + eta_rot_N2 + eta_vib_O2 +
%             			0.3*(eta_el_O2 + eta_el_N2 + eta_el_Ar)
%
% The calculation are based on BOLSIG+. References:
%
%	* http://www.cpat.ups-tlse.fr/BOLSIG-Electron-Boltzmann-equation
%	* G.J.M. Hagelaar and L.C. Pitchford, "Solving the Boltzmann equation
%     to obtain electron transport coefficients and rate coefficients for
%     fluid models", Plasma Sources Sci. Technol. 14 722-733 (2005)
%
% We use the US Standard Atmosphere. We replaced original 2.5e19cm-3 at
% 0 km altitude by 2.688e19 cm-3 which is our reference number density at
% temperature 273 K.
%
% This MATLAB function is provided as is and we do not guarantee or
% assume any responsibility for the accuracy, completeness,
% or usefulness of any information provided herein and derived from
% subsequent calculations.
%
% Please send your comments to Jeremy A Riousset (riousset@psu.edu)
% and Victor P Pasko (vpasko@psu.edu)
%
% Created by Jeremy A Riousset and Victor P Pasko, April 21, 2009
%

if i<1 || i>15
    error('Unknown process');
end

%% Conversion to Td
% We use the US Standard Atmosphere. We replaced original 2.5e19cm-3 at 0 km
% altitude by 2.688e19 cm-3 which is our reference number density at temperature 273 K.
% * statm(1,:) altitude in km
% * statm(2,:) neutral density in cm^-3

statm=[
    0e+0	 2.688e+19
    5e+0	 1.53e+19
    1e+1	 8.59e+18
    1.5e+1	 4.05e+18
    2e+1	 1.85e+18
    2.5e+1	 8.33e+17
    3e+1	 3.83e+17
    3.5e+1	 1.76e+17
    4e+1	 8.31e+16
    4.5e+1	 4.088e+16
    5e+1	 2.13e+16
    5.5e+1	 1.181e+16
    6e+1	 6.33e+15
    6.4e+1	 3.93e+15
    6.8e+1	 2.39e+15
    7.2e+1	 1.39e+15
    7.6e+1	 7.72e+14
    8e+1	 4.03e+14
    8.4e+1	 1.99e+14
    8.8e+1	 9.48e+13
    9.2e+1	 4.37e+13
    9.6e+1	 2.07e+13
    1e+2	 1.04e+13
    1.08e+2	 3.18e+12
    1.14e+2	 1.43e+12
    1.2e+2	 6.61e+11
    1.26e+2	 3.4e+11
    1.32e+2	 1.91e+11
    1.4e+2	 9.7e+10
    1.5e+2	 4.92e+10
    ];
if(0<=h && h<150), %means input is altitude in km
    N = 10^(interp1q(statm(:,1),log10(statm(:,2)),h));
else%means input is neutral density
    N = h;
end %if

%% Output
% The data table used for reference contains the following data;
% # E (Td)
% # nu.ela.O2
% # nu.ela.N2 (momentum)
% # nu.ela.Ar (momentum)
% # nu.rot.O2
% # nu.rot.N2
% # nu.vib.O2
% # nu.vib.N2 (eta.V)
% # nu.el.O2
% # nu.el.N2
% # nu.el.Ar
% # nu.iz.O2
% # nu.iz.N2
% # nu.iz.Ar
% # nu.att.O2
% # eta.T

%% Method
% It is based on BOLSIG+ calculations with:
% * 78.11% N2
% * 20.91% O2
% * 0.98% Ar
% * T.g = 273K
% * n.g = 2.688e19 cm^-3 (used to convert Td into V/m)

T.g   = 273; %_K
n.g   = statm(2,1); %_cm^-3
table = [0.0372	0.00127133	0.0150935	0.000120412	4.17111e-05	0.941435	0	0	0	0	0	0	0	0	0.0420376	0.957962
    0.04179	0.00120914	0.0142572	0.000106066	7.51926e-05	0.918484	3.28234e-12	0	0	0	0	0	0	0	0.065868	0.934132
    0.04694	0.001152	0.0134903	9.33956e-05	0.00012797	0.886292	1.34972e-10	0	0	0	0	0	0	0	0.0988441	0.901156
    0.05273	0.00109567	0.0127414	8.20538e-05	0.000205217	0.844002	9.11643e-10	0	0	0	0	0	0	0	0.141874	0.858126
    0.05924	0.00103671	0.0119747	7.17404e-05	0.000309551	0.79164	5.7058e-09	8.74261e-12	0	0	0	0	0	0	0.194968	0.805032
    0.06654	0.000974736	0.0111722	6.23025e-05	0.000440597	0.730629	3.19216e-08	9.81165e-11	0	0	0	0	0	0	0.256721	0.743279
    0.07475	0.000908657	0.0103449	5.36719e-05	0.000592237	0.663402	1.59984e-07	1.02738e-09	0	0	0	0	0	0	0.324699	0.675301
    0.08397	0.000842982	0.00952093	4.588e-05	0.000755205	0.593409	7.28443e-07	9.81084e-09	0	0	0	0	0	0	0.395425	0.604575
    0.09432	0.000781563	0.00875218	3.88862e-05	0.000916883	0.524126	3.07507e-06	8.48161e-08	0	0	0	0	0	0	0.465381	0.534619
    0.106	0.00073272	0.00811448	3.27485e-05	0.00106756	0.458713	1.22447e-05	6.48295e-07	0	0	0	0	0	0	0.531327	0.468672
    0.119	0.000706228	0.00770783	2.74129e-05	0.00120253	0.399779	3.21654e-05	0	0	0	0	0	0	0	0.590545	0.409455
    0.1337	0.000732635	0.00781347	2.28206e-05	0.00133416	0.349677	0.000213122	2.71736e-05	0	0	0	0	0	0	0.64018	0.359793
    0.1502	0.000853916	0.00881693	1.88165e-05	0.00150835	0.310375	0.000978142	0.000149437	0	0	0	0	0	0	0.6773	0.322551
    0.1687	0.00119542	0.0118369	1.51677e-05	0.00185253	0.2858	0.00455955	0.000745656	0	0	0	0	0	0	0.693995	0.305259
    0.1895	0.00191601	0.0182692	1.17842e-05	0.0025636	0.280495	0.0179293	0.00296649	0	0	0	0	0	0	0.675849	0.321185
    0.2129	0.00263256	0.0245787	9.57535e-06	0.00333017	0.283001	0.042759	0.00700049	6.12663e-18	0	0	0	0	0	0.636689	0.356311
    0.2391	0.00299321	0.0276492	8.55331e-06	0.00378557	0.279024	0.0709122	0.0114879	8.42828e-14	0	0	0	0	0	0.604139	0.384373
    0.2686	0.0031619	0.0289771	7.99829e-06	0.00405775	0.270282	0.101757	0.0162169	1.25693e-13	0	0	0	0	0	0.57554	0.408244
    0.3018	0.0032459	0.0295199	7.63903e-06	0.00423513	0.259453	0.135298	0.0213433	9.52476e-13	0	0	0	0	0	0.546898	0.431759
    0.339	0.00327602	0.0295762	7.37207e-06	0.00435845	0.246943	0.171875	0.0265989	1.0461e-11	0	0	0	0	0	0.517365	0.456036
    0.3808	0.00327858	0.0293982	7.17616e-06	0.00444605	0.233849	0.210834	0.03201	1.13008e-10	0	0	0	0	0	0.486177	0.481813
    0.4277	0.0032589	0.0289939	7.01573e-06	0.00450576	0.219832	0.252454	0.0373285	1.00527e-09	0	0	0	0	0	0.45362	0.509051
    0.4805	0.00323622	0.0285687	6.91712e-06	0.00453318	0.206339	0.294908	0.0426731	7.11564e-09	0	0	0	0	0	0.419735	0.537592
    0.5397	0.00317479	0.0277918	6.79951e-06	0.00454078	0.191469	0.340577	0.0470792	4.36533e-08	0	0	0	0	0	0.385361	0.56756
    0.6063	0.00312319	0.0270866	6.76838e-06	0.00451582	0.177181	0.385807	0.0515834	2.22002e-07	0	0	0	0	0	0.350696	0.597721
    0.6811	0.00305296	0.0262307	6.74984e-06	0.00446757	0.162811	0.431774	0.0553787	9.21821e-07	0	0	0	0	0	0.316277	0.628344
    0.7651	0.00298316	0.0253551	6.80982e-06	0.00438979	0.148712	0.476468	0.0588352	3.61189e-06	0	0	0	0	0	0.283247	0.657916
    0.8594	0.00290409	0.0244012	6.9206e-06	0.00430076	0.134923	0.521501	0.061659	1.21918e-05	0	0	0	0	0	0.250292	0.68804
    0.9654	0.00281841	0.0233626	7.08199e-06	0.00419318	0.121425	0.564611	0.063631	3.53075e-05	0	0	0	0	0	0.219916	0.716428
    1.084	0.00275132	0.0224649	7.4174e-06	0.00407087	0.109142	0.603477	0.0657121	9.44903e-05	0	0	0	0	0	0.19228	0.741941
    1.218	0.0026863	0.0215415	7.86182e-06	0.00394579	0.0973905	0.640373	0.0675582	0.000226358	0	0	0	0	0	0.16627	0.766013
    1.368	0.00264505	0.0207813	8.52137e-06	0.00382098	0.0869268	0.671775	0.0703489	0.000487969	0	0	0	0	0	0.143205	0.786104
    1.537	0.00264899	0.0203082	9.55638e-06	0.00369847	0.0778597	0.696072	0.0757916	0.00100383	0	0	0	0	0	0.122607	0.800898
    1.727	0.00266809	0.0199378	1.0745e-05	0.00358104	0.0701253	0.713296	0.0836232	0.00176875	0	0	0	0	0	0.104989	0.81015
    1.94	0.0027526	0.019974	1.24706e-05	0.00346248	0.0640145	0.718749	0.0986054	0.00308502	0	0	0	0	0	0.0893449	0.80989
    2.179	0.00283812	0.0200609	1.42129e-05	0.00334438	0.0590837	0.715316	0.118564	0.00469581	0	0	0	0	0	0.0760826	0.802066
    2.448	0.00294542	0.0202944	1.62559e-05	0.00320695	0.0550937	0.698531	0.14904	0.0069331	0	0	0	0	0	0.063939	0.782168
    2.75	0.00305408	0.0206201	1.82212e-05	0.00305161	0.0524224	0.669264	0.187935	0.00941196	0	0	0	0	0	0.0542227	0.751254
    3.089	0.00312701	0.0207654	1.99882e-05	0.00286873	0.0500276	0.628957	0.236654	0.0121523	0	0	0	0	0	0.0454277	0.709411
    3.47	0.00315941	0.0207451	2.1437e-05	0.00265372	0.047923	0.578433	0.294591	0.0149062	0	0	0	0	0	0.0375666	0.657408
    3.897	0.0031291	0.0203968	2.22559e-05	0.00242427	0.0459625	0.523452	0.35611	0.0172512	0	0	0	0	0	0.0312517	0.600563
    4.378	0.00303769	0.0197313	2.25109e-05	0.00217489	0.0438977	0.464396	0.421854	0.0191502	0	0	0	0	0	0.0257357	0.539005
    4.918	0.00289049	0.018761	2.21728e-05	0.00192036	0.0417353	0.404986	0.488199	0.0204359	0	0	0	0	0	0.0210498	0.476446
    5.524	0.00269249	0.0175211	2.1304e-05	0.00166925	0.0394523	0.34781	0.552679	0.0210395	0	0	0	0	0	0.0171146	0.415478
    6.206	0.0024639	0.0160992	2.00134e-05	0.0014336	0.0371035	0.295395	0.612636	0.0209706	0	0	0	0	0	0.0138781	0.358806
    6.971	0.00221832	0.0145655	1.84309e-05	0.00121788	0.0347407	0.248509	0.66715	0.0203351	0	0	0	0	0	0.0112444	0.307371
    7.831	0.00196777	0.0129932	1.6682e-05	0.00102386	0.0324058	0.207032	0.716242	0.0192602	0	0	0	0	0	0.00905853	0.261217
    8.796	0.00172249	0.0114476	1.48647e-05	0.000852625	0.03013	0.171102	0.759549	0.0178921	0	0	0	0	0	0.00728887	0.220638
    9.881	0.00149018	0.00997399	1.30658e-05	0.000705453	0.0279478	0.140666	0.797005	0.0163261	0	0	0	0	0	0.00587263	0.185694
    11.1	0.0012755	0.0086043	1.13557e-05	0.000579997	0.0258604	0.115075	0.829202	0.0146748	0	0	0	0	0	0.00471613	0.155809
    12.47	0.00108399	0.00736925	9.78365e-06	0.000474746	0.0239063	0.0937553	0.856562	0.0130527	0	0	0	0	0	0.00378593	0.130515
    14.01	0.000913981	0.00626153	8.35437e-06	0.000386401	0.0220528	0.0760668	0.879814	0.0114788	1.55031e-10	0	0	0	0	0.00301729	0.109133
    15.73	0.000766341	0.00529431	7.09058e-06	0.000313507	0.0203573	0.0615187	0.899301	0.0100256	9.44751e-09	0	0	0	0	0.00241631	0.091265
    17.67	0.000640095	0.00446059	5.99089e-06	0.000254262	0.0188057	0.0498135	0.915366	0.00870961	1.00541e-07	0	0	0	0	0.00194425	0.0765931
    19.85	0.000531193	0.003737	5.03371e-06	0.00020505	0.0173508	0.0401134	0.928952	0.00752885	7.97403e-07	7.06165e-13	9.67035e-13	0	0	0.00157551	0.0642014
    22.3	0.000440228	0.0031273	4.2237e-06	0.000164977	0.0160431	0.032244	0.940214	0.00651396	4.88941e-06	1.14468e-11	2.59937e-11	0	0	0.00124293	0.0539795
    25.05	0.000364855	0.00261897	3.55113e-06	0.000132973	0.0148463	0.0259794	0.949333	0.00568095	2.43863e-05	1.27833e-10	3.5671e-10	0	0	0.00101543	0.0456577
    28.14	0.000301175	0.00218935	2.9874e-06	0.000106444	0.0137689	0.0208036	0.956833	0.00508627	0.000100803	1.07333e-09	3.45046e-09	0	0	0.000807339	0.0387285
    31.61	0.000250004	0.00184103	2.5514e-06	8.53653e-05	0.0127832	0.0167073	0.962382	0.00493094	0.000353123	7.04058e-09	2.57309e-08	3.06478e-11	2.29368e-13	0.000664813	0.0332547
    35.51	0.000208546	0.00155932	2.23999e-06	6.83574e-05	0.0118801	0.0134471	0.965606	0.00557185	0.00106852	3.71477e-08	1.52515e-07	8.6652e-10	1.45429e-11	0.000587392	0.0291578
    39.89	0.000176005	0.00133772	2.07289e-06	5.44759e-05	0.0110451	0.0108632	0.965415	0.00768787	0.0028529	1.62589e-07	7.42604e-07	9.95045e-09	1.91593e-10	0.000565144	0.0266408
    44.81	0.000152106	0.0011765	2.08659e-06	4.32694e-05	0.0102768	0.00890086	0.959649	0.0123499	0.00683009	6.03902e-07	3.03915e-06	8.27675e-08	1.70897e-09	0.000615903	0.0263058
    50.34	0.000135625	0.00106514	2.30053e-06	3.40324e-05	0.00950156	0.00744342	0.945541	0.0208001	0.0147054	1.92495e-06	1.05686e-05	5.29985e-07	1.14843e-08	0.00075868	0.0288343
    56.54	0.000126279	0.00100042	2.73419e-06	2.64742e-05	0.00872797	0.00647228	0.919526	0.0343037	0.0287877	5.35645e-06	3.1806e-05	2.74547e-06	6.20846e-08	0.000986396	0.0352852
    63.52	0.00012383	0.000979603	3.39427e-06	2.04198e-05	0.00792913	0.0059559	0.877494	0.0540434	0.0520429	1.32941e-05	8.46189e-05	1.18366e-05	2.7817e-07	0.00129785	0.0468421
    71.35	0.000125324	0.000979627	4.19005e-06	1.53965e-05	0.00708003	0.00572591	0.819021	0.0790804	0.0860611	2.93077e-05	0.000198166	4.25548e-05	1.0346e-06	0.00163575	0.0634817
    80.15	0.000129787	0.000993778	5.05504e-06	1.14686e-05	0.00620498	0.00573483	0.743908	0.10826	0.132183	5.86e-05	0.000417565	0.000131752	3.29864e-06	0.00195786	0.0852303
    90.03	0.000134005	0.000999254	5.84829e-06	8.33168e-06	0.00530734	0.00580406	0.657486	0.138172	0.188618	0.000106405	0.000793128	0.000349262	8.9606e-06	0.00220783	0.110328
    101.1	0.000136884	0.000990911	6.48044e-06	6.00102e-06	0.00445294	0.00585078	0.566163	0.165499	0.252141	0.000178233	0.00137977	0.000821114	2.14985e-05	0.00235171	0.13679
    113.6	0.00013787	0.000967637	6.93451e-06	4.25226e-06	0.00364883	0.00581905	0.472904	0.188941	0.320878	0.000279632	0.00223848	0.00173256	4.60844e-05	0.00239535	0.163614
    127.6	0.000136082	0.000924093	7.16293e-06	2.96603e-06	0.00292187	0.00565748	0.385747	0.205862	0.389187	0.000412406	0.00340204	0.00331834	8.93215e-05	0.00233269	0.188288
    143.4	0.000131384	0.000864527	7.15125e-06	2.04933e-06	0.00232127	0.00534287	0.3112	0.21503	0.451479	0.000574802	0.00487693	0.00582835	0.000158147	0.00218268	0.208795
    161	0.000125293	0.000798044	6.98387e-06	1.44299e-06	0.00182727	0.00494693	0.247472	0.218131	0.50743	0.000767513	0.00669114	0.00955665	0.000260292	0.00198628	0.225604
    180.9	0.000118383	0.000729127	6.72217e-06	9.91357e-07	0.00140848	0.00449386	0.192908	0.216602	0.556893	0.000989965	0.00887101	0.0148038	0.00040333	0.00177188	0.239103
    203.2	0.000110949	0.000661132	6.3741e-06	6.95639e-07	0.00108745	0.00401093	0.149993	0.210824	0.596764	0.00123418	0.0113974	0.0217678	0.00059134	0.00155074	0.248524
    228.3	0.00010353	0.000596303	5.98787e-06	4.92058e-07	0.000838247	0.00352509	0.116522	0.202145	0.62772	0.00149494	0.0142581	0.0306247	0.000826888	0.00133863	0.254478
    256.4	9.65722e-05	0.000537149	5.59628e-06	3.45215e-07	0.00063885	0.00306331	0.0893721	0.192015	0.651232	0.00176833	0.0174735	0.0415385	0.00111182	0.00114724	0.257846
    288	9.0121e-05	0.000483928	5.20829e-06	2.48827e-07	0.000487981	0.00263394	0.0685918	0.180854	0.666852	0.00204472	0.0210125	0.054525	0.00144226	0.00097695	0.258626
    323.6	8.42489e-05	0.00043649	4.83602e-06	1.77324e-07	0.000372348	0.00224283	0.0526597	0.169308	0.675563	0.00231433	0.0248506	0.0695198	0.00181364	0.000829545	0.257297
    363.5	7.90451e-05	0.000394978	4.48813e-06	1.29082e-07	0.000283567	0.00189687	0.0402983	0.15788	0.678264	0.00257215	0.0289699	0.0864325	0.00221958	0.000704152	0.254274
    408.3	7.44609e-05	0.000358608	4.16646e-06	9.33168e-08	0.000215672	0.00159323	0.0308348	0.14681	0.675674	0.00280944	0.0333409	0.105036	0.00264896	0.000599328	0.249834
    458.6	7.04618e-05	0.000326848	3.87213e-06	6.90773e-08	0.000164606	0.00133079	0.0236832	0.136225	0.66856	0.00301811	0.0379075	0.125105	0.0030935	0.000512085	0.244237
    515.2	6.70087e-05	0.000299139	3.60618e-06	5.03504e-08	0.000124827	0.00110548	0.0180852	0.126307	0.657884	0.00319538	0.0426339	0.146312	0.00354254	0.000440398	0.237816
    578.7	6.40689e-05	0.000275	3.36704e-06	3.73945e-08	9.50936e-05	0.00091422	0.0138823	0.11707	0.644195	0.00333654	0.0474543	0.168341	0.00398624	0.000381922	0.230732
    650.1	6.16058e-05	0.000253893	3.15375e-06	2.77324e-08	7.24023e-05	0.000752872	0.010658	0.108524	0.628163	0.00344045	0.0523433	0.190977	0.0044148	0.000334605	0.223182
    730.3	5.95808e-05	0.00023539	2.9655e-06	2.06372e-08	5.51114e-05	0.000617464	0.00819286	0.100646	0.610435	0.0035089	0.0572433	0.213884	0.00482356	0.000296616	0.215347
    820.4	5.8002e-05	0.000219225	2.80127e-06	1.55014e-08	4.20128e-05	0.000504609	0.00631449	0.0934451	0.591601	0.00353901	0.0621083	0.236695	0.00520427	0.00026663	0.207402
    921.5	5.67612e-05	0.000204728	2.65797e-06	1.15918e-08	3.20515e-05	0.00041056	0.00487658	0.0868373	0.571982	0.00353784	0.0669303	0.259334	0.00555194	0.000243166	0.199414
    1035	5.59261e-05	0.000192226	2.53618e-06	8.78296e-09	2.44808e-05	0.000333256	0.0037743	0.0808352	0.551933	0.00350601	0.0716039	0.281648	0.0058664	0.000225413	0.191491
    1163	5.54766e-05	0.000180921	2.43413e-06	6.62435e-09	1.87329e-05	0.000269627	0.0029306	0.0753463	0.53183	0.00344678	0.0762073	0.303353	0.00614604	0.00021239	0.183714
    1306	5.53523e-05	0.000171063	2.35013e-06	5.01124e-09	1.43541e-05	0.00021744	0.00228122	0.070383	0.511951	0.00336669	0.0806308	0.324332	0.00639105	0.000203569	0.176171
    1467	5.56224e-05	0.00016224	2.28672e-06	3.81831e-09	1.09958e-05	0.000174811	0.00177789	0.0658418	0.492485	0.00326804	0.0849193	0.344501	0.00660246	0.000198501	0.168884
    1648	5.6235e-05	0.000154451	2.24208e-06	2.91426e-09	8.44829e-06	0.000140265	0.00139075	0.0616925	0.473564	0.00315722	0.0890631	0.363789	0.00678459	0.000196942	0.161886
    1852	5.72155e-05	0.000147453	2.21475e-06	2.22481e-09	6.50291e-06	0.000112293	0.00109129	0.0578635	0.455265	0.00303417	0.0930021	0.382283	0.00693727	0.000198612	0.155174
    2080	5.85613e-05	0.000141263	2.20646e-06	1.70804e-09	4.99926e-06	8.96352e-05	0.00085575	0.0543205	0.437751	0.00290518	0.0967601	0.399841	0.00706696	0.000203638	0.14879
    2336	6.02735e-05	0.000135832	2.21811e-06	1.30949e-09	3.86188e-06	7.14415e-05	0.000674575	0.0510078	0.420969	0.00277204	0.100359	0.416557	0.00717465	0.000212139	0.142698
    2625	6.23599e-05	0.000131069	2.25052e-06	1.00542e-09	2.98253e-06	5.68558e-05	0.00053195	0.0478685	0.405115	0.00263829	0.103744	0.432355	0.00726714	0.000224627	0.136942
    2948	6.48561e-05	0.000126865	2.30466e-06	7.7609e-10	2.30556e-06	4.51683e-05	0.000419857	0.0448793	0.390037	0.00250541	0.106956	0.447382	0.00733685	0.000241387	0.131468
    3312	6.77649e-05	0.000123294	2.38286e-06	5.99989e-10	1.786e-06	3.57774e-05	0.000332032	0.0420011	0.375872	0.00237582	0.109977	0.461544	0.00740278	0.000263801	0.126306
    3720	7.11098e-05	0.000120238	2.48859e-06	4.62613e-10	1.38274e-06	2.82717e-05	0.000262478	0.0392039	0.362526	0.00225054	0.112805	0.474984	0.00745186	0.000292205	0.121418
    ];
%% Convert E
% * to match the requirements of the interpolation
% * Convert E to Td using the density at altitude h

S = size(E);
if(S(1)==1)
    E = E';
end

LoBnd = table(1,1);
HiBnd = table(end,1);
E_Td = E/N/(1e-17*1e2);

if min(E_Td)<LoBnd
    warning('MATLAB:paramAmbiguous','Reduced field < 3.72e-2 Td, results may be inaccurate');
end
if max(E_Td)>HiBnd
    warning('MATLAB:paramAmbiguous','Reduced field > 3.72e+3 Td, results may be inaccurate');
end

nu = interp1q(log10(table(:,1)),table(:,i+1),log10(E_Td));
nu(isfinite(nu)==0)=0;
eta = ...
    (E_Td< LoBnd) .* table(1,i+1)       + ...
    (E_Td>=LoBnd) .* (E_Td<=HiBnd) .*nu + ...
    (E_Td> HiBnd) .* table(end,i+1)     ;
end