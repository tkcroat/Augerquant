'''
AESquantparams needed for Auger quant on smooth-differentiated peaks 
This file contains AES peak centers and positions of desired background regions for typical Auger semi-quantitative analyses.
Input ROIs are cut and pasted from the excel spreadsheet named 
EDXS_quant_parameters.xls.  
Code:
'''
Si1AESpeak =94
SAESpeak =154
ClAESpeak =186
CAESpeak =276
CaAESpeak =296
OAESpeak =513
Fe1AESpeak =600
Fe2AESpeak =654
FeAESpeak =707
MgAESpeak =1185
AlAESpeak =1390
SiAESpeak =1610
NAESpeak =387
TiAESpeak =389
Ti2AESpeak =422
NaAESpeak =966
CsAESpeak =573
Cs2AESpeak =561

# Fe3 is main peak, Fe2 at 93% and Fe1 at 70%
# Corresponding positive peak positions
Si1pospeak =79
Spospeak =147
Clpospeak =179
Cpospeak =252
Capospeak =288
Opospeak =505
Fe1pospeak =582
Fe2pospeak =642
Fepospeak =697
Mgpospeak =1180
Alpospeak =1384
Sipospeak =1600
Npospeak =374
Tipospeak =379
Ti2pospeak =415
Napospeak =957 # guesstimate
Cspospeak =566
Cs2pospeak =552

# Search width for finding Auger peaks (typically varies by element)
Si1searchwidth =6
Ssearchwidth =6
Clsearchwidth =6
Csearchwidth =6
Casearchwidth =6
Osearchwidth =6
Fe1searchwidth =6
Fe2searchwidth =6
Fesearchwidth =6
Mgsearchwidth =6
Alsearchwidth =6
Sisearchwidth =8
Nsearchwidth =6
Tisearchwidth =6
Ti2searchwidth =6
Nasearchwidth =6 
Cssearchwidth =6
Cs2searchwidth =6

#Auger background checks for S7D7
b121AESpeak =121	
b121pospeak =111	
b121searchwidth =9
b205AESpeak =205	
b205pospeak =195	
b205searchwidth =9
b320AESpeak =320	
b320pospeak =310	
b320searchwidth =9
b345AESpeak =345
b345pospeak =335	
b345searchwidth =9
b450AESpeak =450	
b450pospeak =440	
b450searchwidth =9
b540AESpeak =540	
b540pospeak =530	
b540searchwidth =9
b565AESpeak =565	
b565pospeak =555	
b565searchwidth =9
b760AESpeak =760	
b760pospeak =750	
b760searchwidth =9
b1115AESpeak =1115	
b1115pospeak =1105	
b1115searchwidth =9
b1220AESpeak =1215	
b1220pospeak =1205	
b1220searchwidth =9
b1305AESpeak =1305	
b1305pospeak =1295	
b1305searchwidth =9
b1415AESpeak =1415	
b1415pospeak =1405	
b1415searchwidth =9
b1550AESpeak =1550	
b1550pospeak =1540	
b1550searchwidth =9
b1660AESpeak =1660	
b1660pospeak =1650	
b1660searchwidth =9

# compare peak-to-peak width with usual value
Si1peakwidth =15
Speakwidth =6.6
Clpeakwidth =7
Cpeakwidth =24.1
Capeakwidth =7
Opeakwidth =7.6
Fe1peakwidth =18.3
Fe2peakwidth =12.6
Fepeakwidth =9.9
Mgpeakwidth =6.8
Alpeakwidth =7
Sipeakwidth =10
Npeakwidth =13
Tipeakwidth =10
Ti2peakwidth =7
Napeakwidth =9 # guesstimate
Cspeakwidth =7
Cs2peakwidth =9

