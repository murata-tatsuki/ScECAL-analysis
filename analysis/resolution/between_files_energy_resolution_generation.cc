#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1S.h"
#include "TH2D.h"
#include "TH2S.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include "TLegend.h"
#include "TCanvas.h"
#include "weightedSSA.cpp"
#include "showerRaddius.cpp"
#include "langaus.C"        // langaus(chargeH, &fitFunc, &peakP, &peakPError);
// #include "../RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/murata_t/scecal/ScECAL_CR/analyseCode/RawtoRoot/include/EBUdecode.h"
// #include "/megraid01/users/data_beamtest_SPS2022/analysis/ECAL_Analysis/include/EBUdecode.h"

using namespace std;


// エネルギー分布を求めるマクロ
// calibを使う


const unordered_set<int> bad_channels = {
  // bad channel の CellID をここに列挙
  1, 5, 13, 14, 33, 10004, 10022, 10030, 10031, 10032, 10033, 10034, 10035, 20029, 20033, 20034, 20035, 30003, 30005, 30009, 30013, 30015, 30023, 30025, 30027, 30028, 30030, 30031, 30032, 30033, 30034, 30035, 40001, 40002, 40003, 40005, 40007, 40009, 40011, 40013, 40015, 40017, 40024, 40025, 40026, 40027, 40028, 40029, 40030, 40032, 40034, 40035, 50003, 50004, 50009, 110010, 110022, 110025, 110030, 110031, 110032, 110033, 110034, 110035, 120029, 120030, 120032, 120035, 130024, 130028, 130030, 130032, 130033, 130034, 130035, 140009, 140015, 140034, 150003, 200001, 200007, 200034, 210024, 210030, 210031, 210032, 210033, 210034, 210035, 220030, 220033, 220034, 220035, 230001, 230012, 230021, 230024, 230025, 230027, 230028, 230030, 230031, 230032, 230033, 230034, 230035, 240018, 240026, 240027, 240029, 240030, 240034, 240035, 250000, 250009, 300008, 310010, 310012, 310028, 310032, 310034, 310035, 330001, 330027, 330029, 330031, 330033, 330034, 330035, 340022, 340032, 340034, 340035, 350005, 350009, 400013, 400014, 400016, 400018, 400020, 400022, 400024, 400025, 400028, 410004, 410005, 410006, 410007, 410009, 410010, 410011, 410012, 410017, 410019, 410020, 410021, 410022, 410023, 410025, 410026, 410027, 410028, 410029, 410030, 410031, 410032, 410033, 410034, 410035, 420000, 420002, 420014, 420020, 420025, 420027, 420033, 420034, 420035, 430017, 430021, 430023, 430024, 430025, 430026, 430027, 430028, 430029, 430030, 430031, 430032, 430033, 430034, 430035, 440000, 440001, 440003, 440005, 440006, 440012, 440016, 440019, 440021, 440023, 440024, 440025, 440026, 440027, 440028, 440029, 440030, 440031, 440032, 440033, 440034, 440035, 450001, 450002, 450003, 450008, 450015, 450017, 450018, 500001, 500002, 500003, 500022, 510011, 510013, 510017, 510023, 510025, 510026, 510027, 510028, 510029, 510030, 510031, 510032, 510033, 510034, 510035, 520003, 520028, 520030, 520033, 520034, 520035, 530001, 530003, 530007, 530009, 530011, 530012, 530013, 530017, 530020, 530021, 530022, 530023, 530024, 530025, 530026, 530027, 530028, 530029, 530030, 530031, 530032, 530033, 530034, 530035, 540011, 540032, 540033, 540034, 540035, 550009, 550010, 600000, 600016, 600025, 600027, 600029, 610021, 610022, 610025, 610029, 610030, 610031, 610032, 610033, 610034, 610035, 620000, 620033, 620034, 620035, 630001, 630003, 630005, 630006, 630007, 630009, 630011, 630013, 630014, 630015, 630016, 630017, 630018, 630019, 630020, 630021, 630022, 630023, 630024, 630025, 630026, 630027, 630028, 630029, 630030, 630031, 630032, 630033, 630034, 630035, 640000, 640001, 640002, 640003, 640005, 640007, 640008, 640010, 640011, 640012, 640013, 640015, 640017, 640018, 640019, 640021, 640022, 640023, 640024, 640025, 640026, 640027, 640028, 640029, 640030, 640032, 640033, 640034, 640035, 650007, 650008, 650009, 650010, 650017, 700000, 700001, 700035, 710001, 710003, 710007, 710009, 710013, 710017, 710019, 710021, 710022, 710023, 710025, 710027, 710028, 710031, 710032, 710034, 710035, 720028, 720029, 720030, 720031, 720032, 720033, 720034, 720035, 730003, 730005, 730006, 730007, 730008, 730009, 730011, 730012, 730014, 730015, 730016, 730017, 730019, 730021, 730022, 730023, 730024, 730025, 730026, 730027, 730029, 730030, 730031, 730032, 730033, 730034, 730035, 740001, 740002, 740007, 740008, 740012, 740027, 740034, 740035, 800001, 800002, 800003, 800004, 800005, 800006, 800016, 800020, 800022, 800030, 800032, 810001, 810003, 810007, 810009, 810010, 810011, 810012, 810013, 810014, 810015, 810017, 810021, 810023, 810025, 810026, 810027, 810028, 810029, 810030, 810031, 810032, 810033, 810034, 810035, 820016, 820017, 820033, 820034, 820035, 830000, 830002, 830003, 830004, 830005, 830007, 830008, 830009, 830011, 830012, 830013, 830015, 830016, 830017, 830018, 830019, 830021, 830023, 830024, 830025, 830027, 830028, 830029, 830030, 830031, 830032, 830033, 830034, 830035, 840000, 840001, 840002, 840003, 840004, 840005, 840007, 840008, 840010, 840012, 840014, 840019, 840020, 840023, 840024, 840025, 840026, 840027, 840028, 840029, 840030, 840031, 840032, 840033, 840034, 840035, 850001, 850002, 850003, 850008, 850009, 850019, 850021, 850022, 850024, 850025, 850027, 850029, 900002, 900015, 910005, 910007, 910021, 910029, 910033, 910034, 910035, 920015, 920016, 920017, 920030, 920033, 920035, 930009, 930015, 930016, 930017, 930018, 930019, 930023, 930024, 930025, 930026, 930027, 930028, 930029, 930030, 930031, 930032, 930033, 930034, 930035, 940001, 940003, 940005, 940029, 940030, 940031, 940033, 940034, 940035, 1000002, 1000004, 1000015, 1000018, 1000020, 1000024, 1000025, 1000026, 1000028, 1000029, 1000030, 1000032, 1000035, 1010001, 1010002, 1010003, 1010005, 1010007, 1010008, 1010012, 1010014, 1010015, 1010017, 1010018, 1010019, 1010020, 1010021, 1010022, 1010023, 1010024, 1010025, 1010027, 1010028, 1010029, 1010030, 1010031, 1010032, 1010033, 1010034, 1010035, 1020002, 1020015, 1020024, 1020032, 1020033, 1020034, 1020035, 1030001, 1030003, 1030004, 1030005, 1030006, 1030007, 1030008, 1030009, 1030010, 1030011, 1030012, 1030013, 1030015, 1030016, 1030017, 1030018, 1030019, 1030020, 1030021, 1030023, 1030024, 1030025, 1030026, 1030027, 1030028, 1030029, 1030030, 1030031, 1030032, 1030033, 1030034, 1030035, 1040000, 1040001, 1040002, 1040003, 1040004, 1040005, 1040006, 1040007, 1040009, 1040010, 1040011, 1040012, 1040013, 1040014, 1040015, 1040016, 1040018, 1040019, 1040021, 1040022, 1040023, 1040024, 1040025, 1040026, 1040027, 1040028, 1040029, 1040030, 1040032, 1040033, 1040034, 1040035, 1050005, 1050009, 1050013, 1100004, 1100005, 1100007, 1100011, 1100024, 1100032, 1110006, 1110014, 1110015, 1110022, 1110024, 1110025, 1110030, 1110031, 1110032, 1110033, 1110034, 1110035, 1120001, 1120015, 1120017, 1120029, 1130003, 1130004, 1130005, 1130007, 1130009, 1130011, 1130015, 1130021, 1130023, 1130025, 1130027, 1130028, 1130029, 1130030, 1130031, 1130032, 1130033, 1130034, 1130035, 1140000, 1140004, 1140007, 1140010, 1140022, 1140025, 1140030, 1140031, 1140032, 1140033, 1140034, 1140035, 1150010, 1150011, 1200003, 1200015, 1200016, 1200018, 1200026, 1210003, 1210007, 1210009, 1210013, 1210015, 1210017, 1210018, 1210019, 1210020, 1210021, 1210022, 1210023, 1210025, 1210027, 1210029, 1210030, 1210031, 1210032, 1210033, 1210034, 1210035, 1220002, 1220024, 1220026, 1220033, 1220034, 1220035, 1230001, 1230003, 1230005, 1230007, 1230015, 1230017, 1230019, 1230021, 1230022, 1230023, 1230025, 1230026, 1230027, 1230029, 1230030, 1230031, 1230032, 1230033, 1230034, 1230035, 1240000, 1240001, 1240002, 1240003, 1240005, 1240006, 1240008, 1240009, 1240010, 1240011, 1240013, 1240014, 1240015, 1240016, 1240017, 1240019, 1240021, 1240022, 1240023, 1240024, 1240025, 1240026, 1240027, 1240028, 1240029, 1240030, 1240032, 1240033, 1240034, 1240035, 1250001, 1250002, 1250006, 1250010, 1250013, 1300001, 1300002, 1300003, 1300004, 1300006, 1300007, 1300009, 1300012, 1300020, 1310004, 1310009, 1310017, 1310019, 1310025, 1310026, 1310027, 1310028, 1310030, 1310031, 1310033, 1310034, 1310035, 1320027, 1320028, 1320029, 1320030, 1320031, 1320032, 1320033, 1320034, 1320035, 1330011, 1330015, 1330017, 1330018, 1330019, 1330020, 1330021, 1330022, 1330023, 1330024, 1330025, 1330026, 1330027, 1330028, 1330029, 1330031, 1330032, 1330033, 1330034, 1330035, 1340007, 1340019, 1340023, 1340028, 1340030, 1340032, 1340033, 1340034, 1340035, 1350006, 1350008, 1350009, 1400005, 1400016, 1400018, 1400020, 1400022, 1400024, 1400026, 1400028, 1400030, 1400032, 1410001, 1410002, 1410004, 1410005, 1410006, 1410007, 1410009, 1410010, 1410013, 1410015, 1410018, 1410019, 1410020, 1410021, 1410022, 1410023, 1410025, 1410027, 1410028, 1410029, 1410030, 1410031, 1410032, 1410033, 1410034, 1410035, 1420000, 1420001, 1420002, 1420015, 1420016, 1420017, 1420031, 1420033, 1420034, 1420035, 1430000, 1430001, 1430003, 1430005, 1430008, 1430010, 1430011, 1430013, 1430015, 1430018, 1430019, 1430021, 1430022, 1430023, 1430025, 1430027, 1430028, 1430029, 1430030, 1430031, 1430032, 1430033, 1430034, 1430035, 1440001, 1440002, 1440003, 1440004, 1440005, 1440006, 1440007, 1440008, 1440009, 1440010, 1440011, 1440012, 1440013, 1440014, 1440015, 1440016, 1440017, 1440018, 1440019, 1440020, 1440021, 1440022, 1440023, 1440024, 1440025, 1440026, 1440027, 1440028, 1440029, 1440030, 1440031, 1440032, 1440033, 1440034, 1440035, 1450000, 1450001, 1450002, 1450004, 1450006, 1450007, 1450022, 1450024, 1450026, 1450027, 1450029, 1500003, 1500006, 1510001, 1510017, 1510029, 1510032, 1510033, 1510034, 1520031, 1520034, 1530000, 1530005, 1530021, 1530022, 1530023, 1530025, 1530027, 1530028, 1530029, 1530030, 1530031, 1530033, 1530034, 1530035, 1540001, 1540011, 1540027, 1540028, 1540032, 1540034, 1540035, 1550009, 1550025, 1600000, 1600012, 1600020, 1600026, 1600030, 1600031, 1600032, 1610000, 1610005, 1610009, 1610011, 1610013, 1610015, 1610017, 1610019, 1610021, 1610022, 1610025, 1610027, 1610029, 1610030, 1610031, 1610032, 1610033, 1610034, 1610035, 1620009, 1620012, 1620013, 1620014, 1620015, 1620016, 1620018, 1620019, 1620020, 1620021, 1620022, 1620024, 1620027, 1620028, 1620030, 1620033, 1620034, 1620035, 1630005, 1630007, 1630009, 1630011, 1630013, 1630015, 1630017, 1630019, 1630023, 1630024, 1630025, 1630026, 1630027, 1630029, 1630030, 1630031, 1630032, 1630033, 1630034, 1630035, 1640002, 1640005, 1640007, 1640009, 1640011, 1640013, 1640015, 1640017, 1640019, 1640021, 1640022, 1640023, 1640026, 1640027, 1640028, 1640029, 1640033, 1640034, 1650003, 1650009, 1650010, 1650011, 1650022, 1650026, 1650027, 1650028, 1650029, 1700001, 1700003, 1700006, 1710004, 1710005, 1710009, 1710013, 1710014, 1710015, 1710019, 1710022, 1710023, 1710025, 1710026, 1710027, 1710029, 1710031, 1710032, 1710033, 1710034, 1710035, 1720029, 1720030, 1720032, 1720033, 1720034, 1720035, 1730001, 1730005, 1730007, 1730009, 1730011, 1730013, 1730016, 1730017, 1730018, 1730020, 1730027, 1730031, 1730032, 1730033, 1730035, 1740000, 1740001, 1740006, 1740026, 1740034, 1750009, 1750010, 1750021, 1800016, 1800020, 1800022, 1800028, 1800029, 1800035, 1810011, 1810015, 1810017, 1810023, 1810025, 1810028, 1810029, 1810030, 1810031, 1810032, 1810033, 1810034, 1810035, 1820016, 1820026, 1820032, 1820033, 1820034, 1820035, 1830001, 1830003, 1830004, 1830005, 1830006, 1830007, 1830009, 1830010, 1830011, 1830012, 1830013, 1830014, 1830015, 1830017, 1830019, 1830020, 1830021, 1830022, 1830023, 1830024, 1830025, 1830027, 1830028, 1830029, 1830030, 1830031, 1830032, 1830033, 1830034, 1830035, 1840001, 1840002, 1840003, 1840004, 1840005, 1840007, 1840008, 1840009, 1840010, 1840011, 1840013, 1840014, 1840015, 1840016, 1840017, 1840018, 1840020, 1840021, 1840022, 1840023, 1840024, 1840025, 1840026, 1840027, 1840028, 1840029, 1840030, 1840032, 1840033, 1840034, 1840035, 1850000, 1850001, 1850003, 1850004, 1850006, 1850007, 1850009, 1850010, 1850024, 1850025, 1900000, 1900001, 1900004, 1900005, 1900006, 1900026, 1910003, 1910007, 1910009, 1910015, 1910017, 1910019, 1910021, 1910022, 1910023, 1910025, 1910026, 1910027, 1910029, 1910030, 1910031, 1910032, 1910033, 1910034, 1910035, 1920016, 1920025, 1920028, 1920029, 1920030, 1920031, 1920032, 1920034, 1920035, 1930009, 1930013, 1930016, 1930018, 1930020, 1930022, 1930023, 1930025, 1930027, 1930028, 1930029, 1930031, 1930032, 1930033, 1930035, 1940000, 1940003, 1940017, 1940026, 1940032, 1940033, 1940034, 1940035, 1950005, 1950009, 1950010, 1950013, 2000000, 2000018, 2000020, 2000022, 2000028, 2010001, 2010002, 2010006, 2010007, 2010008, 2010009, 2010010, 2010011, 2010013, 2010014, 2010015, 2010016, 2010017, 2010018, 2010019, 2010020, 2010021, 2010022, 2010023, 2010025, 2010026, 2010027, 2010028, 2010029, 2010030, 2010031, 2010032, 2010033, 2010034, 2010035, 2020006, 2020015, 2020017, 2020028, 2020029, 2020032, 2020033, 2020034, 2020035, 2030000, 2030001, 2030005, 2030006, 2030007, 2030013, 2030015, 2030017, 2030019, 2030023, 2030024, 2030025, 2030027, 2030029, 2030030, 2030031, 2030032, 2030033, 2030034, 2030035, 2040000, 2040001, 2040002, 2040003, 2040004, 2040007, 2040008, 2040010, 2040011, 2040015, 2040016, 2040017, 2040018, 2040019, 2040020, 2040021, 2040022, 2040023, 2040024, 2040025, 2040026, 2040027, 2040028, 2040029, 2040030, 2040031, 2040032, 2040033, 2040034, 2040035, 2050002, 2050006, 2050009, 2050010, 2050011, 2050012, 2050021, 2050023, 2050029, 2100000, 2100002, 2100003, 2100005, 2100006, 2100023, 2110019, 2110025, 2110027, 2110028, 2110032, 2110034, 2110035, 2120017, 2120031, 2120032, 2120033, 2120034, 2120035, 2130002, 2130007, 2130009, 2130010, 2130011, 2130013, 2130017, 2130019, 2130020, 2130021, 2130023, 2130024, 2130025, 2130026, 2130027, 2130029, 2130030, 2130031, 2130032, 2130033, 2130034, 2130035, 2140019, 2140026, 2140027, 2140030, 2140032, 2140033, 2140034, 2140035, 2150009, 2150010, 2150021, 2200020, 2200022, 2200024, 2210017, 2210019, 2210023, 2210025, 2210027, 2210028, 2210029, 2210030, 2210031, 2210032, 2210033, 2210034, 2210035, 2220000, 2220001, 2220033, 2220034, 2220035, 2230002, 2230007, 2230008, 2230011, 2230015, 2230017, 2230020, 2230023, 2230025, 2230027, 2230028, 2230030, 2230031, 2230032, 2230033, 2230034, 2230035, 2240000, 2240001, 2240002, 2240003, 2240005, 2240006, 2240007, 2240009, 2240011, 2240012, 2240013, 2240014, 2240015, 2240017, 2240018, 2240019, 2240020, 2240021, 2240022, 2240023, 2240024, 2240025, 2240026, 2240027, 2240028, 2240029, 2240030, 2240032, 2240034, 2240035, 2250007, 2250022, 2250026, 2250028, 2300000, 2300001, 2300003, 2300004, 2300006, 2300007, 2300013, 2300033, 2300034, 2310005, 2310007, 2310011, 2310013, 2310015, 2310019, 2310021, 2310023, 2310024, 2310025, 2310026, 2310027, 2310028, 2310029, 2310030, 2310031, 2310032, 2310033, 2310034, 2310035, 2320028, 2320029, 2320031, 2320032, 2320033, 2320034, 2330000, 2330001, 2330003, 2330005, 2330009, 2330010, 2330011, 2330012, 2330015, 2330016, 2330017, 2330018, 2330021, 2330022, 2330023, 2330024, 2330025, 2330026, 2330027, 2330028, 2330029, 2330031, 2330032, 2330033, 2330035, 2340027, 2340030, 2340031, 2340032, 2340033, 2340034, 2340035, 2350008, 2350009, 2350010, 2350021, 2400016, 2400020, 2400022, 2410005, 2410007, 2410010, 2410013, 2410019, 2410023, 2410025, 2410027, 2410029, 2410030, 2410031, 2410032, 2410033, 2410034, 2410035, 2420000, 2420011, 2420015, 2420017, 2420033, 2420034, 2420035, 2430000, 2430001, 2430002, 2430003, 2430005, 2430007, 2430009, 2430011, 2430013, 2430015, 2430017, 2430019, 2430020, 2430021, 2430022, 2430023, 2430024, 2430027, 2430028, 2430029, 2430030, 2430031, 2430032, 2430033, 2430034, 2430035, 2440000, 2440001, 2440002, 2440003, 2440005, 2440007, 2440009, 2440010, 2440011, 2440012, 2440013, 2440014, 2440015, 2440016, 2440017, 2440019, 2440020, 2440021, 2440022, 2440023, 2440024, 2440025, 2440026, 2440027, 2440028, 2440029, 2440030, 2440031, 2440032, 2440033, 2440035, 2450000, 2450001, 2450013, 2450022, 2450023, 2450028, 2450029, 2500000, 2500001, 2500002, 2500003, 2500004, 2500005, 2500006, 2500007, 2500013, 2500014, 2500015, 2500022, 2500034, 2500035, 2510004, 2510010, 2510015, 2510017, 2510019, 2510020, 2510021, 2510022, 2510023, 2510024, 2510025, 2510026, 2510027, 2510028, 2510029, 2510030, 2510031, 2510032, 2510033, 2510034, 2510035, 2520032, 2520033, 2520034, 2530001, 2530003, 2530005, 2530009, 2530011, 2530015, 2530016, 2530019, 2530020, 2530021, 2530022, 2530023, 2530025, 2530026, 2530027, 2530028, 2530029, 2530030, 2530031, 2530032, 2530033, 2530034, 2530035, 2540001, 2540003, 2540008, 2540010, 2540013, 2540014, 2540015, 2540024, 2540025, 2540027, 2540030, 2540031, 2540033, 2540034, 2540035, 2550008, 2550009, 2550015, 2550021, 2600001, 2600003, 2600005, 2600012, 2600018, 2600020, 2600022, 2600024, 2600026, 2600030, 2600032, 2610001, 2610003, 2610004, 2610005, 2610007, 2610011, 2610013, 2610014, 2610015, 2610017, 2610019, 2610020, 2610021, 2610022, 2610023, 2610024, 2610025, 2610027, 2610029, 2610030, 2610031, 2610032, 2610033, 2610034, 2610035, 2620001, 2620015, 2620017, 2620033, 2620034, 2620035, 2630001, 2630003, 2630005, 2630006, 2630007, 2630009, 2630011, 2630012, 2630013, 2630015, 2630017, 2630019, 2630021, 2630022, 2630023, 2630024, 2630025, 2630027, 2630028, 2630029, 2630030, 2630031, 2630032, 2630033, 2630034, 2630035, 2640000, 2640001, 2640002, 2640003, 2640005, 2640006, 2640007, 2640009, 2640010, 2640011, 2640013, 2640015, 2640016, 2640017, 2640019, 2640021, 2640022, 2640023, 2640024, 2640025, 2640026, 2640027, 2640028, 2640029, 2640030, 2640031, 2640032, 2640033, 2640034, 2640035, 2650006, 2650021, 2650023, 2650024, 2650027, 2650028, 2650029, 2700000, 2700002, 2700004, 2700009, 2700012, 2700030, 2700032, 2710000, 2710001, 2710002, 2710003, 2710007, 2710009, 2710013, 2710015, 2710016, 2710019, 2710021, 2710025, 2710026, 2710029, 2710030, 2710032, 2710033, 2710034, 2710035, 2720027, 2720032, 2720033, 2720035, 2730001, 2730005, 2730006, 2730007, 2730009, 2730011, 2730012, 2730013, 2730015, 2730016, 2730019, 2730022, 2730023, 2730025, 2730027, 2730028, 2730029, 2730031, 2730032, 2730033, 2730034, 2730035, 2740009, 2740013, 2740016, 2740030, 2740032, 2740033, 2740034, 2750009, 2750010, 2750011, 2750023, 2800005, 2800018, 2800021, 2800022, 2800032, 2810006, 2810013, 2810018, 2810025, 2810027, 2810028, 2810029, 2810030, 2810031, 2810032, 2810033, 2810034, 2810035, 2820006, 2820016, 2820017, 2820026, 2820033, 2820034, 2820035, 2830001, 2830003, 2830005, 2830007, 2830009, 2830010, 2830011, 2830013, 2830014, 2830015, 2830016, 2830017, 2830019, 2830021, 2830023, 2830025, 2830027, 2830028, 2830029, 2830030, 2830031, 2830032, 2830033, 2830034, 2830035, 2840000, 2840001, 2840002, 2840003, 2840005, 2840013, 2840015, 2840016, 2840019, 2840021, 2840022, 2840023, 2840024, 2840025, 2840026, 2840027, 2840028, 2840029, 2840030, 2840032, 2840033, 2840034, 2840035, 2850009, 2850019, 2900000, 2900001, 2900007, 2910000, 2910007, 2910019, 2910021, 2910022, 2910023, 2910024, 2910027, 2910029, 2910031, 2910032, 2910033, 2910034, 2910035, 2920023, 2920032, 2920033, 2920034, 2920035, 2930016, 2930019, 2930022, 2930023, 2930024, 2930025, 2930026, 2930027, 2930028, 2930029, 2930030, 2930031, 2930032, 2930033, 2930034, 2930035, 2940015, 2940017, 2940022, 2940027, 2940032, 2940034, 2940035, 2950009, 2950012, 3000000, 3000001, 3000002, 3000003, 3000004, 3000005, 3000006, 3000007, 3000008, 3000009, 3000010, 3000011, 3000012, 3000013, 3000014, 3000015, 3000016, 3000017, 3000018, 3000019, 3000020, 3000021, 3000022, 3000023, 3000024, 3000025, 3000026, 3000027, 3000028, 3000029, 3000030, 3000031, 3000032, 3000033, 3000034, 3000035, 3010000, 3010001, 3010002, 3010003, 3010004, 3010005, 3010006, 3010007, 3010008, 3010009, 3010010, 3010011, 3010012, 3010013, 3010014, 3010015, 3010016, 3010017, 3010018, 3010019, 3010020, 3010021, 3010022, 3010023, 3010024, 3010025, 3010026, 3010027, 3010028, 3010029, 3010030, 3010031, 3010032, 3010033, 3010034, 3010035, 3020000, 3020001, 3020002, 3020003, 3020004, 3020005, 3020006, 3020007, 3020008, 3020009, 3020010, 3020011, 3020012, 3020013, 3020014, 3020015, 3020016, 3020017, 3020018, 3020019, 3020020, 3020021, 3020022, 3020023, 3020024, 3020025, 3020026, 3020027, 3020028, 3020029, 3020030, 3020031, 3020032, 3020033, 3020034, 3020035, 3030000, 3030001, 3030002, 3030003, 3030004, 3030005, 3030006, 3030007, 3030008, 3030009, 3030010, 3030011, 3030012, 3030013, 3030014, 3030015, 3030016, 3030017, 3030018, 3030019, 3030020, 3030021, 3030022, 3030023, 3030024, 3030025, 3030026, 3030027, 3030028, 3030029, 3030030, 3030031, 3030032, 3030033, 3030034, 3030035, 3040000, 3040001, 3040002, 3040003, 3040004, 3040005, 3040006, 3040007, 3040008, 3040009, 3040010, 3040011, 3040012, 3040013, 3040014, 3040015, 3040016, 3040017, 3040018, 3040019, 3040020, 3040021, 3040022, 3040023, 3040024, 3040025, 3040026, 3040027, 3040028, 3040029, 3040030, 3040031, 3040032, 3040033, 3040034, 3040035, 3050000, 3050001, 3050002, 3050003, 3050004, 3050005, 3050006, 3050007, 3050008, 3050009, 3050010, 3050011, 3050012, 3050013, 3050014, 3050015, 3050016, 3050017, 3050018, 3050019, 3050020, 3050021, 3050022, 3050023, 3050024, 3050025, 3050026, 3050027, 3050028, 3050029, 3100000, 3100001, 3100002, 3100003, 3100004, 3100005, 3100006, 3100007, 3100008, 3100009, 3100010, 3100011, 3100012, 3100013, 3100014, 3100015, 3100016, 3100017, 3100018, 3100019, 3100020, 3100021, 3100022, 3100023, 3100024, 3100025, 3100026, 3100027, 3100028, 3100029, 3100030, 3100031, 3100032, 3100033, 3100034, 3100035, 3110000, 3110001, 3110002, 3110003, 3110004, 3110005, 3110006, 3110007, 3110008, 3110009, 3110010, 3110011, 3110012, 3110013, 3110014, 3110015, 3110016, 3110017, 3110018, 3110019, 3110020, 3110021, 3110022, 3110023, 3110024, 3110025, 3110026, 3110027, 3110028, 3110029, 3110030, 3110031, 3110032, 3110033, 3110034, 3110035, 3120000, 3120001, 3120002, 3120003, 3120004, 3120005, 3120006, 3120007, 3120008, 3120009, 3120010, 3120011, 3120012, 3120013, 3120014, 3120015, 3120016, 3120017, 3120018, 3120019, 3120020, 3120021, 3120022, 3120023, 3120024, 3120025, 3120026, 3120027, 3120028, 3120029, 3120030, 3120031, 3120032, 3120033, 3120034, 3120035, 3130000, 3130001, 3130002, 3130003, 3130004, 3130005, 3130006, 3130007, 3130008, 3130009, 3130010, 3130011, 3130012, 3130013, 3130014, 3130015, 3130016, 3130017, 3130018, 3130019, 3130020, 3130021, 3130022, 3130023, 3130024, 3130025, 3130026, 3130027, 3130028, 3130029, 3130030, 3130031, 3130032, 3130033, 3130034, 3130035, 3140000, 3140001, 3140002, 3140003, 3140004, 3140005, 3140006, 3140007, 3140008, 3140009, 3140010, 3140011, 3140012, 3140013, 3140014, 3140015, 3140016, 3140017, 3140018, 3140019, 3140020, 3140021, 3140022, 3140023, 3140024, 3140025, 3140026, 3140027, 3140028, 3140029, 3140030, 3140031, 3140032, 3140033, 3140034, 3140035, 3150000, 3150001, 3150002, 3150003, 3150004, 3150005, 3150006, 3150007, 3150008, 3150009, 3150010, 3150011, 3150012, 3150013, 3150014, 3150015, 3150016, 3150017, 3150018, 3150019, 3150020, 3150021, 3150022, 3150023, 3150024, 3150025, 3150026, 3150027, 3150028, 3150029
};

bool cog_position_check(double x, double y, double range=200){
  return (abs(x)<range && abs(y)<range);
}

#include <vector>
#include <cmath>
#include <iostream>

/*
double CalculateShowerRadiusRMS(const vector<double>& x, const vector<double>& y, const vector<double>& energy){
  // データのサイズチェック
  if (x.empty() || x.size() != y.size() || x.size() != energy.size()) {
    std::cerr << "Warning: Input vector sizes mismatch or empty." << std::endl;
    return -1.0; 
  }

  double E_tot = 0.0;
  double sum_X = 0.0;
  double sum_Y = 0.0;
  size_t n_hits = x.size();

  // 1. 総エネルギーとエネルギー重心(Center of Gravity)の計算
  for (size_t i = 0; i < n_hits; ++i) {
    if (energy[i] <= 0) continue; // ノイズ等による負のエネルギーや0を排除
    E_tot += energy[i];
    sum_X += x[i] * energy[i];
    sum_Y += y[i] * energy[i];
  }

  // 総エネルギーがゼロ（ヒットなし）の場合はエラーを返す
  if (E_tot <= 0.0) return 0.0;

  double X_CoG = sum_X / E_tot;
  double Y_CoG = sum_Y / E_tot;

  // 2. 各軸の分散(Variance)の計算
  double sum_var_X = 0.0;
  double sum_var_Y = 0.0;

  for (size_t i = 0; i < n_hits; ++i) {
    if (energy[i] <= 0) continue;
    sum_var_X += energy[i] * std::pow(x[i] - X_CoG, 2);
    sum_var_Y += energy[i] * std::pow(y[i] - Y_CoG, 2);
  }

  double var_x = sum_var_X / E_tot;
  double var_y = sum_var_Y / E_tot;

  // 3. RMS半径 (R_RMS = sqrt(sigma_x^2 + sigma_y^2)) の計算
  double R_num_RMS = std::sqrt(var_x + var_y);

  // 4. 2次元ガウシアンを仮定してR90へ換算 (R90 ≒ 1.518 * R_RMS)
  double R90 = 1.518 * R_num_RMS;

  return R90;
}
*/

bool thresholdCheck(double Edep, double threshold, double thresholdSigma, double MIP){
  double E_th_i = threshold * (0.308 / MIP); // チャンネルごとの閾値
  double sigma_i = thresholdSigma * (0.308 / MIP); // ノイズ幅
  
  // ヒット保持確率を計算 (TMath::Erf を使用)
  double prob = 0.5 * (1.0 + TMath::Erf((Edep - E_th_i) / (sqrt(2.0) * sigma_i)));

  // 一様乱数を用いて確率的にヒットを残す
  return gRandom->Uniform(0,1) < prob;

  // return E_th_i < Edep;
}

int main(int argc, char* argv[])
{ 
    if(argc < 3){                                                     //エラー処理
        cout << "usage: ./between_files  output.root Nenegry E1...En E1Nfiles...EnNfiles input1.root input2.root  ...  cog_range b_only_best figure_path " << endl;
        return 1;
    }
  gROOT->SetBatch(kTRUE);

    //cout << argv[1] << endl;
  cout << "=====>  " << argv[1] << endl;

  int Nenegry = atoi(argv[2]);
  double energy[Nenegry];
  int energy_files[Nenegry];
  for(int i=0;i<Nenegry;i++){
    energy[i] = atoi(argv[3+i]);
    if(energy[i]<1) energy[i] = 0.5;
    energy_files[i] = atoi(argv[3+Nenegry+i]);
  }
  
  
  int rawfilenum = argc - 2 - 2*Nenegry - 3;
  int irawfilenum = 3 + 2*Nenegry;
  int cog_range = atoi(argv[argc-3]);
  int bool_only_best = atoi(argv[argc-2]);


  // double cycleID, triggerID;
  // vector<int> *cellIDs = nullptr;
  // vector<int> *BCIDs = nullptr;
  // vector<int> *hitTags = nullptr;
  // vector<int> *gainTags = nullptr;
  // vector<double> *charges = nullptr;
  // vector<double> *times = nullptr;
  // vector<double> *temp = nullptr;
  // vector<double> *posX = nullptr;
  // vector<double> *posY = nullptr;
  // vector<double> *posZ = nullptr;

  // double cycleID_ssa, triggerID_ssa;
  // vector<int> *cellIDs_ssa = nullptr;
  // vector<int> *BCIDs_ssa = nullptr;
  // vector<int> *hitTags_ssa = nullptr;
  // vector<int> *gainTags_ssa = nullptr;
  // vector<double> *charges_ssa = nullptr;
  // vector<double> *times_ssa = nullptr;
  // vector<double> *temp_ssa = nullptr;
  // vector<double> *posX_ssa = nullptr;
  // vector<double> *posY_ssa = nullptr;
  // vector<double> *posZ_ssa = nullptr;
  // vector<int> *ssaTag_ssa = nullptr;

  double EnergyDep;
  vector<int> *cellID = nullptr;
  vector<double> *Hit_Energy = nullptr;
  vector<double> *Hit_X = nullptr;
  vector<double> *Hit_Y = nullptr;
  vector<double> *Hit_Z = nullptr;
  


  const int layerNu = 32;
  const int chipNu = 6;
  const int channelNu = 36;
  const int rowNu = 5;
  const int colNu = 42;
  const int doublelayerNu = 2;

  const int range_maximum = 12000;
  const int range_maximum_PS = 500;
  const double binWidth_PS = 0.1;



  TFile fileout(argv[1],"RECREATE");
  fileout.mkdir("all__channel");
  fileout.mkdir("best_channel");
  // gain_histo->Write();



  // dataを入れるもの
  TH1F* energy_deposition_bestChannels[Nenegry];
  TH1F* energy_deposition[Nenegry];
  TH2F* hitmap_bestChannels[Nenegry][layerNu];
  TH2F* hitmap[Nenegry][layerNu];
  TH1F* number_of_hits_bestChannels[Nenegry];
  TH1F* number_of_hits[Nenegry];
  TH2F* hit_vs_e[Nenegry];
  TH2F* hit_vs_e_bestChannels[Nenegry];
  TH2F* layer_vs_edep[Nenegry];
  TH2F* layer_vs_edep_bestChannels[Nenegry];
  TH2F* cog[Nenegry][layerNu];
  TH1F* shower_radiusRMS[Nenegry];
  TH1F* shower_radius90[Nenegry];
  TH2F* layer_vs_rRMS[Nenegry];
  TH2F* layer_vs_r90[Nenegry];
  TH1F* edep_1hit[Nenegry][2];
  // TH2F* edep_vs_rRMS[Nenegry];
  // TH2F* edep_vs_r90[Nenegry];
  for (int i = 0; i < Nenegry; ++i) {
    string histo_title = energy[i]<1 ? "0.5 GeV energy deposition" : Form("%g GeV energy deposition",energy[i]);
    int nbin = energy[i]>5 ? range_maximum : range_maximum_PS/binWidth_PS;
    int max_range = energy[i]>5 ? range_maximum : range_maximum_PS;
    energy_deposition[i] = new TH1F(Form("energy_deposition_%d",i),Form("%s",histo_title.c_str()),nbin,0,max_range);
    energy_deposition[i]->SetXTitle("Reconstructed Energy [MeV]");
    energy_deposition_bestChannels[i] = new TH1F(Form("energy_deposition_bestChannels_%d",i),Form("%s",histo_title.c_str()),nbin,0,max_range);
    energy_deposition_bestChannels[i]->SetXTitle("Detected Energy [MeV]");
    
    number_of_hits[i] = new TH1F(Form("nhits_%d",i),Form("%g GeV number of hits",energy[i]),3000,0,3000);
    number_of_hits[i]->SetXTitle("number of hits");
    number_of_hits_bestChannels[i] = new TH1F(Form("nhits_bestChannels_%d",i),Form("%g GeV number of hits",energy[i]),3000,0,3000);
    number_of_hits_bestChannels[i]->SetXTitle("number of hits");

    hit_vs_e[i] = new TH2F(Form("hit_vs_e_%d",i),Form("%g GeV nhits vs. edep",energy[i]), nbin,0,max_range, 3000,0,3000);
    hit_vs_e[i]->SetXTitle("Detected Energy [MeV]");
    hit_vs_e[i]->SetYTitle("number of hits");
    hit_vs_e_bestChannels[i] = new TH2F(Form("hit_vs_e_bestChannels_%d",i),Form("%g GeV nhits vs. edep",energy[i]), nbin,0,max_range, 3000,0,3000);
    hit_vs_e_bestChannels[i]->SetXTitle("Detected Energy [MeV]");
    hit_vs_e_bestChannels[i]->SetYTitle("number of hits");

    double edep_max = energy[i]>5 ? 1000 : 100;
    layer_vs_edep[i] = new TH2F(Form("layer_vs_edep_%d",i),Form("%g GeV layer vs. edep",energy[i]), layerNu,0,layerNu, edep_max,0,edep_max);
    layer_vs_edep[i]->SetXTitle("layer");
    layer_vs_edep[i]->SetYTitle("Detected Energy [MeV]");
    layer_vs_edep_bestChannels[i] = new TH2F(Form("layer_vs_edep_bestChannels_%d",i),Form("%g GeV layer vs. edep",energy[i]), layerNu,0,layerNu, edep_max,0,edep_max);
    layer_vs_edep_bestChannels[i]->SetXTitle("layer");
    layer_vs_edep_bestChannels[i]->SetYTitle("Detected Energy [MeV]");
    
    for (int ilayer = 0; ilayer < layerNu; ++ilayer) {
      hitmap_bestChannels[i][ilayer] = new TH2F(Form("hitmap_bestChannels_%dGeV_layer%d",i,ilayer),Form("hitmap %g GeV layer %d",energy[i],ilayer),226,-113,113,226,-113,113);
      hitmap[i][ilayer] = new TH2F(Form("hitmap_%dGeV_layer%d",i,ilayer),Form("hitmap %g GeV layer %d",energy[i],ilayer),226,-113,113,226,-113,113);
      cog[i][ilayer] = new TH2F(Form("cog_%dGeV_layer_%d",i,ilayer),Form("center of gravity %g GeV layer %d",energy[i],ilayer),226,-113,113,226,-113,113);
    }

    shower_radiusRMS[i] = new TH1F(Form("shower_radiusRMS_%d",i),Form("%g GeV shower radius (RMS)",energy[i]),2000,0,200);
    shower_radiusRMS[i]->SetXTitle("radius [mm]");
    shower_radius90[i] = new TH1F(Form("shower_radius90_%d",i),Form("%g GeV shower radius R90",energy[i]),2000,0,200);
    shower_radius90[i]->SetXTitle("radius [mm]");
    layer_vs_r90[i] = new TH2F(Form("layer_vs_r90_%d",i),Form("%g GeV layer vs. r90",energy[i]), layerNu,0,layerNu, 2000,0,200);
    layer_vs_r90[i]->SetXTitle("layer");
    layer_vs_r90[i]->SetYTitle("radius [mm]");
    layer_vs_rRMS[i] = new TH2F(Form("layer_vs_rRMS_%d",i),Form("%g GeV layer vs. rRMS",energy[i]), layerNu,0,layerNu, 2000,0,200);
    layer_vs_rRMS[i]->SetXTitle("layer");
    layer_vs_rRMS[i]->SetYTitle("radius [mm]");
    // edep_vs_rRMS[i] = new TH2F(Form("edep_vs_rRMS_%d",i),Form("%g GeV edep vs. rRMS",energy[i]), nbin,0,max_range, 2000,0,200);
    // edep_vs_rRMS[i]->SetXTitle("layer");
    // edep_vs_rRMS[i]->SetYTitle("radius [mm]");
    // edep_vs_r90[i] = new TH2F(Form("edep_vs_r90_%d",i),Form("%g GeV edep vs. rRMS",energy[i]), nbin,0,max_range, 2000,0,200);
    // edep_vs_r90[i]->SetXTitle("layer");
    // edep_vs_r90[i]->SetYTitle("radius [mm]");

    edep_1hit[i][1] = new TH1F(Form("edep_1hit_15um_%d",i),Form("%g GeV Edep of 1 hit (15 um)",energy[i]), 12000,0,12);
    edep_1hit[i][1]->SetXTitle("[MeV]");
    edep_1hit[i][0] = new TH1F(Form("edep_1hit_10um_%d",i),Form("%g GeV Edep of 1 hit (10 um)",energy[i]), 12000,0,12);
    edep_1hit[i][0]->SetXTitle("[MeV]");
  }
  
  

  // MIP チェック
  cout<<"=========== MIP Extract ============"<<endl;
    TFile *fileinMIP;
    fileinMIP = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/all_auto_muon_v4_trackfit.root");
    TTree *treeMIP = (TTree*) fileinMIP->Get("MIP_Fit");
    if(!treeMIP)  {cout<<"!!! GET MIP FILE FAILED !!!"<<endl;}
    int entry_max_MIP = treeMIP->GetEntries();
    Int_t _mipCellID, _NDF;
    Double_t _landauMPV, _ChiSqr;
    treeMIP->SetBranchAddress("CellID",&_mipCellID);
    treeMIP->SetBranchAddress("LandauMPV",&_landauMPV);
    treeMIP->SetBranchAddress("ChiSquare",&_ChiSqr);
    treeMIP->SetBranchAddress("NDF",&_NDF);
    double MIPPeakValue[layerNu][chipNu][channelNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          MIPPeakValue[i_layer][i_chip][i_channel] = 1;
        }
      }
    }
    double meanMIP15  = 0 ; double meanMIP10  = 0;
    double chnNoMIP15 = 0 ; double chnNoMIP10 = 0;
    for(int entry=0; entry!=treeMIP->GetEntries(); entry++){
      treeMIP->GetEntry(entry);
      int layerID = _mipCellID/1e5;
      if(layerID>=30)  continue;  
      int chipID  = (_mipCellID%100000) /1e4;
      int chanID  = _mipCellID%100;
      if(_landauMPV<0) continue;
      if (layerID>3 && layerID<28) {
        // if(_ChiSqr/(double)_NDF>2) continue;
        chnNoMIP10++;
        meanMIP10 +=_landauMPV;
      } else {
        // if(_ChiSqr/(double)_NDF>1.7) continue;
        chnNoMIP15++;
        meanMIP15 +=_landauMPV;
      }
      MIPPeakValue[layerID][chipID][chanID] = _landauMPV;
    }
    meanMIP10 /= chnNoMIP10;
    meanMIP15 /= chnNoMIP15;
    cout<<" channel count(15um:10um) : "<<chnNoMIP15<<" , "<<chnNoMIP10<<endl;
    cout<<" mean of MIPPeakValue : "<<meanMIP15<<" , "<<meanMIP10<<endl;

    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if (il>3 && il<28) {
            MIPPeakValue[il][ip][ic] = (MIPPeakValue[il][ip][ic]==1) ? meanMIP10 : MIPPeakValue[il][ip][ic];
          } else {
            MIPPeakValue[il][ip][ic] = (MIPPeakValue[il][ip][ic]==1) ? meanMIP15 : MIPPeakValue[il][ip][ic];
          }
        }
      }
    }
  fileinMIP->Close();
  fileout.cd();



  
  cout<<"=========== Threshold Extract ============"<<endl;
    TFile *ThresholdFile = new TFile("/megraid01/users/murata_t/scecal/ScECAL_BeamTest/analysis/result/threshold/threshold.root","READ");
    TTree* ThresholdTree = (TTree*)ThresholdFile->Get("threshold");
    Int_t _thresholdCellID;
    Double_t _threshold, _thresholdSigma, _thoresold_with_pedestal;
    ThresholdTree->SetBranchAddress("CellID",&_thresholdCellID);
    ThresholdTree->SetBranchAddress("threshold",&_threshold);
    ThresholdTree->SetBranchAddress("sigma",&_thresholdSigma);

    double thresholdValue[layerNu][chipNu][channelNu];
    double thresholdSigmaValue[layerNu][chipNu][channelNu];
    double thresholdValue_chip[layerNu][chipNu];
    double thresholdSigmaValue_chip[layerNu][chipNu];
    double thresholdValue_chipCount[layerNu][chipNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          thresholdValue[i_layer][i_chip][i_channel] = 1;
          thresholdSigmaValue[i_layer][i_chip][i_channel] = 0;
        }
        thresholdValue_chip[i_layer][i_chip] = 1;
        thresholdSigmaValue_chip[i_layer][i_chip] = 0;
        thresholdValue_chipCount[i_layer][i_chip] = 0;
      }
    }

    double meanThreshold15=0; double meanThreshold10=0;
    double sigmaThreshold15=0; double sigmaThreshold10=0;
    double chnNoThreshold15=0; double chnNoThreshold10=0;
    for(int entry=0; entry!=ThresholdTree->GetEntries(); entry++){
      ThresholdTree->GetEntry(entry);
      int layerID = _thresholdCellID/1e5;
      if(layerID>=30)  continue;  
      int chipID  = (_thresholdCellID%100000) /1e4;
      int chanID   = _thresholdCellID%100;
      if (layerID>3 && layerID<28) {
        // if(_ChiSqr/(double)_NDF>2) continue;
        chnNoThreshold10++;
        meanThreshold10  += _threshold;
        sigmaThreshold10 += _thresholdSigma;
        // meanThreshold10 += _thoresold_with_pedestal - thres;
        // meanThreshold10 += _thoresold_with_pedestal -_threshold + thres;
      } else {
        // if(_ChiSqr/(double)_NDF>1.7) continue;
        chnNoThreshold15++;
        meanThreshold15  += _threshold;
        sigmaThreshold15 += _thresholdSigma;
        // meanThreshold15 += _thoresold_with_pedestal - thres;
        // meanThreshold15 += _thoresold_with_pedestal -_threshold + thres;
      }
        // if (_threshold>50) {
        thresholdValue_chip[layerID][chipID] += _threshold;
        thresholdSigmaValue_chip[layerID][chipID] += _thresholdSigma;
        thresholdValue_chipCount[layerID][chipID]++;
        // }
        thresholdValue[layerID][chipID][chanID] = _threshold;
        thresholdSigmaValue[layerID][chipID][chanID] = _thresholdSigma;
        // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal - thres;
        // thresholdValue_with_pedestal[layerID][chipID][chanID] = _thoresold_with_pedestal -_threshold + thres;
    }
    meanThreshold10 /= chnNoThreshold10;
    meanThreshold15 /= chnNoThreshold15;
    sigmaThreshold10 /= chnNoThreshold10;
    sigmaThreshold15 /= chnNoThreshold15;
    cout << " channel count(15um:10um) : " << chnNoThreshold15 << " , " << chnNoThreshold10 << endl;
    cout << " mean of Threshold : " << meanThreshold15 << " , " << meanThreshold10 << endl;
    cout << " sigma of Threshold : " << sigmaThreshold15 << " , " << sigmaThreshold10 << endl;

    for (int il=0; il<30; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if (il>3 && il<28) {
            thresholdValue[il][ip][ic] = (thresholdValue[il][ip][ic]==1) ? thresholdValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdValue[il][ip][ic];
            thresholdSigmaValue[il][ip][ic] = (thresholdSigmaValue[il][ip][ic]==0) ? thresholdSigmaValue_chip[il][ip]/thresholdValue_chipCount[il][ip] : thresholdSigmaValue[il][ip][ic];
          } 
        }
      }
    }
  ThresholdFile->Close();
  fileout.cd();



  cout<<"=========== PedestalExtract ============"<<endl;
    TFile* pedFile = new TFile("/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/pedestal2023_SPS.root");
    TTree* pedTree = (TTree*)pedFile->Get("ChnLevel"); 
    if(!pedTree)  {cout<<"!!! GET PED FILE FAILED !!!"<<endl;}
    vector<int>* _pedCellID = nullptr;
    vector<double>* _pedMeanHigh = nullptr;
    vector<double>* _pedSigmaHigh = nullptr;
    vector<double>* _pedMeanLow = nullptr;
    vector<double>* _pedSigmaLow = nullptr;
    pedTree ->SetBranchAddress("CellID",&_pedCellID);
    pedTree ->SetBranchAddress("PedHighMean",&_pedMeanHigh);
    pedTree ->SetBranchAddress("PedHighSig",&_pedSigmaHigh);
    pedTree ->SetBranchAddress("PedLowMean",&_pedMeanLow);
    pedTree ->SetBranchAddress("PedLowSig",&_pedSigmaLow);

    double meanPedHigh15  = 0 ; double meanPedHigh10  = 0;
    double sigmaPedHigh15 = 0 ; double sigmaPedHigh10 = 0;
    double chnNoPedHigh15 = 0 ; double chnNoPedHigh10 = 0;
    double meanPedLow15   = 0 ; double meanPedLow10   = 0;
    double sigmaPedLow15  = 0 ; double sigmaPedLow10  = 0;
    double chnNoPedLow15  = 0 ; double chnNoPedLow10  = 0;
    double pedMeanHigh[layerNu][chipNu][channelNu];
    double pedSigmaHigh[layerNu][chipNu][channelNu];
    double pedMeanLow[layerNu][chipNu][channelNu];
    double pedSigmaLow[layerNu][chipNu][channelNu];
    TH2F *ped_map[layerNu];
    for (int i_layer = 0; i_layer < layerNu; ++i_layer) {
      for (int i_chip = 0; i_chip < chipNu; ++i_chip) {
        for (int i_channel = 0; i_channel < channelNu; ++i_channel) {
          pedMeanHigh[i_layer][i_chip][i_channel] = 1;
          pedSigmaHigh[i_layer][i_chip][i_channel] = 1;
          pedMeanLow[i_layer][i_chip][i_channel] = 1;
          pedSigmaLow[i_layer][i_chip][i_channel] = 1;
        }
      }
    }
    for(int entry=0; entry!=pedTree->GetEntries(); ++entry){
      pedTree->GetEntry(entry);
      //cout<<" event : "<<_pedCellID->size()<<" "<<_pedMeanCharges->size()<<" "<<_pedMeanTimes->size()<<endl;
      for(size_t ihit=0; ihit!=_pedCellID->size(); ++ihit)
      {
        int layerID = _pedCellID->at(ihit)/1e5;
        int chipID = _pedCellID->at(ihit)%100000/1e4;
        int channelID = _pedCellID->at(ihit)%100;        

        if (layerID>3 && layerID<28) {
          // if(_ChiSqr/(double)_NDF>2) continue;
          chnNoPedHigh10++;
          chnNoPedLow10++;
          meanPedHigh10 += _pedMeanHigh->at(ihit);
          meanPedLow10  += _pedMeanLow ->at(ihit);
          sigmaPedHigh10 += _pedSigmaHigh->at(ihit);
          sigmaPedLow10  += _pedSigmaLow ->at(ihit);
        } else {
           // if(_ChiSqr/(double)_NDF>1.7) continue;
          chnNoPedHigh15++;
          chnNoPedLow15++;
          meanPedHigh15 += _pedMeanHigh->at(ihit);
          meanPedLow15  += _pedMeanLow ->at(ihit);
          sigmaPedHigh15 += _pedSigmaHigh->at(ihit);
          sigmaPedLow15  += _pedSigmaLow ->at(ihit);
        }
        pedMeanHigh[layerID][chipID][channelID] = _pedMeanHigh->at(ihit);
        pedSigmaHigh[layerID][chipID][channelID] = _pedSigmaHigh->at(ihit);
        pedMeanLow[layerID][chipID][channelID] = _pedMeanLow->at(ihit);
        pedSigmaLow[layerID][chipID][channelID] = _pedSigmaLow->at(ihit);
      }
    }
    meanPedHigh15 /= chnNoPedHigh15;
    meanPedLow15  /= chnNoPedLow15;
    meanPedHigh10 /= chnNoPedHigh10;
    meanPedLow10  /= chnNoPedLow10;
    sigmaPedHigh15 /= chnNoPedHigh15;
    sigmaPedLow15  /= chnNoPedLow15;
    sigmaPedHigh10 /= chnNoPedHigh10;
    sigmaPedLow10  /= chnNoPedLow10;
    cout << " HIGH GAIN PEDESTAL ... ..."  << endl;
    cout << " channel count(15um:10um) : " << chnNoPedHigh15 << " , " << chnNoPedHigh10 << endl;
    cout << " mean of pedestal High : "    << meanPedHigh15  << " , " << meanPedHigh10  << endl;
    cout << " sigma of pedestal High: "    << sigmaPedHigh15 << " , " << sigmaPedHigh10 << endl;
    cout << " LOW GAIN PEDESTAL ... ..."   << endl;
    cout << " channel count(15um:10um) : " << chnNoPedLow15  << " , " << chnNoPedLow10  << endl;
    cout << " mean of pedestal Low : "     << meanPedLow15   << " , " << meanPedLow10   << endl;
    cout << " sigma of pedestal Low: "     << sigmaPedLow15  << " , " << sigmaPedLow10  << endl;

    for (int il=0; il<layerNu; il++) {
      for (int ip=0; ip<chipNu; ip++) {
        for (int ic=0; ic<channelNu; ic++) {
          if(pedMeanHigh[il][ip][ic]==1) cout << il << ", " << ip << ", " << ic << endl;
          if (il>3 && il<28) {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh10 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow10 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh10 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow10 : pedSigmaLow[il][ip][ic];
          } else {
            pedMeanHigh[il][ip][ic]  = (pedMeanHigh[il][ip][ic]==1)  ? meanPedHigh15 : pedMeanHigh[il][ip][ic];
            pedMeanLow[il][ip][ic]   = (pedMeanLow[il][ip][ic]==1)    ? meanPedLow15 : pedMeanLow[il][ip][ic];
            pedSigmaHigh[il][ip][ic] = (pedSigmaHigh[il][ip][ic]==1) ? sigmaPedHigh15 : pedSigmaHigh[il][ip][ic];
            pedSigmaLow[il][ip][ic]  = (pedSigmaLow[il][ip][ic]==1)  ? sigmaPedLow15 : pedSigmaLow[il][ip][ic];
          }
        }
      }
    }
  pedFile->Close();
  fileout.cd();
    







  



  // 保存形式
  



  cout << "made save TTree" << endl;



  // vector<string> skipFiles{"ECAL_Run196_20230504_123821.root", "ECAL_Run197_20230504_135137.root", "ECAL_Run194_20230504_082851.root", "ECAL_Run195_20230504_104136.root", "AHCAL_Run133_20230501_081117.root"};


  // data をとってきてる
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    if(Nenegry>1) cout << energy[ienegry] << " GeV" << endl;

    TFile *filein[energy_files[ienegry]];
    TTree *tree[energy_files[ienegry]];
    int entry_max[energy_files[ienegry]];
    for(int i=0; i<energy_files[ienegry]; i++){
      filein[i] = new TFile(argv[irawfilenum+i]);
      tree[i] = (TTree*) filein[i]->Get("MC_Truth");
      entry_max[i] = tree[i]->GetEntries();
    }

    for(int irawfile=0; irawfile<energy_files[ienegry]; irawfile++){
      cout << irawfile << "/" << energy_files[ienegry] << endl;

      // string fileName = argv[irawfilenum+irawfile];
      // auto result = find(skipFiles.begin(), skipFiles.end(), fileName);
      // if (result != skipFiles.end()) continue;
      
      tree[irawfile]->SetBranchAddress("CellID", &cellID);
      tree[irawfile]->SetBranchAddress("Hit_Energy", &Hit_Energy);
      // tree[irawfile]->SetBranchAddress("TotalEnergyDep", &EnergyDep);
      tree[irawfile]->SetBranchAddress("Hit_X", &Hit_X);
      tree[irawfile]->SetBranchAddress("Hit_Y", &Hit_Y);
      tree[irawfile]->SetBranchAddress("Hit_Z", &Hit_Z);

      for(int ientry=0; ientry<entry_max[irawfile]; ientry++){
        // if(ientry%1000==0) cout << ientry << "/" << entry_max[irawfile] << endl;
        tree[irawfile]->GetEntry(ientry);

        double sumEdep = 0;
        double cog_x[32] = {0};
        double cog_y[32] = {0};
        int nhits[32] = {0};
        int n_of_hits[2] = {0};
        double energies[32] = {0};
        EnergyDep = 0;
        vector<double> layerSSA, xSSA, ySSA, eSSA;
        for(int ihit=0;ihit<Hit_Energy->size();ihit++){
          int _layerID = cellID->at(ihit)/1e5;
          if(_layerID>=30) continue;
          if(Hit_Energy->at(ihit)<0) continue;
          EnergyDep += Hit_Energy->at(ihit);
        }
        energy_deposition[ienegry]->Fill(EnergyDep);
        for(int ihit=0;ihit<Hit_Energy->size();ihit++){
          int _layerID = cellID->at(ihit)/1e5;
          int _chipID  = (cellID->at(ihit)%100000) /1e4;
          int _chanID  = cellID->at(ihit)%100;
          if(_layerID>=30) continue;
          if(Hit_Energy->at(ihit)<0) continue;
          layerSSA.push_back(_layerID);
          xSSA.push_back(Hit_X->at(ihit));
          ySSA.push_back(Hit_Y->at(ihit));
          eSSA.push_back(Hit_Energy->at(ihit));
          n_of_hits[0]++;
          const int currentCellId = cellID->at(ihit);
          if(bool_only_best) if(bad_channels.find(currentCellId) != bad_channels.end()) continue;
          sumEdep += Hit_Energy->at(ihit);

          cog_x[_layerID] += Hit_X->at(ihit) * Hit_Energy->at(ihit);
          cog_y[_layerID] += Hit_Y->at(ihit) * Hit_Energy->at(ihit);
          nhits[_layerID]++;
          n_of_hits[1]++;
          energies[_layerID] += Hit_Energy->at(ihit);

          int SiPMType = ((_layerID < 4) || (_layerID >27)) ? 1 : 0;
          // if(MIPPeakValue[_layerID][_chipID][_chanID]<Hit_Energy->at(ihit)) edep_1hit[ienegry][SiPMType]->Fill(Hit_Energy->at(ihit));
          if(thresholdCheck(Hit_Energy->at(ihit), thresholdValue[_layerID][_chipID][_chanID], thresholdSigmaValue[_layerID][_chipID][_chanID], MIPPeakValue[_layerID][_chipID][_chanID])) edep_1hit[ienegry][SiPMType]->Fill(Hit_Energy->at(ihit));
        }
        number_of_hits[ienegry]->Fill(n_of_hits[0]);
        hit_vs_e[ienegry]->Fill(EnergyDep, n_of_hits[0]);
        for(int i=0;i<layerNu;i++) layer_vs_edep[ienegry]->Fill(i, energies[i]);
        // double radius = CalculateShowerRadiusRMS(Hit_X, Hit_Y, Hit_Energy);
        // SSAResult ssa_output = PerformStripSplitSSA(xSSA, ySSA, layerSSA, eSSA);
        // ShowerContainer shower_output = CalculateShowerRadius(ssa_output.x, ssa_output.y, ssa_output.layer, ssa_output.energy);
        // shower_radiusRMS[ienegry]->Fill(shower_output.rms_radius);
        // shower_radius90[ienegry]->Fill(shower_output.r90_radius);
        // // edep_vs_rRMS[ienegry]->Fill(sumEdep, shower_output.rms_radius);
        // // edep_vs_r90[ienegry]->Fill(sumEdep, shower_output.r90_radius);
        // for(int i=0;i<layerNu;i++){
        //   layer_vs_rRMS[ienegry]->Fill(i, shower_output.layers[i].rms_radius);
        //   layer_vs_r90[ienegry]->Fill(i, shower_output.layers[i].r90_radius);
        // }
        
        if( 
          cog_position_check(cog_x[9]/energies[9], cog_y[9]/energies[9], cog_range) &&
          cog_position_check(cog_x[10]/energies[10], cog_y[10]/energies[10], cog_range) 
        ) {
          energy_deposition_bestChannels[ienegry]->Fill(sumEdep);
          number_of_hits_bestChannels[ienegry]->Fill(n_of_hits[1]);
          hit_vs_e_bestChannels[ienegry]->Fill(sumEdep, n_of_hits[1]);
          for(int i=0;i<layerNu;i++) layer_vs_edep_bestChannels[ienegry]->Fill(i, energies[i]);
        }

        for (int ilayer = 0; ilayer < layerNu; ++ilayer) {
          cog[ienegry][ilayer]->Fill(cog_x[ilayer]/energies[ilayer], cog_y[ilayer]/energies[ilayer]);
        }
      }

      EnergyDep = 0;
      cellID = nullptr;
      Hit_Energy = nullptr;
      Hit_X = nullptr;
      Hit_Y = nullptr;
      Hit_Z = nullptr;
      delete filein[irawfile];
    }

    irawfilenum += energy_files[ienegry];
  }

  cout << "fitting histograms ... " << endl;

  fileout.cd();
  gStyle->SetOptFit(1111111);



  TCanvas* C_res[Nenegry];
  TCanvas* C_res_bestChannels[Nenegry];
  TF1 *gaussian = new TF1("gaussian","gaus",200,range_maximum);
  TF1 *gaussian_PS = new TF1("gaussian_PS","gaus",10,range_maximum_PS);
  double energy_mean[Nenegry], energy_sigma[Nenegry], energy_mean_error[Nenegry], energy_sigma_error[Nenegry];
  double energy_mean_bestChannels[Nenegry], energy_sigma_bestChannels[Nenegry], energy_mean_error_bestChannels[Nenegry], energy_sigma_error_bestChannels[Nenegry];
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    C_res[ienegry] = new TCanvas(Form("Edep_%dGeV",(int)energy[ienegry]), Form("%g GeV energy deposit",energy[ienegry]), 1);
    C_res[ienegry]->cd();
    gStyle->SetOptFit(1111111);

    energy_deposition[ienegry]->GetXaxis()->SetRange(200,range_maximum);
    if(energy[ienegry]<=5) energy_deposition[ienegry]->GetXaxis()->SetRange(10/binWidth_PS,range_maximum_PS/binWidth_PS);
    double amp = energy_deposition[ienegry]->GetMaximum();
    double mean = energy_deposition[ienegry]->GetMean();
    double sigma = energy_deposition[ienegry]->GetStdDev();
    double fit_range_low = 15 > mean-3*sigma ? 15 : mean-3*sigma;
    // double mean = energy_deposition[ienegry]->GetMaximumBin();
    // if(energy[ienegry]<=5){
    //   amp = 50;
    //   mean = energy[ienegry] * 60;
    // }
    
    gaussian->SetParameter(0,amp);
    gaussian->SetParameter(1,mean);
    gaussian->SetParameter(2,sigma);
    energy_deposition[ienegry]->Fit("gaussian","","NQ",fit_range_low,mean+3*sigma);
    mean = gaussian->GetParameter(1);
    sigma = gaussian->GetParameter(2);
    fit_range_low = 15 > mean-1.5*sigma ? 15 : mean-1.5*sigma;
    energy_deposition[ienegry]->Fit("gaussian","","NQ",fit_range_low,mean+1.5*sigma);
    // if(energy[ienegry]>5) energy_deposition[ienegry]->Fit("gaussian","","NQ",mean-50*(energy[ienegry]/20),mean+50*(energy[ienegry]/20));
    // else energy_deposition[ienegry]->Fit("gaussian","","NQ",mean-20*energy[ienegry],mean+20*energy[ienegry]);
    // if(energy[ienegry] == 20) energy_deposition[ienegry]->Fit("gaussian","","NQ",1050,1250);

    //////////////////////////////////
    // int histo_range_minimum = gaussian->GetParameter(1) - gaussian->GetParameter(2) * 7 > 15 ? gaussian->GetParameter(1) - gaussian->GetParameter(2) * 7 : 15;
    // int histo_range_maximum = gaussian->GetParameter(1) + gaussian->GetParameter(2) * 7;
    // energy_deposition[ienegry]->GetXaxis()->SetRange(histo_range_minimum,histo_range_maximum);
    //////////////////////////////////

    energy_deposition[ienegry]->Draw();
    gaussian->Draw("same");
    fileout.cd("all__channel");
    C_res[ienegry]->Write();
    fileout.cd();
    energy_mean[ienegry] = gaussian->GetParameter(1);
    energy_sigma[ienegry] = gaussian->GetParameter(2);
    energy_mean_error[ienegry] = gaussian->GetParError(1);
    energy_sigma_error[ienegry] = gaussian->GetParError(2);





    C_res_bestChannels[ienegry] = new TCanvas(Form("Edep_%dGeV_bestChannels",(int)energy[ienegry]), Form("%g GeV energy deposit",energy[ienegry]), 1);
    C_res_bestChannels[ienegry]->cd();
    gStyle->SetOptFit(1111111);

    energy_deposition_bestChannels[ienegry]->GetXaxis()->SetRange(200,range_maximum);
    if(energy[ienegry]<=5) energy_deposition_bestChannels[ienegry]->GetXaxis()->SetRange(10/binWidth_PS,range_maximum_PS/binWidth_PS);
    amp = energy_deposition_bestChannels[ienegry]->GetMaximum();
    mean = energy_deposition_bestChannels[ienegry]->GetMean();
    sigma = energy_deposition_bestChannels[ienegry]->GetStdDev();
    fit_range_low = 15 > mean-3*sigma ? 15 : mean-3*sigma;
    // amp = energy_deposition_bestChannels[ienegry]->GetMaximum();
    // mean = energy_deposition_bestChannels[ienegry]->GetMaximumBin();
    // if(energy[ienegry]<=5){
    //   amp = 50;
    //   mean = energy[ienegry] * 60;
    // }

    gaussian->SetParameter(0,amp);
    gaussian->SetParameter(1,mean);
    gaussian->SetParameter(2,sigma);
    energy_deposition_bestChannels[ienegry]->Fit("gaussian","","NQ",fit_range_low,mean+3*sigma);
    mean = gaussian->GetParameter(1);
    sigma = gaussian->GetParameter(2);
    fit_range_low = 15 > mean-1.5*sigma ? 15 : mean-1.5*sigma;
    energy_deposition_bestChannels[ienegry]->Fit("gaussian","","NQ",fit_range_low,mean+1.5*sigma);
    // if(energy[ienegry]>5) energy_deposition_bestChannels[ienegry]->Fit("gaussian","","NQ",mean-50*(energy[ienegry]/20),mean+50*(energy[ienegry]/20));
    // else energy_deposition_bestChannels[ienegry]->Fit("gaussian","","NQ",mean-20*energy[ienegry],mean+20*energy[ienegry]);
    // if(energy[ienegry] == 20) energy_deposition[ienegry]->Fit("gaussian","","NQ",1050,1250);
    energy_deposition_bestChannels[ienegry]->Draw();
    gaussian->Draw("same");
    fileout.cd("best_channel");
    C_res_bestChannels[ienegry]->Write();
    fileout.cd();
    energy_mean_bestChannels[ienegry] = gaussian->GetParameter(1);
    energy_sigma_bestChannels[ienegry] = gaussian->GetParameter(2);
    energy_mean_error_bestChannels[ienegry] = gaussian->GetParError(1);
    energy_sigma_error_bestChannels[ienegry] = gaussian->GetParError(2);
  }


  fileout.cd("all__channel");

  TCanvas *linear_plots = new TCanvas("linear_plots", "linear_plots", 2560,1440); 
  linear_plots->Divide(2,2);
  TCanvas *resolution_plots = new TCanvas("resolution_plots", "resolution_plots", 2560,1440); 
  resolution_plots->Divide(2,2);

  TGraphErrors *Eres = new TGraphErrors();
  TGraphErrors *Elinear = new TGraphErrors();
  Eres->SetTitle("Energy Resolution");
  Eres->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres->GetYaxis()->SetTitle("Energy Resolution");
  Elinear->SetTitle("Energy Linearity (all channels)");
  Elinear->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Elinear->GetYaxis()->SetTitle("Detected Energy [MeV]");

  TGraphErrors *Eres_lowEnergy = new TGraphErrors();
  Eres_lowEnergy->SetTitle("Energy Resolution (all channels)");
  Eres_lowEnergy->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres_lowEnergy->GetYaxis()->SetTitle("Detected Energy Resolution sigma(E)/E");
  int Eres_lowEnergy_count = 0;

  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    if(energy[ienegry] != 200){
      Elinear->SetPoint(ienegry, energy[ienegry], energy_mean[ienegry]);
      Elinear->SetPointError(ienegry, 0, energy_mean_error[ienegry]);
    }
    Eres->SetPoint(ienegry, energy[ienegry], energy_sigma[ienegry]/energy_mean[ienegry]);
    double res_error = sqrt( pow(energy_sigma_error[ienegry]/energy_mean[ienegry],2) + pow(energy_sigma[ienegry]*energy_mean_error[ienegry]/energy_mean[ienegry]/energy_mean[ienegry],2));
    Eres->SetPointError(ienegry, 0, res_error);

    if(energy[ienegry] == 80 || energy[ienegry]>150) continue;
    Eres_lowEnergy->SetPoint(Eres_lowEnergy_count, energy[ienegry], energy_sigma[ienegry]/energy_mean[ienegry]);
    Eres_lowEnergy->SetPointError(Eres_lowEnergy_count++, 0, res_error);
  }
  TF1 *res_func = new TF1("res_func","[0]+[1]/sqrt(x)",1,300);
  res_func->SetParameter(0,0);
  res_func->SetParameter(1,0.2);
  Eres->Fit("res_func");

  TF1 *res_func_low_e = new TF1("res_func_low_e","sqrt([0]^2/x+[1]^2+[2]^2/x^2)",1,300);
  res_func_low_e->SetParameter(0,0.2);
  res_func_low_e->SetParameter(1,0.02);
  res_func_low_e->SetParameter(2,0.2);
  Eres_lowEnergy->Fit("res_func_low_e");

  TF1 *res_func_low = new TF1("res_func_low","sqrt([0]^2/x+[1]^2)",1,300);
  res_func_low->SetParameter(0,0.2);
  res_func_low->SetParameter(1,0);
  Eres_lowEnergy->Fit("res_func_low");
  resolution_plots->cd();
  resolution_plots->cd(1);
  resolution_plots->SetGrid();
  gStyle->SetOptFit(1111111);
  Eres_lowEnergy->SetMinimum(0);
  Eres_lowEnergy->SetMaximum(0.6);
  Eres_lowEnergy->Draw("AP");
  TLegend *legend_resolution = new TLegend(0.5, 0.4, 0.9, 0.55);
  legend_resolution->AddEntry(res_func_low, Form("#frac{%.3f}{#sqrt{E}} #oplus %.3f", abs(res_func_low->GetParameter(0)), abs(res_func_low->GetParameter(1))) , "l");
  legend_resolution->Draw("same");

  Elinear->Write("Elinear");
  Eres->Write("Eres");
  res_func->Write();

  Eres_lowEnergy->Write("Eres_lowEnergy");
  res_func_low->Write();
  res_func_low_e->Write();

  TGraphErrors *Eres_recoE = new TGraphErrors();
  Eres_recoE->SetTitle("Energy Resolution (all channels)");
  Eres_recoE->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres_recoE->GetYaxis()->SetTitle("Reconstructed Energy Resolution sigma(E)/E");
  int Eres_recoE_count = 0;
  TF1 *e_linear_ps = new TF1("e_linear_ps","[0]+[1]*x",0,300);
  Elinear->Fit("e_linear_ps","","",0.7,7);
  double reco_inter_ps = e_linear_ps->GetParameter(0);
  double reco_slope_ps = e_linear_ps->GetParameter(1);
  TF1 *e_linear_sps = new TF1("e_linear_sps","[0]+[1]*x",0,300);
  Elinear->Fit("e_linear_sps","","",7,110);
  double reco_inter_sps = e_linear_sps->GetParameter(0);
  double reco_slope_sps = e_linear_sps->GetParameter(1);
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    double res_error = sqrt( pow(energy_sigma_error[ienegry]/energy_mean[ienegry],2) + pow(energy_sigma[ienegry]*energy_mean_error[ienegry]/energy_mean[ienegry]/energy_mean[ienegry],2));
    if(energy[ienegry] == 20 || energy[ienegry] == 80 || energy[ienegry]>150) continue;
    double slope = energy[ienegry] < 7 ? reco_slope_ps : reco_slope_sps;
    double intecept = energy[ienegry] < 7 ? reco_inter_ps : reco_inter_sps;
    double reco_mean = (energy_mean[ienegry] - intecept) / slope;
    double reco_sigma = energy_sigma[ienegry] / slope;
    Eres_recoE->SetPoint(Eres_recoE_count, energy[ienegry], reco_sigma/reco_mean);
    Eres_recoE->SetPointError(Eres_recoE_count++, 0, res_error);
  }
  Eres_recoE->Fit("res_func_low");
  Eres_recoE->Write("Eres_recoE");
  e_linear_ps->Write("e_linear_ps");
  e_linear_sps->Write("e_linear_sps");

  resolution_plots->cd();
  resolution_plots->cd(2);
  resolution_plots->SetGrid();
  gStyle->SetOptFit(1111111);
  Eres_recoE->SetMinimum(0);
  Eres_recoE->SetMaximum(0.6);
  Eres_recoE->Draw("AP");
  TLegend *legend_resolution_recoE = new TLegend(0.5, 0.4, 0.9, 0.55);
  legend_resolution_recoE->AddEntry(res_func_low, Form("#frac{%.3f}{#sqrt{E}} #oplus %.3f", abs(res_func_low->GetParameter(0)), abs(res_func_low->GetParameter(1))) , "l");
  legend_resolution_recoE->Draw("same");

  fileout.cd();
  

  TCanvas *energy_distributions[3];
  for(int i=0;i<3;i++){
    energy_distributions[i] = new TCanvas(Form("energy_distributions_%d",i), Form("energy_distributions_%d",i), 2560,1440); 
    energy_distributions[i]->Divide(4,2);
  }

  vector<pair<double, int>> energy_itr;
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    energy_itr.push_back({energy[ienegry], ienegry});
  }
  sort(energy_itr.begin(), energy_itr.end());
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    energy_distributions[icanvas]->cd();
    energy_distributions[icanvas]->cd(idivision);
    gStyle->SetOptFit(1111111);
    energy_deposition[ienegry]->Draw();
  }
  for(int i=0;i<3;i++){
    energy_distributions[i]->Write();
    energy_distributions[i]->SaveAs(Form("%s/Edep_%d.png",argv[argc-1],i));
    delete energy_distributions[i];
  }


  linear_plots->cd();
  linear_plots->cd(1);
  gStyle->SetOptFit(1111111);
  Elinear->Draw("AP");
  e_linear_ps->SetLineColor(kBlue);
  e_linear_ps->Draw("same");
  e_linear_sps->SetLineColor(kRed);
  e_linear_sps->Draw("same");
  TLegend *legend_linear = new TLegend(0.4, 0.1, 0.8, 0.4);
  legend_linear->AddEntry(e_linear_ps, Form("ps fit  : %.3f + %.3f x", e_linear_ps->GetParameter(0), e_linear_ps->GetParameter(1)) , "l");
  legend_linear->AddEntry(e_linear_sps, Form("sps fit : %.3f + %.3f x", e_linear_sps->GetParameter(0), e_linear_sps->GetParameter(1)) , "l");
  legend_linear->Draw("same");



  // best channels の保存をしている

  fileout.cd("best_channel");

  TGraphErrors *Eres_bestChannels = new TGraphErrors();
  TGraphErrors *Elinear_bestChannels = new TGraphErrors();
  Eres_bestChannels->SetTitle("Energy Resolution");
  Eres_bestChannels->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres_bestChannels->GetYaxis()->SetTitle("Energy Resolution");
  Eres_bestChannels->SetMinimum(0);
  Elinear_bestChannels->SetTitle("Energy Linearity (best channels)");
  Elinear_bestChannels->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Elinear_bestChannels->GetYaxis()->SetTitle("Detected Energy [MeV]");

  TGraphErrors *Eres_lowEnergy_bestChannels = new TGraphErrors();
  Eres_lowEnergy_bestChannels->SetTitle("Energy Resolution (best channels)");
  Eres_lowEnergy_bestChannels->SetMinimum(0);
  Eres_lowEnergy_bestChannels->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres_lowEnergy_bestChannels->GetYaxis()->SetTitle("Energy Resolution sigma(E)/E");
  int Eres_lowEnergy_count_bestChannels = 0;

  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    if(energy[ienegry] != 200){
      Elinear_bestChannels->SetPoint(ienegry, energy[ienegry], energy_mean_bestChannels[ienegry]);
      Elinear_bestChannels->SetPointError(ienegry, 0, energy_mean_error_bestChannels[ienegry]);
    }
    Eres_bestChannels->SetPoint(ienegry, energy[ienegry], energy_sigma_bestChannels[ienegry]/energy_mean_bestChannels[ienegry]);
    double res_error_bestChannels = sqrt( pow(energy_sigma_error_bestChannels[ienegry]/energy_mean_bestChannels[ienegry],2) + pow(energy_sigma_bestChannels[ienegry]*energy_mean_error_bestChannels[ienegry]/energy_mean_bestChannels[ienegry]/energy_mean_bestChannels[ienegry],2));
    Eres_bestChannels->SetPointError(ienegry, 0, res_error_bestChannels);

    if(energy[ienegry] == 80 || energy[ienegry]>150) continue;
    Eres_lowEnergy_bestChannels->SetPoint(Eres_lowEnergy_count_bestChannels, energy[ienegry], energy_sigma_bestChannels[ienegry]/energy_mean_bestChannels[ienegry]);
    Eres_lowEnergy_bestChannels->SetPointError(Eres_lowEnergy_count_bestChannels++, 0, res_error_bestChannels);
  }
  TF1 *res_func_bestChannels = new TF1("res_func_bestChannels","[0]+[1]/sqrt(x)",1,300);
  res_func_bestChannels->SetParameter(0,0);
  res_func_bestChannels->SetParameter(1,0.2);
  Eres_bestChannels->Fit("res_func_bestChannels");
  Elinear_bestChannels->Fit("pol1");

  TF1 *res_func_low_e_bestChannels = new TF1("res_func_low_e_bestChannels","sqrt([0]^2/x+[1]^2+[2]^2/x^2)",1,300);
  res_func_low_e_bestChannels->SetParameter(0,0.2);
  res_func_low_e_bestChannels->SetParameter(1,0.02);
  res_func_low_e_bestChannels->SetParameter(2,0.2);
  Eres_lowEnergy_bestChannels->Fit("res_func_low_e_bestChannels");

  TF1 *res_func_low_bestChannels = new TF1("res_func_low_bestChannels","sqrt([0]^2/x+[1]^2)",1,300);
  res_func_low_bestChannels->SetParameter(0,0.2);
  res_func_low_bestChannels->SetParameter(1,0);
  Eres_lowEnergy_bestChannels->Fit("res_func_low_bestChannels");
  resolution_plots->cd();
  resolution_plots->cd(3);
  resolution_plots->SetGrid();
  Eres_lowEnergy_bestChannels->SetMinimum(0);
  Eres_lowEnergy_bestChannels->SetMaximum(0.6);
  gStyle->SetOptFit(1111111);
  Eres_lowEnergy_bestChannels->Draw("AP");
  TLegend *legend_bestChannels_resolution = new TLegend(0.5, 0.4, 0.9, 0.55);
  legend_bestChannels_resolution->AddEntry(res_func_low_bestChannels, Form("#frac{%.3f}{#sqrt{E}} #oplus %.3f", abs(res_func_low_bestChannels->GetParameter(0)), abs(res_func_low_bestChannels->GetParameter(1))) , "l");
  legend_bestChannels_resolution->Draw("same");

  Elinear_bestChannels->Write("Elinear_bestChannels");
  Eres_bestChannels->Write("Eres");
  res_func_bestChannels->Write();

  Eres_lowEnergy_bestChannels->Write("Eres_lowEnergy_bestChannel");
  res_func_low_bestChannels->Write();
  res_func_low_e_bestChannels->Write();


  TGraphErrors *Eres_recoE_bestChannels = new TGraphErrors();
  Eres_recoE_bestChannels->SetTitle("Energy Resolution (best channels)");
  Eres_recoE_bestChannels->GetXaxis()->SetTitle("Beam Energy [GeV]");
  Eres_recoE_bestChannels->GetYaxis()->SetTitle("Reconstructed Energy Resolution sigma(E)/E");
  int Eres_recoE_count_bestChannels = 0;
  TF1 *e_linear_ps_bestChannels = new TF1("e_linear_ps_bestChannels","[0]+[1]*x",0,300);
  Elinear_bestChannels->Fit("e_linear_ps_bestChannels","","",0.7,7);
  double reco_inter_ps_bestChannels = e_linear_ps_bestChannels->GetParameter(0);
  double reco_slope_ps_bestChannels = e_linear_ps_bestChannels->GetParameter(1);
  TF1 *e_linear_sps_bestChannels = new TF1("e_linear_sps_bestChannels","[0]+[1]*x",0,300);
  Elinear_bestChannels->Fit("e_linear_sps_bestChannels","","",7,110);
  double reco_inter_sps_bestChannels = e_linear_sps_bestChannels->GetParameter(0);
  double reco_slope_sps_bestChannels = e_linear_sps_bestChannels->GetParameter(1);
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    double res_error_bestChannels = sqrt( pow(energy_sigma_error_bestChannels[ienegry]/energy_mean_bestChannels[ienegry],2) + pow(energy_sigma_bestChannels[ienegry]*energy_mean_error_bestChannels[ienegry]/energy_mean_bestChannels[ienegry]/energy_mean_bestChannels[ienegry],2));
    if(energy[ienegry] == 20 || energy[ienegry] == 80 || energy[ienegry]>150) continue;
    double slope = energy[ienegry] < 7 ? reco_slope_ps_bestChannels : reco_slope_sps_bestChannels;
    double intecept = energy[ienegry] < 7 ? reco_inter_ps_bestChannels : reco_inter_sps_bestChannels;
    double reco_mean = (energy_mean_bestChannels[ienegry] - intecept) / slope;
    double reco_sigma = energy_sigma_bestChannels[ienegry] / slope;
    Eres_recoE_bestChannels->SetPoint(Eres_recoE_count_bestChannels, energy[ienegry], reco_sigma/reco_mean);
    Eres_recoE_bestChannels->SetPointError(Eres_recoE_count_bestChannels++, 0, res_error_bestChannels);
  }
  Eres_recoE_bestChannels->Fit("res_func_low_bestChannels");
  Eres_recoE_bestChannels->Write("Eres_recoE_bestChannels");
  e_linear_ps_bestChannels->Write("e_linear_ps_bestChannels");
  e_linear_sps_bestChannels->Write("e_linear_sps_bestChannels");

  resolution_plots->cd();
  resolution_plots->cd(4);
  gStyle->SetOptFit(1111111);
  resolution_plots->SetGrid();
  Eres_recoE_bestChannels->SetMinimum(0);
  Eres_recoE_bestChannels->SetMaximum(0.6);
  Eres_recoE_bestChannels->Draw("AP");
  TLegend *legend_bestChannels_resolution_recoE = new TLegend(0.5, 0.4, 0.9, 0.55);
  legend_bestChannels_resolution_recoE->AddEntry(res_func_low_bestChannels, Form("#frac{%.3f}{#sqrt{E}} #oplus %.3f", abs(res_func_low_bestChannels->GetParameter(0)), abs(res_func_low_bestChannels->GetParameter(1))) , "l");
  legend_bestChannels_resolution_recoE->Draw("same");
  

  fileout.cd();

  TCanvas *energy_distributions_bestChannels[3];
  for(int i=0;i<3;i++){
    energy_distributions_bestChannels[i] = new TCanvas(Form("energy_distributions_bestChannels_%d",i), Form("energy_distributions_bestChannels_%d",i), 2560,1440); 
    energy_distributions_bestChannels[i]->Divide(4,2);
  }

  vector<pair<double, int>> energy_itr_bestChannels;
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    energy_itr_bestChannels.push_back({energy[ienegry], ienegry});
  }
  sort(energy_itr_bestChannels.begin(), energy_itr_bestChannels.end());
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    energy_distributions_bestChannels[icanvas]->cd();
    energy_distributions_bestChannels[icanvas]->cd(idivision);
    gStyle->SetOptFit(1111111);
    energy_deposition_bestChannels[ienegry]->Draw();
  }
  for(int i=0;i<3;i++){
    energy_distributions_bestChannels[i]->Write();
    energy_distributions_bestChannels[i]->SaveAs(Form("%s/Edep_best_%d.png",argv[argc-1],i));
    delete energy_distributions_bestChannels[i];
  }


  /*


  linear_plots->cd();
  linear_plots->cd(2);
  gStyle->SetOptFit(1111111);
  Elinear_bestChannels->Draw("AP");
  e_linear_ps_bestChannels->SetLineColor(kBlue);
  e_linear_ps_bestChannels->Draw("same");
  e_linear_sps_bestChannels->SetLineColor(kRed);
  e_linear_sps_bestChannels->Draw("same");
  TLegend *legend_linear_bestChannels = new TLegend(0.4, 0.1, 0.8, 0.4);
  legend_linear_bestChannels->AddEntry(e_linear_ps_bestChannels, Form("ps fit  : %.3f + %.3f x", e_linear_ps_bestChannels->GetParameter(0), e_linear_ps_bestChannels->GetParameter(1)) , "l");
  legend_linear_bestChannels->AddEntry(e_linear_sps_bestChannels, Form("sps fit : %.3f + %.3f x", e_linear_sps_bestChannels->GetParameter(0), e_linear_sps_bestChannels->GetParameter(1)) , "l");
  legend_linear_bestChannels->Draw("same");




  linear_plots->Write("linear_plots");
  resolution_plots->Write("resolution_plots");




  linear_plots->SaveAs(Form("%s/linearity.png",argv[argc-1]));
  resolution_plots->SaveAs(Form("%s/resolution.png",argv[argc-1]));









  fileout.cd("all__channel");
  TCanvas *C_cog[Nenegry];
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    C_cog[ienegry] = new TCanvas(Form("C_cog_%d",ienegry), Form("center of gravity %g GeV",energy[ienegry]), 2560,1440); 
    C_cog[ienegry]->Divide(6,5);
    C_cog[ienegry]->cd();
    for (int ilayer = 0; ilayer < 30; ++ilayer) {
      C_cog[ienegry]->cd(ilayer+1);
      cog[ienegry][ilayer]->Draw("colz");
    }
    C_cog[ienegry]->Write(Form("C_cog_%dGeV",(int)energy[ienegry]));
    delete C_cog[ienegry];
  }

  TCanvas *C_nhits[3];
  for(int i=0;i<3;i++){
    C_nhits[i] = new TCanvas(Form("C_nhits_%d",i), Form("C_nhits_%d",i), 2560,1440); 
    C_nhits[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_nhits[icanvas]->cd();
    C_nhits[icanvas]->cd(idivision);
    gStyle->SetOptFit(1111111);
    number_of_hits[ienegry]->Draw();
  }
  for(int i=0;i<3;i++){
    C_nhits[i]->Write();
    C_nhits[i]->SaveAs(Form("%s/numberof_hits_%d.png",argv[argc-1],i));
    delete C_nhits[i];
  }

  TCanvas* C_ledep[3];
  for(int i=0;i<3;i++){
    C_ledep[i] = new TCanvas(Form("C_ledep_%d",i), Form("C_ledep_%d",i), 2560,1440); 
    C_ledep[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_ledep[icanvas]->cd();
    C_ledep[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogz(1);
    layer_vs_edep[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_ledep[i]->Write();
    C_ledep[i]->SaveAs(Form("%s/layer_vs_edep_%d.png",argv[argc-1],i));
    delete C_ledep[i];
  }

  TCanvas* C_rRMS[3];
  for(int i=0;i<3;i++){
    C_rRMS[i] = new TCanvas(Form("C_rRMS_%d",i), Form("C_rRMS_%d",i), 2560,1440); 
    C_rRMS[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;
    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_rRMS[icanvas]->cd();
    C_rRMS[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    shower_radiusRMS[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_rRMS[i]->Write();
    C_rRMS[i]->SaveAs(Form("%s/shower_rRMS%d.png",argv[argc-1],i));
    delete C_rRMS[i];
  }

  TCanvas* C_r90[3];
  for(int i=0;i<3;i++){
    C_r90[i] = new TCanvas(Form("C_r90_%d",i), Form("C_r90_%d",i), 2560,1440); 
    C_r90[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;
    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_r90[icanvas]->cd();
    C_r90[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);
    shower_radius90[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_r90[i]->Write();
    C_r90[i]->SaveAs(Form("%s/shower_r90%d.png",argv[argc-1],i));
    delete C_r90[i];
  }

  TCanvas* C_lrRMS[3];
  for(int i=0;i<3;i++){
    C_lrRMS[i] = new TCanvas(Form("C_lrRMS_%d",i), Form("C_lrRMS_%d",i), 2560,1440); 
    C_lrRMS[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;
    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_lrRMS[icanvas]->cd();
    C_lrRMS[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogz(1);
    layer_vs_rRMS[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_lrRMS[i]->Write();
    C_lrRMS[i]->SaveAs(Form("%s/layer_vs_rRMS_%d.png",argv[argc-1],i));
    delete C_lrRMS[i];
  }

  TCanvas* C_lr90[3];
  for(int i=0;i<3;i++){
    C_lr90[i] = new TCanvas(Form("C_lr90_%d",i), Form("C_lr90_%d",i), 2560,1440); 
    C_lr90[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;
    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_lr90[icanvas]->cd();
    C_lr90[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogz(1);
    layer_vs_r90[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_lr90[i]->Write();
    C_lr90[i]->SaveAs(Form("%s/layer_vs_r90_%d.png",argv[argc-1],i));
    delete C_lr90[i];
  }
  


  fileout.cd();
  fileout.cd("best_channel");
  TCanvas *C_nhits_bestChannels[3];
  for(int i=0;i<3;i++){
    C_nhits_bestChannels[i] = new TCanvas(Form("C_nhits_bestChannels_%d",i), Form("C_nhits_bestChannels_%d",i), 2560,1440); 
    C_nhits_bestChannels[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_nhits_bestChannels[icanvas]->cd();
    C_nhits_bestChannels[icanvas]->cd(idivision);
    gStyle->SetOptFit(1111111);
    number_of_hits_bestChannels[ienegry]->Draw();
  }
  for(int i=0;i<3;i++){
    C_nhits_bestChannels[i]->Write();
    C_nhits_bestChannels[i]->SaveAs(Form("%s/numberof_hits_bestChannels_%d.png",argv[argc-1],i));
    delete C_nhits_bestChannels[i];
  }





  fileout.cd();
  TCanvas* C_e_summary[Nenegry];
  for(int ienegry=0; ienegry<Nenegry; ienegry++){
    C_e_summary[ienegry] = new TCanvas(Form("Edep_%dGeV_summary",(int)energy[ienegry]), Form("%g GeV",energy[ienegry]), 2560,1440); 
    C_e_summary[ienegry]->cd();
    C_e_summary[ienegry]->Divide(4,2);

    C_e_summary[ienegry]->cd(1);
    hit_vs_e[ienegry]->Draw("colz");
    C_e_summary[ienegry]->cd(2);
    gStyle->SetOptFit(1111111);
    energy_deposition[ienegry]->Draw();
    C_e_summary[ienegry]->cd(3);
    number_of_hits[ienegry]->Draw();
    C_e_summary[ienegry]->cd(4);
    gPad->SetLogz(1);
    layer_vs_edep[ienegry]->Draw("colz");
    

    C_e_summary[ienegry]->cd(5);
    hit_vs_e_bestChannels[ienegry]->Draw("colz");
    C_e_summary[ienegry]->cd(6);
    gStyle->SetOptFit(1111111);
    energy_deposition_bestChannels[ienegry]->Draw();
    C_e_summary[ienegry]->cd(7);
    number_of_hits_bestChannels[ienegry]->Draw();
    C_e_summary[ienegry]->cd(8);
    gPad->SetLogz(1);
    layer_vs_edep_bestChannels[ienegry]->Draw("colz");

    C_e_summary[ienegry]->Write();
    C_e_summary[ienegry]->SaveAs(Form("%s/summary_%dGeV.png",argv[argc-1],(int)energy[ienegry]));
    delete C_e_summary[ienegry];
  }


  fileout.cd("best_channel");
  TCanvas* C_ledep_bestChannels[3];
  for(int i=0;i<3;i++){
    C_ledep_bestChannels[i] = new TCanvas(Form("C_ledep_bestChannels_%d",i), Form("C_ledep_bestChannels_%d",i), 2560,1440); 
    C_ledep_bestChannels[i]->Divide(4,2);
  }
  for(int iienegry=0; iienegry<Nenegry; iienegry++){
    int ienegry = energy_itr_bestChannels[iienegry].second;

    int icanvas = iienegry/8;
    int idivision = iienegry%8 + 1;
    C_ledep_bestChannels[icanvas]->cd();
    C_ledep_bestChannels[icanvas]->cd(idivision);
    gStyle->SetOptStat(0);
    gPad->SetLogz(1);
    layer_vs_edep_bestChannels[ienegry]->Draw("colz");
  }
  for(int i=0;i<3;i++){
    C_ledep_bestChannels[i]->Write();
    C_ledep_bestChannels[i]->SaveAs(Form("%s/layer_vs_edep_best_%d.png",argv[argc-1],i));
    delete C_ledep_bestChannels[i];
  }


  */


  fileout.cd();
  TCanvas* C_edep_1hit[3][2];
  TCanvas* C_edep_1hit_zoom[3][2];
  for(int sipmtype=0;sipmtype<2;sipmtype++){
    const int sipmtype_int = sipmtype == 1 ? 15 : 10;
    for(int i=0;i<3;i++){
      C_edep_1hit[i][sipmtype] = new TCanvas(Form("C_edep_1hit_%dum_%d",sipmtype_int,i), Form("C_edep_1hit_%dum_%d",sipmtype_int,i), 2560,1440); 
      C_edep_1hit[i][sipmtype]->Divide(4,2);
      C_edep_1hit_zoom[i][sipmtype] = new TCanvas(Form("C_edep_1hit_%dum_zoom_%d",sipmtype_int,i), Form("C_edep_1hit_%dum_zoom_%d",sipmtype_int,i), 2560,1440); 
      C_edep_1hit_zoom[i][sipmtype]->Divide(4,2);
    }
    for(int iienegry=0; iienegry<Nenegry; iienegry++){
      int ienegry = energy_itr_bestChannels[iienegry].second;
  
      int icanvas = iienegry/8;
      int idivision = iienegry%8 + 1;
  
      C_edep_1hit[icanvas][sipmtype]->cd();
      C_edep_1hit[icanvas][sipmtype]->cd(idivision);
      // gStyle->SetOptStat(0);
      // gPad->SetLogz(1);
      edep_1hit[ienegry][sipmtype]->Draw("colz");
  
      C_edep_1hit_zoom[icanvas][sipmtype]->cd();
      C_edep_1hit_zoom[icanvas][sipmtype]->cd(idivision);
      // gStyle->SetOptStat(0);
      // gPad->SetLogz(1);
      edep_1hit[ienegry][sipmtype]->GetXaxis()->SetRangeUser(0,2);
      edep_1hit[ienegry][sipmtype]->Draw("colz");
    }
    for(int i=0;i<3;i++){
      C_edep_1hit[i][sipmtype]->Write();
      C_edep_1hit[i][sipmtype]->SaveAs(Form("%s/edep_1hit_%dum_%d.png",argv[argc-1],sipmtype_int,i));
  
      C_edep_1hit_zoom[i][sipmtype]->Write();
      C_edep_1hit_zoom[i][sipmtype]->SaveAs(Form("%s/edep_1hit_%dum_zoom_%d.png",argv[argc-1],sipmtype_int,i));

      delete C_edep_1hit[i][sipmtype];
      delete C_edep_1hit_zoom[i][sipmtype];
    }
  }
  





}