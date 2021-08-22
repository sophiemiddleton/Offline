#

import sys, glob



#------------------------------------------------------------------------------
def fix_filename(fn,first_run,first_subrun):

    w    = fn.split('.')
    w[4] = "%06i_%08i"%(first_run,first_subrun)
    fn1  = '.'.join(w);

    return fn1
    
fix_data = [
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00000000.art.event_list",   7371,     0,  2011,   134,  2007,   272,  2008,   410,  1345),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00000503.art.event_list",   7485,   503,  2012,   638,  2003,   774,  2019,   910,  1451),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00001005.art.event_list",   7373,  1005,  2003,  1142,  2012,  1279,  2004,  1416,  1354),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00001505.art.event_list",   7271,  1505,  2019,  1646,  2019,  1784,  2017,  1922,  1216),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00002006.art.event_list",   7285,  2006,  2004,  2141,  2016,  2276,  2012,  2420,  1253),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00002507.art.event_list",   7278,  2507,  2009,  2644,  2004,  2780,  2015,  2919,  1250),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00003010.art.event_list",   7319,  3010,  2001,  3147,  2001,  3282,  2017,  3421,  1300),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00003511.art.event_list",   7447,  3511,  2001,  3645,  2008,  3780,  2005,  3916,  1433),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00004012.art.event_list",   7313,  4012,  2019,  4154,  2004,  4289,  2009,  4424,  1281),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00004513.art.event_list",   7382,  4513,  2001,  4647,  2010,  4784,  2007,  4927,  1364),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00005015.art.event_list",   7360,  5015,  2005,  5152,  2015,  5290,  2006,  5426,  1334),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00005516.art.event_list",   7593,  5516,  2011,  5649,  2028,  5780,  2002,  5912,  1552),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00006017.art.event_list",   7356,  6017,  2010,  6151,  2010,  6286,  2006,  6430,  1330),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00006517.art.event_list",   7241,  6517,  2001,  6653,  2003,  6793,  2002,  6932,  1235),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00007019.art.event_list",   7398,  7019,  2010,  7156,  2001,  7289,  2015,  7428,  1372),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00007521.art.event_list",   7428,  7521,  2013,  7656,  2002,  7796,  2003,  7926,  1410),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00008022.art.event_list",   7435,  8022,  2008,  8156,  2011,  8292,  2004,  8428,  1412),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00008522.art.event_list",   7415,  8522,  2017,  8659,  2007,  8792,  2006,  8930,  1385),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00009027.art.event_list",   7493,  9027,  2004,  9165,  2012,  9297,  2011,  9427,  1466),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581292/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00009528.art.event_list",   6983,  9528,  2013,  9669,  2010,  9807,  2009,  9943,   951),

("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00010008.art.event_list",   7431, 10008,  2017, 10148,  2007, 10282,  2004, 10415,  1403),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00010511.art.event_list",   7335, 10511,  2001, 10647,  2009, 10785,  2016, 10924,  1309),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00011011.art.event_list",   7376, 11011,  2006, 11146,  2004, 11280,  2006, 11417,  1360),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00011512.art.event_list",   7417, 11512,  2007, 11651,  2005, 11785,  2014, 11919,  1391),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00012013.art.event_list",   7192, 12013,  2010, 12154,  2011, 12298,  2007, 12436,  1164),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00012514.art.event_list",   7294, 12514,  2011, 12655,  2002, 12789,  2007, 12933,  1274),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00013016.art.event_list",   7409, 13016,  2018, 13151,  2003, 13286,  2009, 13424,  1379),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00013519.art.event_list",   7485, 13519,  2007, 13656,  2014, 13793,  2002, 13924,  1462),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00014019.art.event_list",   7443, 14019,  2013, 14156,  2010, 14295,  2015, 14426,  1405),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00014521.art.event_list",   7428, 14521,  2012, 14660,  2004, 14794,  2011, 14930,  1401),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00015022.art.event_list",   7323, 15022,  2013, 15159,  2001, 15295,  2009, 15437,  1300),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00015524.art.event_list",   7412, 15524,  2001, 15655,  2005, 15794,  2002, 15931,  1404),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00016025.art.event_list",   7418, 16025,  2002, 16164,  2001, 16297,  2004, 16431,  1411),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00016527.art.event_list",   7553, 16527,  2011, 16651,  2008, 16782,  2007, 16922,  1527),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00017027.art.event_list",   7225, 17027,  2015, 17167,  2005, 17308,  2012, 17448,  1193),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00017531.art.event_list",   7383, 17531,  2003, 17665,  2008, 17801,  2003, 17935,  1369),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00018031.art.event_list",   7316, 18031,  2016, 18170,  2003, 18309,  2018, 18444,  1279),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00018531.art.event_list",   7273, 18531,  2009, 18669,  2012, 18807,  2001, 18945,  1251),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00019032.art.event_list",   7297, 19032,  2013, 19173,  2015, 19305,  2010, 19444,  1259),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35581344/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00019533.art.event_list",   6828, 19533,  2005, 19672,  2007, 19814,  2008, 19952,   808),

("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00020006.art.event_list",   7382, 20006,  2014, 20146,  2004, 20280,  2011, 20413,  1353),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00020507.art.event_list",   7330, 20507,  2009, 20648,  2007, 20785,  2002, 20917,  1312),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00021012.art.event_list",   7288, 21012,  2009, 21148,  2016, 21285,  2001, 21427,  1262),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00021517.art.event_list",   7442, 21517,  2010, 21655,  2010, 21789,  2016, 21924,  1406),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00022020.art.event_list",   7506, 22020,  2014, 22156,  2006, 22294,  2009, 22423,  1477),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00022521.art.event_list",   7388, 22521,  2010, 22655,  2002, 22795,  2012, 22932,  1364),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00023025.art.event_list",   7230, 23025,  2011, 23163,  2001, 23309,  2010, 23446,  1208),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00023529.art.event_list",   7332, 23529,  2005, 23666,  2006, 23802,  2004, 23938,  1317),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00024031.art.event_list",   7361, 24031,  2013, 24170,  2017, 24310,  2018, 24444,  1313),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00024537.art.event_list",   7437, 24537,  2007, 24676,  2006, 24808,  2001, 24942,  1423),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00025040.art.event_list",   7490, 25040,  2010, 25173,  2011, 25303,  2007, 25441,  1462),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00025544.art.event_list",   7488, 25544,  2005, 25681,  2016, 25814,  2018, 25949,  1449),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00026050.art.event_list",   7501, 26050,  2011, 26194,  2008, 26325,  2009, 26458,  1473),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00026556.art.event_list",   7500, 26556,  2006, 26695,  2009, 26827,  2011, 26958,  1474),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00027059.art.event_list",   7518, 27059,  2011, 27193,  2002, 27322,  2012, 27460,  1493),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00027560.art.event_list",   7360, 27560,  2006, 27696,  2004, 27833,  2014, 27971,  1336),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00028064.art.event_list",   7274, 28064,  2002, 28195,  2009, 28333,  2002, 28482,  1261),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00028567.art.event_list",   7255, 28567,  2009, 28711,  2003, 28847,  2002, 28984,  1241),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00029070.art.event_list",   7365, 29070,  2002, 29209,  2013, 29346,  2012, 29481,  1338),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663671/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00029573.art.event_list",   6434, 29573,  2013, 29709,  2004, 29846,  2013, 29979,   404),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00030006.art.event_list",   7443, 30006,  2012, 30139,  2009, 30276,  2012, 30410,  1410),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00030512.art.event_list",   7314, 30512,  2001, 30649,  2001, 30784,  2006, 30924,  1306),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00031014.art.event_list",   7347, 31014,  2002, 31156,  2004, 31290,  2011, 31425,  1330),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00031518.art.event_list",   7253, 31518,  2011, 31658,  2003, 31797,  2007, 31936,  1232),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00032020.art.event_list",   7242, 32020,  2005, 32158,  2009, 32297,  2001, 32440,  1227),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00032524.art.event_list",   7320, 32524,  2009, 32664,  2005, 32801,  2012, 32938,  1294),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00033025.art.event_list",   7315, 33025,  2012, 33160,  2007, 33297,  2006, 33433,  1290),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00033526.art.event_list",   7323, 33526,  2016, 33671,  2006, 33807,  2002, 33941,  1299),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00034028.art.event_list",   7326, 34028,  2012, 34168,  2012, 34310,  2003, 34442,  1299),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00034530.art.event_list",   7296, 34530,  2002, 34670,  2002, 34806,  2017, 34944,  1275),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00035034.art.event_list",   7417, 35034,  2005, 35169,  2010, 35303,  2003, 35440,  1399),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00035536.art.event_list",   7321, 35536,  2002, 35674,  2013, 35812,  2005, 35952,  1301),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00036037.art.event_list",   7283, 36037,  2009, 36171,  2008, 36316,  2004, 36454,  1262),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00036539.art.event_list",   7288, 36539,  2012, 36677,  2011, 36819,  2013, 36957,  1252),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00037042.art.event_list",   7294, 37042,  2002, 37179,  2005, 37320,  2019, 37460,  1268),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00037543.art.event_list",   7512, 37543,  2003, 37677,  2014, 37812,  2008, 37947,  1487),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00038045.art.event_list",   7316, 38045,  2013, 38185,  2003, 38320,  2011, 38461,  1289),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00038548.art.event_list",   7494, 38548,  2009, 38685,  2002, 38815,  2014, 38948,  1469),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00039049.art.event_list",   7290, 39049,  2015, 39186,  2005, 39326,  2005, 39461,  1265),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35663963/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00039549.art.event_list",   6683, 39549,  2007, 39684,  2019, 39820,  2005, 39957,   652),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00040002.art.event_list",   7333, 40002,  2001, 40140,  2005, 40283,  2002, 40415,  1325),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00040505.art.event_list",   7508, 40505,  2013, 40638,  2004, 40772,  2001, 40907,  1490),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00041006.art.event_list",   7475, 41006,  2013, 41143,  2003, 41277,  2005, 41408,  1454),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00041508.art.event_list",   7154, 41508,  2003, 41648,  2003, 41790,  2002, 41930,  1146),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00042010.art.event_list",   7369, 42010,  2009, 42144,  2017, 42279,  2006, 42419,  1337),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00042512.art.event_list",   7334, 42512,  2014, 42649,  2013, 42791,  2004, 42925,  1303),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00043014.art.event_list",   7349, 43014,  2010, 43151,  2009, 43291,  2002, 43423,  1328),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00043515.art.event_list",   7346, 43515,  2011, 43648,  2004, 43787,  2015, 43923,  1316),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00044018.art.event_list",   7393, 44018,  2010, 44153,  2016, 44287,  2003, 44427,  1364),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00044519.art.event_list",   7238, 44519,  2012, 44661,  2001, 44800,  2003, 44938,  1222),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00045024.art.event_list",   7404, 45024,  2007, 45162,  2016, 45301,  2014, 45434,  1367),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00045525.art.event_list",   7465, 45525,  2014, 45658,  2012, 45796,  2005, 45931,  1434),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00046026.art.event_list",   7270, 46026,  2005, 46164,  2002, 46301,  2014, 46445,  1249),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00046527.art.event_list",   7393, 46527,  2011, 46665,  2003, 46798,  2005, 46932,  1374),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00047029.art.event_list",   7358, 47029,  2003, 47171,  2010, 47311,  2012, 47445,  1333),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00047532.art.event_list",   7252, 47532,  2016, 47667,  2014, 47806,  2001, 47949,  1221),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00048034.art.event_list",   7335, 48034,  2016, 48172,  2007, 48307,  2006, 48447,  1306),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00048536.art.event_list",   7429, 48536,  2003, 48675,  2019, 48810,  2003, 48950,  1404),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00049039.art.event_list",   7349, 49039,  2009, 49174,  2007, 49314,  2008, 49454,  1325),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/35664238/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00049539.art.event_list",   6804, 49539,  2006, 49674,  2009, 49813,  2019, 49952,   770),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00050002.art.event_list",   7569, 50002,  2017, 50134,  2007, 50266,  2015, 50404,  1530),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00050506.art.event_list",   7325, 50506,  2009, 50643,  2010, 50780,  2005, 50914,  1301),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00051007.art.event_list",   7343, 51007,  2012, 51142,  2010, 51275,  2007, 51418,  1314),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00051512.art.event_list",   7391, 51512,  2010, 51651,  2006, 51783,  2016, 51922,  1359),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00052015.art.event_list",   7479, 52015,  2014, 52148,  2011, 52284,  2004, 52423,  1450),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00052519.art.event_list",   7322, 52519,  2005, 52658,  2013, 52797,  2004, 52931,  1300),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00053019.art.event_list",   7330, 53019,  2019, 53155,  2007, 53291,  2011, 53433,  1293),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00053523.art.event_list",   7370, 53523,  2008, 53658,  2005, 53797,  2010, 53934,  1347),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00054027.art.event_list",   7437, 54027,  2011, 54169,  2005, 54302,  2011, 54440,  1410),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00054534.art.event_list",   7379, 54534,  2005, 54673,  2002, 54811,  2005, 54946,  1367),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00055036.art.event_list",   7406, 55036,  2004, 55172,  2006, 55315,  2005, 55447,  1391),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00055539.art.event_list",   7334, 55539,  2015, 55674,  2004, 55813,  2005, 55952,  1310),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00056039.art.event_list",   7254, 56039,  2008, 56177,  2017, 56317,  2007, 56455,  1222),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00056541.art.event_list",   7336, 56541,  2012, 56679,  2010, 56810,  2011, 56953,  1303),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00057044.art.event_list",   7482, 57044,  2014, 57180,  2012, 57319,  2003, 57454,  1453),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00057550.art.event_list",   7470, 57550,  2012, 57685,  2020, 57819,  2005, 57956,  1433),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00058056.art.event_list",   7460, 58056,  2010, 58188,  2010, 58325,  2001, 58463,  1439),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00058560.art.event_list",   7469, 58560,  2001, 58698,  2002, 58832,  2013, 58967,  1453),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00059063.art.event_list",   7307, 59063,  2009, 59198,  2008, 59336,  2015, 59480,  1275),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186551/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00059565.art.event_list",   6429, 59565,  2012, 59702,  2008, 59836,  2006, 59976,   403),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00060007.art.event_list",   7553, 60007,  2012, 60138,  2008, 60273,  2002, 60405,  1531),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00060512.art.event_list",   7472, 60512,  2001, 60647,  2016, 60785,  2009, 60914,  1446),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00061016.art.event_list",   7383, 61016,  2009, 61152,  2003, 61285,  2009, 61426,  1362),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00061520.art.event_list",   7310, 61520,  2016, 61657,  2010, 61795,  2011, 61935,  1273),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00062023.art.event_list",   7342, 62023,  2019, 62161,  2011, 62301,  2008, 62439,  1304),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00062527.art.event_list",   7497, 62527,  2019, 62662,  2011, 62795,  2008, 62931,  1459),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00063030.art.event_list",   7297, 63030,  2022, 63170,  2008, 63308,  2009, 63447,  1258),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00063531.art.event_list",   7400, 63531,  2011, 63672,  2006, 63803,  2010, 63938,  1373),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00064033.art.event_list",   7321, 64033,  2011, 64173,  2003, 64313,  2012, 64450,  1295),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00064537.art.event_list",   7520, 64537,  2008, 64670,  2009, 64809,  2009, 64942,  1494),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00065041.art.event_list",   7428, 65041,  2012, 65173,  2010, 65309,  2008, 65448,  1398),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00065544.art.event_list",   7375, 65544,  2004, 65681,  2014, 65812,  2006, 65951,  1351),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00066044.art.event_list",   7242, 66044,  2006, 66182,  2006, 66322,  2007, 66462,  1223),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00066544.art.event_list",   7342, 66544,  2014, 66680,  2016, 66818,  2010, 66957,  1302),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00067045.art.event_list",   7381, 67045,  2002, 67187,  2006, 67319,  2010, 67453,  1363),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00067547.art.event_list",   7313, 67547,  2012, 67685,  2009, 67828,  2004, 67963,  1288),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00068051.art.event_list",   7433, 68051,  2006, 68186,  2003, 68325,  2005, 68458,  1419),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00068555.art.event_list",   7445, 68555,  2001, 68688,  2009, 68821,  2005, 68957,  1430),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00069058.art.event_list",   7530, 69058,  2011, 69197,  2010, 69325,  2016, 69460,  1493),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186606/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00069560.art.event_list",   6468, 69560,  2012, 69694,  2013, 69833,  2001, 69971,   442),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00000/dig.mu2e.cry31s41b0.Mu2eII.001002_00070003.art.event_list",   7456, 70003,  2007, 70140,  2014, 70277,  2005, 70406,  1430),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00001/dig.mu2e.cry31s41b0.Mu2eII.001002_00070504.art.event_list",   7403, 70504,  2012, 70639,  2004, 70774,  2015, 70910,  1372),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00002/dig.mu2e.cry31s41b0.Mu2eII.001002_00071004.art.event_list",   7316, 71004,  2005, 71142,  2014, 71285,  2003, 71422,  1294),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00003/dig.mu2e.cry31s41b0.Mu2eII.001002_00071509.art.event_list",   7413, 71509,  2012, 71648,  2012, 71786,  2002, 71918,  1387),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00004/dig.mu2e.cry31s41b0.Mu2eII.001002_00072013.art.event_list",   7463, 72013,  2001, 72144,  2010, 72277,  2001, 72413,  1451),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00005/dig.mu2e.cry31s41b0.Mu2eII.001002_00072514.art.event_list",   7359, 72514,  2015, 72653,  2005, 72788,  2003, 72925,  1336),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00006/dig.mu2e.cry31s41b0.Mu2eII.001002_00073018.art.event_list",   7512, 73018,  2013, 73152,  2002, 73284,  2008, 73414,  1489),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00007/dig.mu2e.cry31s41b0.Mu2eII.001002_00073518.art.event_list",   7339, 73518,  2002, 73654,  2008, 73794,  2016, 73929,  1313),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00008/dig.mu2e.cry31s41b0.Mu2eII.001002_00074021.art.event_list",   7404, 74021,  2008, 74161,  2002, 74293,  2019, 74427,  1375),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00009/dig.mu2e.cry31s41b0.Mu2eII.001002_00074525.art.event_list",   7303, 74525,  2004, 74660,  2011, 74798,  2013, 74941,  1275),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00010/dig.mu2e.cry31s41b0.Mu2eII.001002_00075028.art.event_list",   7490, 75028,  2005, 75166,  2008, 75300,  2015, 75433,  1462),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00011/dig.mu2e.cry31s41b0.Mu2eII.001002_00075529.art.event_list",   7381, 75529,  2006, 75669,  2004, 75804,  2012, 75939,  1359),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00012/dig.mu2e.cry31s41b0.Mu2eII.001002_00076032.art.event_list",   7402, 76032,  2010, 76171,  2004, 76307,  2014, 76445,  1374),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00013/dig.mu2e.cry31s41b0.Mu2eII.001002_00076534.art.event_list",   7365, 76534,  2013, 76670,  2002, 76807,  2010, 76941,  1340),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00014/dig.mu2e.cry31s41b0.Mu2eII.001002_00077034.art.event_list",   7360, 77034,  2005, 77171,  2009, 77308,  2007, 77443,  1339),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00015/dig.mu2e.cry31s41b0.Mu2eII.001002_00077536.art.event_list",   7344, 77536,  2009, 77671,  2005, 77809,  2013, 77947,  1317),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00016/dig.mu2e.cry31s41b0.Mu2eII.001002_00078037.art.event_list",   7510, 78037,  2008, 78171,  2005, 78302,  2010, 78440,  1487),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00017/dig.mu2e.cry31s41b0.Mu2eII.001002_00078541.art.event_list",   7483, 78541,  2003, 78672,  2015, 78809,  2012, 78945,  1453),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00018/dig.mu2e.cry31s41b0.Mu2eII.001002_00079044.art.event_list",   7402, 79044,  2001, 79178,  2007, 79315,  2001, 79451,  1393),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cry31s00b0.s4_helix_filter/outstage/11186681/00/00019/dig.mu2e.cry31s41b0.Mu2eII.001002_00079546.art.event_list",   6739, 79546,  2001, 79678,  2002, 79816,  2019, 79955,   717)
]
#------------------------------------------------------------------------------
def fix_fcl(dir,index):

    # list_of_fcls = glob.glob('tmp/Mu2eII/fcl/cry31s41b0.s5_reco_stn_00/*00000000.fcl')
    list_of_fcls = glob.glob(dir+'/*.fcl')

    list_of_fcls.sort()

    print(list_of_fcls)

    first_run = 1002

    nfcls = len(list_of_fcls);

    for i in range(0,nfcls):
        # for fcl_file in list_of_fcls :
        fcl_file = list_of_fcls[i]

        print('------- fcl: ',fcl_file);

        lines = open(fcl_file).readlines()

        # print (lines)

        first_subrun = fix_data[i][2*index+2  ]
        max_events   = fix_data[i][2*index+2+1]

        f2 = open(fcl_file+'.new','w')
    
        for l in  lines:

            l1 = l.strip()

            w = l1.split()
            nw  = len(w)
            
            # print (w)

            if (nw>0):

                if (w[0] == 'source.maxEvents'):
                    l1 = 'source.maxEvents : %5i'%max_events
                elif (w[0] == 'source.skipEvents'):
                    l1 = 'source.firstRun: 1002'
                    f2.write(l1+'\n');
                    l1 = 'source.firstSubRun: %5i'%first_subrun
                elif (w[0] == 'services.TFileService.fileName'):
                    fn = fix_filename(w[2],first_run,first_subrun)
                    l1 = w[0]+'            : '+fn
                elif (w[0] == 'services.TFileService.fileName:'):
                    w[0] = w[0].replace(':','')
                    fn = fix_filename(w[1],first_run,first_subrun)
                    l1 = w[0]+'            : '+fn
                elif (w[0] == 'physics.filters.InitStntuple.histFileName'):
                    fn = fix_filename(w[2],first_run,first_subrun)
                    l1 = w[0]+' : '+fn
                elif (w[0] == 'outputs.defaultOutput.fileName'):
                    fn = fix_filename(w[2],first_run,first_subrun)
                    l1 = w[0]+'            : '+fn

            
            f2.write(l1+'\n');

        #
        f2.close()
    

#------------------------------------------------------------------------------
# main program, 
#------------------------------------------------------------------------------
if (__name__ == '__main__'):

    i = int(sys.argv[1])

    # dir = 'tmp/Mu2eII/fcl/cry31s41b0.s5_reco_stn_%02i'%i
    dir = 'tmp/Mu2eII/fcl/cry31s41b0.s5_reco_stn_%02i'%i
    print('dir=',dir,'i=',i)

    fix_fcl(dir,i)

    sys.exit(0);
