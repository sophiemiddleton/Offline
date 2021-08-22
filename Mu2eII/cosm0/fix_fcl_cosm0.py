#

import sys, glob



#------------------------------------------------------------------------------
def fix_filename(fn,first_run,first_subrun):

    w    = fn.split('.')
    w[4] = "%06i_%08i"%(first_run,first_subrun)
    fn1  = '.'.join(w);

    return fn1
    
fix_data = [
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00000/dig.mu2e.cosm0s41b0.Mu2eII.001002_00000000.art.event_list",   5352,     0,  2005,   187,  2006,   379,  1341),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00001/dig.mu2e.cosm0s41b0.Mu2eII.001002_00000502.art.event_list",   5354,   502,  2001,   696,  2016,   883,  1337),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00002/dig.mu2e.cosm0s41b0.Mu2eII.001002_00001003.art.event_list",   5325,  1003,  2001,  1189,  2006,  1377,  1318),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00003/dig.mu2e.cosm0s41b0.Mu2eII.001002_00001503.art.event_list",   5404,  1503,  2009,  1688,  2005,  1873,  1390),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00004/dig.mu2e.cosm0s41b0.Mu2eII.001002_00002006.art.event_list",   5363,  2006,  2006,  2202,  2005,  2383,  1352),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00005/dig.mu2e.cosm0s41b0.Mu2eII.001002_00002507.art.event_list",   5463,  2507,  2003,  2693,  2005,  2873,  1455),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00006/dig.mu2e.cosm0s41b0.Mu2eII.001002_00003007.art.event_list",   5253,  3007,  2001,  3195,  2004,  3391,  1248),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00007/dig.mu2e.cosm0s41b0.Mu2eII.001002_00003507.art.event_list",   5292,  3507,  2005,  3703,  2012,  3892,  1275),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00008/dig.mu2e.cosm0s41b0.Mu2eII.001002_00004007.art.event_list",   5490,  4007,  2003,  4190,  2003,  4375,  1484),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00009/dig.mu2e.cosm0s41b0.Mu2eII.001002_00004508.art.event_list",   5309,  4508,  2009,  4695,  2007,  4888,  1293),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00010/dig.mu2e.cosm0s41b0.Mu2eII.001002_00005009.art.event_list",   5297,  5009,  2011,  5194,  2001,  5385,  1285),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00011/dig.mu2e.cosm0s41b0.Mu2eII.001002_00005510.art.event_list",   5364,  5510,  2003,  5704,  2006,  5886,  1355),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00012/dig.mu2e.cosm0s41b0.Mu2eII.001002_00006011.art.event_list",   5333,  6011,  2009,  6205,  2002,  6387,  1322),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00013/dig.mu2e.cosm0s41b0.Mu2eII.001002_00006511.art.event_list",   5375,  6511,  2007,  6701,  2006,  6882,  1362),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00014/dig.mu2e.cosm0s41b0.Mu2eII.001002_00007011.art.event_list",   5491,  7011,  2003,  7190,  2007,  7377,  1481),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00015/dig.mu2e.cosm0s41b0.Mu2eII.001002_00007512.art.event_list",   5531,  7512,  2004,  7688,  2010,  7877,  1517),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00016/dig.mu2e.cosm0s41b0.Mu2eII.001002_00008013.art.event_list",   5370,  8013,  2002,  8203,  2009,  8387,  1359),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00017/dig.mu2e.cosm0s41b0.Mu2eII.001002_00008515.art.event_list",   5212,  8515,  2005,  8716,  2005,  8899,  1202),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00018/dig.mu2e.cosm0s41b0.Mu2eII.001002_00009016.art.event_list",   5522,  9016,  2003,  9190,  2001,  9377,  1518),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00019/dig.mu2e.cosm0s41b0.Mu2eII.001002_00009518.art.event_list",   5517,  9518,  2007,  9703,  2009,  9883,  1501),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00020/dig.mu2e.cosm0s41b0.Mu2eII.001002_00010018.art.event_list",   5369, 10018,  2004, 10214,  2006, 10396,  1359),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00021/dig.mu2e.cosm0s41b0.Mu2eII.001002_00010518.art.event_list",   5347, 10518,  2005, 10704,  2008, 10895,  1334),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00022/dig.mu2e.cosm0s41b0.Mu2eII.001002_00011019.art.event_list",   5339, 11019,  2009, 11212,  2008, 11394,  1322),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00023/dig.mu2e.cosm0s41b0.Mu2eII.001002_00011519.art.event_list",   5278, 11519,  2002, 11703,  2004, 11896,  1272),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00024/dig.mu2e.cosm0s41b0.Mu2eII.001002_00012019.art.event_list",   5397, 12019,  2001, 12209,  2002, 12394,  1394),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00025/dig.mu2e.cosm0s41b0.Mu2eII.001002_00012519.art.event_list",   5517, 12519,  2008, 12704,  2008, 12884,  1501),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00026/dig.mu2e.cosm0s41b0.Mu2eII.001002_00013019.art.event_list",   5381, 13019,  2001, 13208,  2005, 13399,  1375),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00027/dig.mu2e.cosm0s41b0.Mu2eII.001002_00013519.art.event_list",   5339, 13519,  2011, 13706,  2002, 13893,  1326),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00028/dig.mu2e.cosm0s41b0.Mu2eII.001002_00014019.art.event_list",   5440, 14019,  2004, 14201,  2003, 14388,  1433),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00029/dig.mu2e.cosm0s41b0.Mu2eII.001002_00014519.art.event_list",   5366, 14519,  2002, 14709,  2007, 14894,  1357),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00030/dig.mu2e.cosm0s41b0.Mu2eII.001002_00015019.art.event_list",   5300, 15019,  2006, 15215,  2009, 15398,  1285),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00031/dig.mu2e.cosm0s41b0.Mu2eII.001002_00015519.art.event_list",   5394, 15519,  2009, 15707,  2007, 15895,  1378),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00032/dig.mu2e.cosm0s41b0.Mu2eII.001002_00016020.art.event_list",   5309, 16020,  2001, 16214,  2001, 16401,  1307),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00033/dig.mu2e.cosm0s41b0.Mu2eII.001002_00016520.art.event_list",   5374, 16520,  2006, 16701,  2007, 16897,  1361),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00034/dig.mu2e.cosm0s41b0.Mu2eII.001002_00017021.art.event_list",   5352, 17021,  2002, 17208,  2009, 17394,  1341),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00035/dig.mu2e.cosm0s41b0.Mu2eII.001002_00017521.art.event_list",   5518, 17521,  2009, 17696,  2012, 17885,  1497),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00036/dig.mu2e.cosm0s41b0.Mu2eII.001002_00018021.art.event_list",   5310, 18021,  2004, 18215,  2012, 18404,  1294),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00037/dig.mu2e.cosm0s41b0.Mu2eII.001002_00018523.art.event_list",   5450, 18523,  2007, 18711,  2004, 18894,  1439),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00038/dig.mu2e.cosm0s41b0.Mu2eII.001002_00019026.art.event_list",   5384, 19026,  2002, 19217,  2009, 19400,  1373),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00039/dig.mu2e.cosm0s41b0.Mu2eII.001002_00019528.art.event_list",   5429, 19528,  2008, 19713,  2005, 19903,  1416),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00040/dig.mu2e.cosm0s41b0.Mu2eII.001002_00020029.art.event_list",   5352, 20029,  2006, 20212,  2003, 20402,  1343),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00041/dig.mu2e.cosm0s41b0.Mu2eII.001002_00020530.art.event_list",   5305, 20530,  2002, 20719,  2006, 20909,  1297),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00042/dig.mu2e.cosm0s41b0.Mu2eII.001002_00021030.art.event_list",   5402, 21030,  2004, 21221,  2015, 21408,  1383),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00043/dig.mu2e.cosm0s41b0.Mu2eII.001002_00021530.art.event_list",   5479, 21530,  2007, 21714,  2001, 21898,  1471),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00044/dig.mu2e.cosm0s41b0.Mu2eII.001002_00022031.art.event_list",   5527, 22031,  2009, 22215,  2018, 22396,  1500),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00045/dig.mu2e.cosm0s41b0.Mu2eII.001002_00022532.art.event_list",   5419, 22532,  2009, 22717,  2006, 22901,  1404),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00046/dig.mu2e.cosm0s41b0.Mu2eII.001002_00023033.art.event_list",   5329, 23033,  2007, 23223,  2002, 23407,  1320),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00047/dig.mu2e.cosm0s41b0.Mu2eII.001002_00023533.art.event_list",   5415, 23533,  2002, 23718,  2001, 23900,  1412),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00048/dig.mu2e.cosm0s41b0.Mu2eII.001002_00024034.art.event_list",   5225, 24034,  2004, 24231,  2004, 24421,  1217),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00049/dig.mu2e.cosm0s41b0.Mu2eII.001002_00024536.art.event_list",   5277, 24536,  2001, 24722,  2004, 24915,  1272),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00050/dig.mu2e.cosm0s41b0.Mu2eII.001002_00025039.art.event_list",   5267, 25039,  2003, 25228,  2010, 25425,  1254),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00051/dig.mu2e.cosm0s41b0.Mu2eII.001002_00025540.art.event_list",   5408, 25540,  2012, 25726,  2009, 25913,  1387),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00052/dig.mu2e.cosm0s41b0.Mu2eII.001002_00026041.art.event_list",   5432, 26041,  2005, 26219,  2008, 26405,  1419),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00053/dig.mu2e.cosm0s41b0.Mu2eII.001002_00026543.art.event_list",   5406, 26543,  2005, 26726,  2010, 26910,  1391),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00054/dig.mu2e.cosm0s41b0.Mu2eII.001002_00027043.art.event_list",   5459, 27043,  2008, 27227,  2006, 27409,  1445),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00055/dig.mu2e.cosm0s41b0.Mu2eII.001002_00027543.art.event_list",   5469, 27543,  2005, 27727,  2001, 27910,  1463),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00056/dig.mu2e.cosm0s41b0.Mu2eII.001002_00028045.art.event_list",   5428, 28045,  2014, 28229,  2005, 28414,  1409),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00057/dig.mu2e.cosm0s41b0.Mu2eII.001002_00028545.art.event_list",   5460, 28545,  2005, 28735,  2002, 28912,  1453),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00058/dig.mu2e.cosm0s41b0.Mu2eII.001002_00029048.art.event_list",   5234, 29048,  2020, 29244,  2009, 29435,  1205),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00059/dig.mu2e.cosm0s41b0.Mu2eII.001002_00029552.art.event_list",   5330, 29552,  2006, 29737,  2005, 29929,  1319),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00060/dig.mu2e.cosm0s41b0.Mu2eII.001002_00030052.art.event_list",   5318, 30052,  2003, 30239,  2009, 30428,  1306),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00061/dig.mu2e.cosm0s41b0.Mu2eII.001002_00030552.art.event_list",   5384, 30552,  2005, 30735,  2003, 30921,  1376),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00062/dig.mu2e.cosm0s41b0.Mu2eII.001002_00031052.art.event_list",   5401, 31052,  2001, 31238,  2006, 31422,  1394),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00063/dig.mu2e.cosm0s41b0.Mu2eII.001002_00031552.art.event_list",   5398, 31552,  2006, 31736,  2013, 31924,  1379),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00064/dig.mu2e.cosm0s41b0.Mu2eII.001002_00032052.art.event_list",   5318, 32052,  2004, 32245,  2004, 32435,  1310),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00065/dig.mu2e.cosm0s41b0.Mu2eII.001002_00032552.art.event_list",   5302, 32552,  2001, 32734,  2008, 32931,  1293),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00066/dig.mu2e.cosm0s41b0.Mu2eII.001002_00033052.art.event_list",   5286, 33052,  2002, 33241,  2009, 33437,  1275),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00067/dig.mu2e.cosm0s41b0.Mu2eII.001002_00033555.art.event_list",   5417, 33555,  2012, 33746,  2003, 33926,  1402),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00068/dig.mu2e.cosm0s41b0.Mu2eII.001002_00034055.art.event_list",   5533, 34055,  2009, 34235,  2002, 34418,  1522),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00069/dig.mu2e.cosm0s41b0.Mu2eII.001002_00034555.art.event_list",   5399, 34555,  2006, 34742,  2011, 34928,  1382),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00070/dig.mu2e.cosm0s41b0.Mu2eII.001002_00035055.art.event_list",   5335, 35055,  2014, 35248,  2003, 35438,  1318),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00071/dig.mu2e.cosm0s41b0.Mu2eII.001002_00035555.art.event_list",   5555, 35555,  2007, 35732,  2012, 35910,  1536),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00072/dig.mu2e.cosm0s41b0.Mu2eII.001002_00036055.art.event_list",   5449, 36055,  2010, 36240,  2010, 36426,  1429),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00073/dig.mu2e.cosm0s41b0.Mu2eII.001002_00036556.art.event_list",   5230, 36556,  2003, 36741,  2006, 36930,  1221),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00074/dig.mu2e.cosm0s41b0.Mu2eII.001002_00037057.art.event_list",   5450, 37057,  2001, 37241,  2005, 37429,  1444),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00075/dig.mu2e.cosm0s41b0.Mu2eII.001002_00037560.art.event_list",   5250, 37560,  2008, 37748,  2007, 37943,  1235),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00076/dig.mu2e.cosm0s41b0.Mu2eII.001002_00038061.art.event_list",   5311, 38061,  2004, 38245,  2008, 38437,  1299),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00077/dig.mu2e.cosm0s41b0.Mu2eII.001002_00038563.art.event_list",   5330, 38563,  2008, 38756,  2003, 38942,  1319),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00078/dig.mu2e.cosm0s41b0.Mu2eII.001002_00039065.art.event_list",   5427, 39065,  2009, 39249,  2011, 39430,  1407),
("/pnfs/mu2e/persistent/users/mu2epro/workflow/Mu2eII.cosm0s00b0.s4_helix_filter/outstage/35162167/00/00079/dig.mu2e.cosm0s41b0.Mu2eII.001002_00039566.art.event_list",   4609, 39566,  2010, 39754,  2005, 39944,   594)
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
    dir = 'tmp/Mu2eII/fcl/cosm0s41b0.s5_reco_stn_%02i'%i
    print('dir=',dir,'i=',i)

    fix_fcl(dir,i)

    sys.exit(0);
