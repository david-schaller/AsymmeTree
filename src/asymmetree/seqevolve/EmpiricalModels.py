# -*- coding: utf-8 -*-

import numpy as np


def _WAG_aa():
    
    S = np.array([[      0.0,  0.551571,  0.509848,  0.738998,   1.02704,  0.908598,   1.58285,   1.41672,  0.316954,  0.193335,  0.397915,  0.906265,  0.893496,  0.210494,   1.43855,   3.37079,   2.12111,  0.113133,  0.240735,   2.00601],
                  [ 0.551571,       0.0,  0.635346,  0.147304,  0.528191,    3.0355,  0.439157,  0.584665,   2.13715,  0.186979,  0.497671,   5.35142,  0.683162,  0.102711,  0.679489,   1.22419,  0.554413,   1.16392,  0.381533,  0.251849],
                  [ 0.509848,  0.635346,       0.0,   5.42942,  0.265256,   1.54364,  0.947198,   1.12556,   3.95629,  0.554236,  0.131528,   3.01201,  0.198221, 0.0961621,  0.195081,   3.97423,   2.03006, 0.0719167,     1.086,  0.196246],
                  [ 0.738998,  0.147304,   5.42942,       0.0, 0.0302949,  0.616783,   6.17416,  0.865584,  0.930676,  0.039437, 0.0848047,  0.479855,  0.103754, 0.0467304,  0.423984,   1.07176,  0.374866,  0.129767,  0.325711,  0.152335],
                  [  1.02704,  0.528191,  0.265256, 0.0302949,       0.0, 0.0988179,  0.021352,  0.306674,  0.248972,  0.170135,  0.384287, 0.0740339,  0.390482,   0.39802,  0.109404,   1.40766,  0.512984,   0.71707,  0.543833,   1.00214],
                  [ 0.908598,    3.0355,   1.54364,  0.616783, 0.0988179,       0.0,   5.46947,  0.330052,   4.29411,  0.113917,  0.869489,    3.8949,   1.54526, 0.0999208,  0.933372,   1.02887,  0.857928,  0.215737,   0.22771,  0.301281],
                  [  1.58285,  0.439157,  0.947198,   6.17416,  0.021352,   5.46947,       0.0,  0.567717,  0.570025,  0.127395,  0.154263,   2.58443,  0.315124, 0.0811339,  0.682355,  0.704939,  0.822765,  0.156557,  0.196303,  0.588731],
                  [  1.41672,  0.584665,   1.12556,  0.865584,  0.306674,  0.330052,  0.567717,       0.0,   0.24941, 0.0304501, 0.0613037,  0.373558,    0.1741,  0.049931,   0.24357,   1.34182,  0.225833,  0.336983,  0.103604,  0.187247],
                  [ 0.316954,   2.13715,   3.95629,  0.930676,  0.248972,   4.29411,  0.570025,   0.24941,       0.0,   0.13819,  0.499462,  0.890432,  0.404141,  0.679371,  0.696198,  0.740169,  0.473307,  0.262569,   3.87344,  0.118358],
                  [ 0.193335,  0.186979,  0.554236,  0.039437,  0.170135,  0.113917,  0.127395, 0.0304501,   0.13819,       0.0,   3.17097,  0.323832,   4.25746,   1.05947, 0.0999288,   0.31944,   1.45816,  0.212483,   0.42017,    7.8213],
                  [ 0.397915,  0.497671,  0.131528, 0.0848047,  0.384287,  0.869489,  0.154263, 0.0613037,  0.499462,   3.17097,       0.0,  0.257555,   4.85402,   2.11517,  0.415844,  0.344739,  0.326622,  0.665309,  0.398618,   1.80034],
                  [ 0.906265,   5.35142,   3.01201,  0.479855, 0.0740339,    3.8949,   2.58443,  0.373558,  0.890432,  0.323832,  0.257555,       0.0,  0.934276,  0.088836,  0.556896,   0.96713,   1.38698,  0.137505,  0.133264,  0.305434],
                  [ 0.893496,  0.683162,  0.198221,  0.103754,  0.390482,   1.54526,  0.315124,    0.1741,  0.404141,   4.25746,   4.85402,  0.934276,       0.0,   1.19063,  0.171329,  0.493905,   1.51612,  0.515706,  0.428437,   2.05845],
                  [ 0.210494,  0.102711, 0.0961621, 0.0467304,   0.39802, 0.0999208, 0.0811339,  0.049931,  0.679371,   1.05947,   2.11517,  0.088836,   1.19063,       0.0,  0.161444,  0.545931,  0.171903,   1.52964,   6.45428,  0.649892],
                  [  1.43855,  0.679489,  0.195081,  0.423984,  0.109404,  0.933372,  0.682355,   0.24357,  0.696198, 0.0999288,  0.415844,  0.556896,  0.171329,  0.161444,       0.0,   1.61328,  0.795384,  0.139405,  0.216046,  0.314887],
                  [  3.37079,   1.22419,   3.97423,   1.07176,   1.40766,   1.02887,  0.704939,   1.34182,  0.740169,   0.31944,  0.344739,   0.96713,  0.493905,  0.545931,   1.61328,       0.0,   4.37802,  0.523742,  0.786993,  0.232739],
                  [  2.12111,  0.554413,   2.03006,  0.374866,  0.512984,  0.857928,  0.822765,  0.225833,  0.473307,   1.45816,  0.326622,   1.38698,   1.51612,  0.171903,  0.795384,   4.37802,       0.0,  0.110864,  0.291148,   1.38823],
                  [ 0.113133,   1.16392, 0.0719167,  0.129767,   0.71707,  0.215737,  0.156557,  0.336983,  0.262569,  0.212483,  0.665309,  0.137505,  0.515706,   1.52964,  0.139405,  0.523742,  0.110864,       0.0,   2.48539,  0.365369],
                  [ 0.240735,  0.381533,     1.086,  0.325711,  0.543833,   0.22771,  0.196303,  0.103604,   3.87344,   0.42017,  0.398618,  0.133264,  0.428437,   6.45428,  0.216046,  0.786993,  0.291148,   2.48539,       0.0,   0.31473],
                  [  2.00601,  0.251849,  0.196246,  0.152335,   1.00214,  0.301281,  0.588731,  0.187247,  0.118358,    7.8213,   1.80034,  0.305434,   2.05845,  0.649892,  0.314887,  0.232739,   1.38823,  0.365369,   0.31473,       0.0]])
    
    freqs = np.array([0.0866279,  0.043972, 0.0390894, 0.0570451, 0.0193078, 0.0367281, 0.0580589, 0.0832518, 0.0244313,  0.048466,  0.086209, 0.0620286, 0.0195027, 0.0384319, 0.0457631, 0.0695179, 0.0610127, 0.0143859, 0.0352742, 0.0708956])
    
    return S, freqs


def _JTT_aa():
    
    S = np.array([[  0.0,  58.0,  54.0,  81.0,  56.0,  57.0, 105.0, 179.0,  27.0,  36.0,  30.0,  35.0,  54.0,  15.0, 194.0, 378.0, 475.0,   9.0,  11.0, 298.0],
                  [ 58.0,   0.0,  45.0,  16.0, 113.0, 310.0,  29.0, 137.0, 328.0,  22.0,  38.0, 646.0,  44.0,   5.0,  74.0, 101.0,  64.0, 126.0,  20.0,  17.0],
                  [ 54.0,  45.0,   0.0, 528.0,  34.0,  86.0,  58.0,  81.0, 391.0,  47.0,  12.0, 263.0,  30.0,  10.0,  15.0, 503.0, 232.0,   8.0,  70.0,  16.0],
                  [ 81.0,  16.0, 528.0,   0.0,  10.0,  49.0, 767.0, 130.0, 112.0,  11.0,   7.0,  26.0,  15.0,   4.0,  15.0,  59.0,  38.0,   4.0,  46.0,  31.0],
                  [ 56.0, 113.0,  34.0,  10.0,   0.0,   9.0,   5.0,  59.0,  69.0,  17.0,  23.0,   7.0,  31.0,  78.0,  14.0, 223.0,  42.0, 115.0, 209.0,  62.0],
                  [ 57.0, 310.0,  86.0,  49.0,   9.0,   0.0, 323.0,  26.0, 597.0,   9.0,  72.0, 292.0,  43.0,   4.0, 164.0,  53.0,  51.0,  18.0,  24.0,  20.0],
                  [105.0,  29.0,  58.0, 767.0,   5.0, 323.0,   0.0, 119.0,  26.0,  12.0,   9.0, 181.0,  18.0,   5.0,  18.0,  30.0,  32.0,  10.0,   7.0,  45.0],
                  [179.0, 137.0,  81.0, 130.0,  59.0,  26.0, 119.0,   0.0,  23.0,   6.0,   6.0,  27.0,  14.0,   5.0,  24.0, 201.0,  33.0,  55.0,   8.0,  47.0],
                  [ 27.0, 328.0, 391.0, 112.0,  69.0, 597.0,  26.0,  23.0,   0.0,  16.0,  56.0,  45.0,  33.0,  40.0, 115.0,  73.0,  46.0,   8.0, 573.0,  11.0],
                  [ 36.0,  22.0,  47.0,  11.0,  17.0,   9.0,  12.0,   6.0,  16.0,   0.0, 229.0,  21.0, 479.0,  89.0,  10.0,  40.0, 245.0,   9.0,  32.0, 961.0],
                  [ 30.0,  38.0,  12.0,   7.0,  23.0,  72.0,   9.0,   6.0,  56.0, 229.0,   0.0,  14.0, 388.0, 248.0, 102.0,  59.0,  25.0,  52.0,  24.0, 180.0],
                  [ 35.0, 646.0, 263.0,  26.0,   7.0, 292.0, 181.0,  27.0,  45.0,  21.0,  14.0,   0.0,  65.0,   4.0,  21.0,  47.0, 103.0,  10.0,   8.0,  14.0],
                  [ 54.0,  44.0,  30.0,  15.0,  31.0,  43.0,  18.0,  14.0,  33.0, 479.0, 388.0,  65.0,   0.0,  43.0,  16.0,  29.0, 226.0,  24.0,  18.0, 323.0],
                  [ 15.0,   5.0,  10.0,   4.0,  78.0,   4.0,   5.0,   5.0,  40.0,  89.0, 248.0,   4.0,  43.0,   0.0,  17.0,  92.0,  12.0,  53.0, 536.0,  62.0],
                  [194.0,  74.0,  15.0,  15.0,  14.0, 164.0,  18.0,  24.0, 115.0,  10.0, 102.0,  21.0,  16.0,  17.0,   0.0, 285.0, 118.0,   6.0,  10.0,  23.0],
                  [378.0, 101.0, 503.0,  59.0, 223.0,  53.0,  30.0, 201.0,  73.0,  40.0,  59.0,  47.0,  29.0,  92.0, 285.0,   0.0, 477.0,  35.0,  63.0,  38.0],
                  [475.0,  64.0, 232.0,  38.0,  42.0,  51.0,  32.0,  33.0,  46.0, 245.0,  25.0, 103.0, 226.0,  12.0, 118.0, 477.0,   0.0,  12.0,  21.0, 112.0],
                  [  9.0, 126.0,   8.0,   4.0, 115.0,  18.0,  10.0,  55.0,   8.0,   9.0,  52.0,  10.0,  24.0,  53.0,   6.0,  35.0,  12.0,   0.0,  71.0,  25.0],
                  [ 11.0,  20.0,  70.0,  46.0, 209.0,  24.0,   7.0,   8.0, 573.0,  32.0,  24.0,   8.0,  18.0, 536.0,  10.0,  63.0,  21.0,  71.0,   0.0,  16.0],
                  [298.0,  17.0,  16.0,  31.0,  62.0,  20.0,  45.0,  47.0,  11.0, 961.0, 180.0,  14.0, 323.0,  62.0,  23.0,  38.0, 112.0,  25.0,  16.0,   0.0]])
    
    freqs = np.array([0.076748, 0.051691, 0.042645, 0.051544, 0.019803, 0.040752,  0.06183, 0.073152, 0.022944, 0.053761, 0.091904, 0.058676, 0.023826, 0.040126, 0.050901, 0.068765, 0.058565, 0.014261, 0.032102, 0.066005])
    
    return S, freqs


def _BLOSUM62_aa():
    
    S = np.array([[           0.0, 0.735790389698, 0.485391055466, 0.543161820899,  1.45999531047, 1.199705704602,   1.1709490428,  1.95588357496, 0.716241444998, 0.605899003687, 0.800016530518, 1.295201266783, 1.253758266664, 0.492964679748, 1.173275900924, 4.325092687057, 1.729178019485, 0.465839367725, 0.718206697586, 2.187774522005],
                  [0.735790389698,            0.0, 1.297446705134, 0.500964408555, 0.227826574209, 3.020833610064,  1.36057419042, 0.418763308518, 1.456141166336, 0.232036445142, 0.622711669692, 5.411115141489, 0.983692987457, 0.371644693209, 0.448133661718,  1.12278310421, 0.914665954563, 0.426382310122, 0.720517441216, 0.438388343772],
                  [0.485391055466, 1.297446705134,            0.0, 3.180100048216, 0.397358949897, 1.839216146992,  1.24048850864, 1.355872344485, 2.414501434208, 0.283017326278, 0.211888159615, 1.593137043457, 0.648441278787, 0.354861249223, 0.494887043702, 2.904101656456, 1.898173634533, 0.191482046247, 0.538222519037, 0.312858797993],
                  [0.543161820899, 0.500964408555, 3.180100048216,            0.0, 0.240836614802, 1.190945703396, 3.761625208368, 0.798473248968, 0.778142664022, 0.418555732462, 0.218131577594, 1.032447924952, 0.222621897958, 0.281730694207, 0.730628272998, 1.582754142065, 0.934187509431, 0.145345046279, 0.261422208965, 0.258129289418],
                  [ 1.45999531047, 0.227826574209, 0.397358949897, 0.240836614802,            0.0,  0.32980150463, 0.140748891814, 0.418203192284, 0.354058109831, 0.774894022794, 0.831842640142, 0.285078800906,  0.76768882348, 0.441337471187, 0.356008498769, 1.197188415094, 1.119831358516, 0.527664418872, 0.470237733696, 1.116352478606],
                  [1.199705704602, 3.020833610064, 1.839216146992, 1.190945703396,  0.32980150463,            0.0, 5.528919177928, 0.609846305383,  2.43534113114, 0.236202451204, 0.580737093181, 3.945277674515, 2.494896077113,  0.14435695975, 0.858570575674, 1.934870924596, 1.277480294596, 0.758653808642,  0.95898974285, 0.530785790125],
                  [  1.1709490428,  1.36057419042,  1.24048850864, 3.761625208368, 0.140748891814, 5.528919177928,            0.0, 0.423579992176, 1.626891056982, 0.186848046932, 0.372625175087, 2.802427151679,  0.55541539747, 0.291409084165, 0.926563934846, 1.769893238937, 1.071097236007, 0.407635648938, 0.596719300346, 0.524253846338],
                  [ 1.95588357496, 0.418763308518, 1.355872344485, 0.798473248968, 0.418203192284, 0.609846305383, 0.423579992176,            0.0, 0.539859124954, 0.189296292376, 0.217721159236, 0.752042440303, 0.459436173579, 0.368166464453, 0.504086599527, 1.509326253224, 0.641436011405, 0.508358924638, 0.308055737035,  0.25334079019],
                  [0.716241444998, 1.456141166336, 2.414501434208, 0.778142664022, 0.354058109831,  2.43534113114, 1.626891056982, 0.539859124954,            0.0, 0.252718447885, 0.348072209797, 1.022507035889, 0.984311525359, 0.714533703928, 0.527007339151,  1.11702976291, 0.585407090225,  0.30124860078, 4.218953969389,  0.20155597175],
                  [0.605899003687, 0.232036445142, 0.283017326278, 0.418555732462, 0.774894022794, 0.236202451204, 0.186848046932, 0.189296292376, 0.252718447885,            0.0, 3.890963773304, 0.406193586642, 3.364797763104, 1.517359325954, 0.388355409206,  0.35754441246,  1.17909119726,  0.34198578754, 0.674617093228, 8.311839405458],
                  [0.800016530518, 0.622711669692, 0.211888159615, 0.218131577594, 0.831842640142, 0.580737093181, 0.372625175087, 0.217721159236, 0.348072209797, 3.890963773304,            0.0, 0.445570274261, 6.030559379572, 2.064839703237, 0.374555687471, 0.352969184527, 0.915259857694,   0.6914746346, 0.811245856323, 2.231405688913],
                  [1.295201266783, 5.411115141489, 1.593137043457, 1.032447924952, 0.285078800906, 3.945277674515, 2.802427151679, 0.752042440303, 1.022507035889, 0.406193586642, 0.445570274261,            0.0, 1.073061184332, 0.266924750511, 1.047383450722, 1.752165917819, 1.303875200799, 0.332243040634,   0.7179934869, 0.498138475304],
                  [1.253758266664, 0.983692987457, 0.648441278787, 0.222621897958,  0.76768882348, 2.494896077113,  0.55541539747, 0.459436173579, 0.984311525359, 3.364797763104, 6.030559379572, 1.073061184332,            0.0,  1.77385516883, 0.454123625103, 0.918723415746, 1.488548053722, 0.888101098152, 0.951682162246, 2.575850755315],
                  [0.492964679748, 0.371644693209, 0.354861249223, 0.281730694207, 0.441337471187,  0.14435695975, 0.291409084165, 0.368166464453, 0.714533703928, 1.517359325954, 2.064839703237, 0.266924750511,  1.77385516883,            0.0, 0.233597909629, 0.540027644824, 0.488206118793, 2.074324893497, 6.747260430801, 0.838119610178],
                  [1.173275900924, 0.448133661718, 0.494887043702, 0.730628272998, 0.356008498769, 0.858570575674, 0.926563934846, 0.504086599527, 0.527007339151, 0.388355409206, 0.374555687471, 1.047383450722, 0.454123625103, 0.233597909629,            0.0, 1.169129577716, 1.005451683149, 0.252214830027, 0.369405319355, 0.496908410676],
                  [4.325092687057,  1.12278310421, 2.904101656456, 1.582754142065, 1.197188415094, 1.934870924596, 1.769893238937, 1.509326253224,  1.11702976291,  0.35754441246, 0.352969184527, 1.752165917819, 0.918723415746, 0.540027644824, 1.169129577716,            0.0,  5.15155629227, 0.387925622098, 0.796751520761, 0.561925457442],
                  [1.729178019485, 0.914665954563, 1.898173634533, 0.934187509431, 1.119831358516, 1.277480294596, 1.071097236007, 0.641436011405, 0.585407090225,  1.17909119726, 0.915259857694, 1.303875200799, 1.488548053722, 0.488206118793, 1.005451683149,  5.15155629227,            0.0, 0.513128126891, 0.801010243199, 2.253074051176],
                  [0.465839367725, 0.426382310122, 0.191482046247, 0.145345046279, 0.527664418872, 0.758653808642, 0.407635648938, 0.508358924638,  0.30124860078,  0.34198578754,   0.6914746346, 0.332243040634, 0.888101098152, 2.074324893497, 0.252214830027, 0.387925622098, 0.513128126891,            0.0, 4.054419006558, 0.266508731426],
                  [0.718206697586, 0.720517441216, 0.538222519037, 0.261422208965, 0.470237733696,  0.95898974285, 0.596719300346, 0.308055737035, 4.218953969389, 0.674617093228, 0.811245856323,   0.7179934869, 0.951682162246, 6.747260430801, 0.369405319355, 0.796751520761, 0.801010243199, 4.054419006558,            0.0,            1.0],
                  [2.187774522005, 0.438388343772, 0.312858797993, 0.258129289418, 1.116352478606, 0.530785790125, 0.524253846338,  0.25334079019,  0.20155597175, 8.311839405458, 2.231405688913, 0.498138475304, 2.575850755315, 0.838119610178, 0.496908410676, 0.561925457442, 2.253074051176, 0.266508731426,            1.0,            0.0]])
    
    freqs = np.array([0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074, 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057, 0.051, 0.013, 0.032, 0.073])
    
    return S, freqs


def _LG_aa():
    
    S = np.array([[      0.0,  0.425093,  0.276818,  0.395144,  2.489084,  0.969894,  1.038545,   2.06604,  0.358858,   0.14983,  0.395337,  0.536518,  1.124035,  0.253701,  1.177651,  4.727182,  2.139501,  0.180717,  0.218959,   2.54787],
                  [ 0.425093,       0.0,  0.751878,  0.123954,  0.534551,  2.807908,   0.36397,  0.390192,  2.426601,  0.126991,  0.301848,  6.326067,  0.484133,  0.052722,  0.332533,  0.858151,  0.578987,  0.593607,   0.31444,  0.170887],
                  [ 0.276818,  0.751878,       0.0,  5.076149,  0.528768,  1.695752,  0.541712,  1.437645,  4.509238,  0.191503,  0.068427,  2.145078,  0.371004,  0.089525,  0.161787,  4.008358,  2.000679,  0.045376,  0.612025,  0.083688],
                  [ 0.395144,  0.123954,  5.076149,       0.0,  0.062556,  0.523386,   5.24387,  0.844926,  0.927114,   0.01069,  0.015076,  0.282959,  0.025548,  0.017416,  0.394456,  1.240275,   0.42586,   0.02989,  0.135107,  0.037967],
                  [ 2.489084,  0.534551,  0.528768,  0.062556,       0.0,  0.084808,  0.003499,  0.569265,  0.640543,  0.320627,  0.594007,  0.013266,   0.89368,  1.105251,  0.075382,  2.784478,   1.14348,  0.670128,  1.165532,  1.959291],
                  [ 0.969894,  2.807908,  1.695752,  0.523386,  0.084808,       0.0,  4.128591,  0.267959,  4.813505,  0.072854,  0.582457,  3.234294,  1.672569,  0.035855,  0.624294,  1.223828,  1.080136,  0.236199,  0.257336,  0.210332],
                  [ 1.038545,   0.36397,  0.541712,   5.24387,  0.003499,  4.128591,       0.0,  0.348847,  0.423881,  0.044265,  0.069673,  1.807177,  0.173735,  0.018811,  0.419409,  0.611973,  0.604545,  0.077852,  0.120037,  0.245034],
                  [  2.06604,  0.390192,  1.437645,  0.844926,  0.569265,  0.267959,  0.348847,       0.0,  0.311484,  0.008705,  0.044261,  0.296636,  0.139538,  0.089586,  0.196961,   1.73999,  0.129836,  0.268491,  0.054679,  0.076701],
                  [ 0.358858,  2.426601,  4.509238,  0.927114,  0.640543,  4.813505,  0.423881,  0.311484,       0.0,  0.108882,  0.366317,  0.697264,  0.442472,  0.682139,  0.508851,  0.990012,  0.584262,  0.597054,  5.306834,  0.119013],
                  [  0.14983,  0.126991,  0.191503,   0.01069,  0.320627,  0.072854,  0.044265,  0.008705,  0.108882,       0.0,  4.145067,  0.159069,  4.273607,  1.112727,  0.078281,  0.064105,  1.033739,   0.11166,  0.232523, 10.649107],
                  [ 0.395337,  0.301848,  0.068427,  0.015076,  0.594007,  0.582457,  0.069673,  0.044261,  0.366317,  4.145067,       0.0,    0.1375,  6.312358,  2.592692,   0.24906,  0.182287,  0.302936,  0.619632,  0.299648,  1.702745],
                  [ 0.536518,  6.326067,  2.145078,  0.282959,  0.013266,  3.234294,  1.807177,  0.296636,  0.697264,  0.159069,    0.1375,       0.0,  0.656604,  0.023918,  0.390322,  0.748683,  1.136863,  0.049906,  0.131932,  0.185202],
                  [ 1.124035,  0.484133,  0.371004,  0.025548,   0.89368,  1.672569,  0.173735,  0.139538,  0.442472,  4.273607,  6.312358,  0.656604,       0.0,  1.798853,  0.099849,   0.34696,  2.020366,  0.696175,  0.481306,  1.898718],
                  [ 0.253701,  0.052722,  0.089525,  0.017416,  1.105251,  0.035855,  0.018811,  0.089586,  0.682139,  1.112727,  2.592692,  0.023918,  1.798853,       0.0,  0.094464,  0.361819,  0.165001,  2.457121,  7.803902,  0.654683],
                  [ 1.177651,  0.332533,  0.161787,  0.394456,  0.075382,  0.624294,  0.419409,  0.196961,  0.508851,  0.078281,   0.24906,  0.390322,  0.099849,  0.094464,       0.0,  1.338132,  0.571468,  0.095131,  0.089613,  0.296501],
                  [ 4.727182,  0.858151,  4.008358,  1.240275,  2.784478,  1.223828,  0.611973,   1.73999,  0.990012,  0.064105,  0.182287,  0.748683,   0.34696,  0.361819,  1.338132,       0.0,  6.472279,  0.248862,  0.400547,  0.098369],
                  [ 2.139501,  0.578987,  2.000679,   0.42586,   1.14348,  1.080136,  0.604545,  0.129836,  0.584262,  1.033739,  0.302936,  1.136863,  2.020366,  0.165001,  0.571468,  6.472279,       0.0,  0.140825,  0.245841,  2.188158],
                  [ 0.180717,  0.593607,  0.045376,   0.02989,  0.670128,  0.236199,  0.077852,  0.268491,  0.597054,   0.11166,  0.619632,  0.049906,  0.696175,  2.457121,  0.095131,  0.248862,  0.140825,       0.0,  3.151815,   0.18951],
                  [ 0.218959,   0.31444,  0.612025,  0.135107,  1.165532,  0.257336,  0.120037,  0.054679,  5.306834,  0.232523,  0.299648,  0.131932,  0.481306,  7.803902,  0.089613,  0.400547,  0.245841,  3.151815,       0.0,  0.249313],
                  [  2.54787,  0.170887,  0.083688,  0.037967,  1.959291,  0.210332,  0.245034,  0.076701,  0.119013, 10.649107,  1.702745,  0.185202,  1.898718,  0.654683,  0.296501,  0.098369,  2.188158,   0.18951,  0.249313,       0.0]])
    
    freqs = np.array([0.079066, 0.055941, 0.041977, 0.053052, 0.012937, 0.040767, 0.071586, 0.057337, 0.022355, 0.062157, 0.099081,   0.0646, 0.022951, 0.042302,  0.04404, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147])
    
    return S, freqs


def _DAYHOFF_aa():
    
    S = np.array([[   0.0,   27.0,   98.0,  120.0,   36.0,   89.0,  198.0,  240.0,   23.0,   65.0,   41.0,   26.0,   72.0,   18.0,  250.0,  409.0,  371.0,    0.0,   24.0,  208.0],
                  [  27.0,    0.0,   32.0,    0.0,   23.0,  246.0,    1.0,    9.0,  240.0,   64.0,   15.0,  464.0,   90.0,   14.0,  103.0,  154.0,   26.0,  201.0,    8.0,   24.0],
                  [  98.0,   32.0,    0.0,  905.0,    0.0,  103.0,  148.0,  139.0,  535.0,   77.0,   34.0,  318.0,    1.0,   14.0,   42.0,  495.0,  229.0,   23.0,   95.0,   15.0],
                  [ 120.0,    0.0,  905.0,    0.0,    0.0,  134.0, 1153.0,  125.0,   86.0,   24.0,    0.0,   71.0,    0.0,    0.0,   13.0,   95.0,   66.0,    0.0,    0.0,   18.0],
                  [  36.0,   23.0,    0.0,    0.0,    0.0,    0.0,    0.0,   11.0,   28.0,   44.0,    0.0,    0.0,    0.0,    0.0,   19.0,  161.0,   16.0,    0.0,   96.0,   49.0],
                  [  89.0,  246.0,  103.0,  134.0,    0.0,    0.0,  716.0,   28.0,  606.0,   18.0,   73.0,  153.0,  114.0,    0.0,  153.0,   56.0,   53.0,    0.0,    0.0,   35.0],
                  [ 198.0,    1.0,  148.0, 1153.0,    0.0,  716.0,    0.0,   81.0,   43.0,   61.0,   11.0,   83.0,   30.0,    0.0,   51.0,   79.0,   34.0,    0.0,   22.0,   37.0],
                  [ 240.0,    9.0,  139.0,  125.0,   11.0,   28.0,   81.0,    0.0,   10.0,    0.0,    7.0,   27.0,   17.0,   15.0,   34.0,  234.0,   30.0,    0.0,    0.0,   54.0],
                  [  23.0,  240.0,  535.0,   86.0,   28.0,  606.0,   43.0,   10.0,    0.0,    7.0,   44.0,   26.0,    0.0,   48.0,   94.0,   35.0,   22.0,   27.0,  127.0,   44.0],
                  [  65.0,   64.0,   77.0,   24.0,   44.0,   18.0,   61.0,    0.0,    7.0,    0.0,  257.0,   46.0,  336.0,  196.0,   12.0,   24.0,  192.0,    0.0,   37.0,  889.0],
                  [  41.0,   15.0,   34.0,    0.0,    0.0,   73.0,   11.0,    7.0,   44.0,  257.0,    0.0,   18.0,  527.0,  157.0,   32.0,   17.0,   33.0,   46.0,   28.0,  175.0],
                  [  26.0,  464.0,  318.0,   71.0,    0.0,  153.0,   83.0,   27.0,   26.0,   46.0,   18.0,    0.0,  243.0,    0.0,   33.0,   96.0,  136.0,    0.0,   13.0,   10.0],
                  [  72.0,   90.0,    1.0,    0.0,    0.0,  114.0,   30.0,   17.0,    0.0,  336.0,  527.0,  243.0,    0.0,   92.0,   17.0,   62.0,  104.0,    0.0,    0.0,  258.0],
                  [  18.0,   14.0,   14.0,    0.0,    0.0,    0.0,    0.0,   15.0,   48.0,  196.0,  157.0,    0.0,   92.0,    0.0,   11.0,   46.0,   13.0,   76.0,  698.0,   12.0],
                  [ 250.0,  103.0,   42.0,   13.0,   19.0,  153.0,   51.0,   34.0,   94.0,   12.0,   32.0,   33.0,   17.0,   11.0,    0.0,  245.0,   78.0,    0.0,    0.0,   48.0],
                  [ 409.0,  154.0,  495.0,   95.0,  161.0,   56.0,   79.0,  234.0,   35.0,   24.0,   17.0,   96.0,   62.0,   46.0,  245.0,    0.0,  550.0,   75.0,   34.0,   30.0],
                  [ 371.0,   26.0,  229.0,   66.0,   16.0,   53.0,   34.0,   30.0,   22.0,  192.0,   33.0,  136.0,  104.0,   13.0,   78.0,  550.0,    0.0,    0.0,   42.0,  157.0],
                  [   0.0,  201.0,   23.0,    0.0,    0.0,    0.0,    0.0,    0.0,   27.0,    0.0,   46.0,    0.0,    0.0,   76.0,    0.0,   75.0,    0.0,    0.0,   61.0,    0.0],
                  [  24.0,    8.0,   95.0,    0.0,   96.0,    0.0,   22.0,    0.0,  127.0,   37.0,   28.0,   13.0,    0.0,  698.0,    0.0,   34.0,   42.0,   61.0,    0.0,   28.0],
                  [ 208.0,   24.0,   15.0,   18.0,   49.0,   35.0,   37.0,   54.0,   44.0,  889.0,  175.0,   10.0,  258.0,   12.0,   48.0,   30.0,  157.0,    0.0,   28.0,    0.0]])
    
    freqs = np.array([0.087127, 0.040904, 0.040432, 0.046872, 0.033474, 0.038255,  0.04953, 0.088612, 0.033618, 0.036886, 0.085357, 0.080482, 0.014753, 0.039772,  0.05068, 0.069577, 0.058542, 0.010494, 0.029916, 0.064718])
    
    return S, freqs


empirical_models = {
                    'WAG': _WAG_aa,
                    'JTT': _JTT_aa,
                    'BLOSUM62': _BLOSUM62_aa,
                    'LG': _LG_aa,
                    'DAYHOFF': _DAYHOFF_aa,
                   }