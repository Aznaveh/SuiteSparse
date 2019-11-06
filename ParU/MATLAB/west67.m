% A  is  67 x 67 
LU = zeros(67,67);
npivots =[]; 
S = zeros(67,67); % n1 = 1
% nf=33
InMatrix=[
1,1, 0.2500000000000000;
1,18, -0.2362845000000000;
1,19, -0.9722222000000000;
1,61, 0.7222222000000000;
1,65, 0.1278394000000000;
2,1, -0.8242248000000000;
2,5, -0.2541193000000000;
2,6, 0.7453416000000000;
3,1, 0.5000000000000000;
3,6, 0.4444444000000000;
3,20, -0.2667757000000000;
3,23, -0.9444444000000000;
4,1, 1.0000000000000000;
4,14, 1.0000000000000000;
4,15, 1.0000000000000000;
4,16, 1.0000000000000000;
4,17, 1.0000000000000000;
5,2, 0.4444444000000000;
5,17, 0.5000000000000000;
5,20, -0.2630706000000000;
5,24, -0.9444444000000000;
6,2, -0.3726708000000000;
6,3, -1.8633540000000000;
6,4, -1.1180120000000000;
6,5, 1.0000000000000000;
6,6, -0.7453416000000000;
6,13, -1.4906830000000000;
7,2, 1.0000000000000000;
7,3, 1.0000000000000000;
7,4, 1.0000000000000000;
7,6, 1.0000000000000000;
7,13, 1.0000000000000000;
8,3, 1.8633540000000000;
8,5, -0.1443354000000000;
8,14, -0.8242248000000000;
9,3, 0.4444444000000000;
9,14, 0.5000000000000000;
9,20, -0.1064573000000000;
9,21, -0.9444444000000000;
10,4, 1.1180120000000000;
10,5, -0.2421498000000000;
10,16, -0.8242248000000000;
11,4, 0.4444444000000000;
11,16, 0.5000000000000000;
11,20, -0.2122056000000000;
11,22, -0.9444444000000000;
12,5, -0.1918557000000000;
12,13, 1.4906830000000000;
12,15, -0.8242248000000000;
13,7, -0.9722222000000000;
13,17, 0.2500000000000000;
13,18, -0.1803900000000000;
13,62, 0.7222222000000000;
13,65, 0.1026879000000000;
14,7, -0.9722222000000000;
14,18, -0.1803900000000000;
14,20, 0.1315353000000000;
14,24, 0.4722222000000000;
14,57, 0.5000000000000000;
15,7, 1.0000000000000000;
15,8, 1.0000000000000000;
15,9, 1.0000000000000000;
15,12, 1.0000000000000000;
15,19, 1.0000000000000000;
16,8, -0.9722222000000000;
16,14, 0.2500000000000000;
16,18, -0.1581626000000000;
16,56, 0.7222222000000000;
16,65, 0.0905272100000000;
17,8, -0.9722222000000000;
17,18, -0.1581626000000000;
17,20, 0.0532286400000000;
17,21, 0.4722222000000000;
17,59, 0.5000000000000000;
18,9, -0.9722222000000000;
18,16, 0.2500000000000000;
18,18, -0.2303917000000000;
18,63, 0.7222222000000000;
18,65, 0.1257342000000000;
19,9, -0.9722222000000000;
19,18, -0.2303917000000000;
19,20, 0.1061028000000000;
19,22, 0.4722222000000000;
19,60, 0.5000000000000000;
20,10, -0.2069954000000000;
20,21, 1.8633540000000000;
20,59, -0.9583187000000000;
21,10, -0.2356469000000000;
21,11, 1.4906830000000000;
21,66, -0.9583187000000000;
22,10, -0.2475675000000000;
22,22, 1.1180120000000000;
22,60, -0.9583187000000000;
23,10, -0.2074873000000000;
23,23, 0.7453416000000000;
23,58, -0.9583187000000000;
24,10, 1.0000000000000000;
24,11, -1.4906830000000000;
24,21, -1.8633540000000000;
24,22, -1.1180120000000000;
24,23, -0.7453416000000000;
24,24, -0.3726708000000000;
25,11, 0.4722222000000000;
25,12, -0.9722222000000000;
25,18, -0.1947711000000000;
25,20, 0.0757454200000000;
25,66, 0.5000000000000000;
26,11, -0.9444444000000000;
26,13, 0.4444444000000000;
26,15, 0.5000000000000000;
26,20, -0.1514908000000000;
27,11, 1.0000000000000000;
27,21, 1.0000000000000000;
27,22, 1.0000000000000000;
27,23, 1.0000000000000000;
27,24, 1.0000000000000000;
28,12, -0.9722222000000000;
28,15, 0.2500000000000000;
28,18, -0.1947711000000000;
28,64, 0.7222222000000000;
28,65, 0.1087668000000000;
29,18, -0.2362845000000000;
29,19, -0.9722222000000000;
29,20, 0.1333878000000000;
29,23, 0.4722222000000000;
29,58, 0.5000000000000000;
30,25, 0.4000000000000000;
30,41, 0.0117829100000000;
30,42, -0.0632597800000000;
30,45, -0.8000000000000000;
30,52, 0.4000000000000000;
31,25, 1.0000000000000000;
31,27, 1.0000000000000000;
31,28, 1.0000000000000000;
31,46, 1.0000000000000000;
31,47, 1.0000000000000000;
32,26, -0.3361556000000000;
32,46, -0.8341818000000000;
32,48, 1.2658229999999999;
33,26, -0.2939196000000000;
33,27, -0.8341818000000000;
33,49, 1.0126580000000001;
34,26, -0.2214815000000000;
34,47, -0.8341818000000000;
34,50, 0.7594937000000000;
35,26, -0.1189860000000000;
35,28, -0.8341818000000000;
35,51, 0.5063291000000000;
36,26, 1.0000000000000000;
36,48, -1.2658229999999999;
36,49, -1.0126580000000001;
36,50, -0.7594937000000000;
36,51, -0.5063291000000000;
36,52, -0.2531646000000000;
37,27, 0.4000000000000000;
37,41, 0.1175679000000000;
37,42, -0.2680186000000000;
37,43, -0.8000000000000000;
37,49, 0.4000000000000000;
38,28, 0.4000000000000000;
38,41, 0.0475943900000000;
38,42, -0.1575082000000000;
38,44, -0.8000000000000000;
38,51, 0.4000000000000000;
39,29, -0.2788416000000000;
39,30, -0.9159533000000000;
39,56, 1.5673980000000001;
40,29, -0.2680186000000000;
40,43, -0.9159533000000000;
40,64, 1.2539180000000001;
41,29, -0.2323717000000000;
41,31, -0.9159533000000000;
41,63, 0.9404389000000000;
42,29, -0.1575082000000000;
42,44, -0.9159533000000000;
42,61, 0.6269592000000000;
43,29, 1.0000000000000000;
43,56, -1.5673980000000001;
43,61, -0.6269592000000000;
43,62, -0.3134796000000000;
43,63, -0.9404389000000000;
43,64, -1.2539180000000001;
44,30, -0.8000000000000000;
44,41, 0.1344622000000000;
44,42, -0.2788416000000000;
44,46, 0.4000000000000000;
44,48, 0.4000000000000000;
45,30, 0.4000000000000000;
45,38, -0.2070986000000000;
45,39, -1.0500000000000000;
45,42, 0.1394208000000000;
45,56, 0.6500000000000000;
45,65, 0.0814744900000000;
46,30, 1.0000000000000000;
46,31, 1.0000000000000000;
46,43, 1.0000000000000000;
46,44, 1.0000000000000000;
46,45, 1.0000000000000000;
47,31, -0.8000000000000000;
47,41, 0.0885926200000000;
47,42, -0.2323717000000000;
47,47, 0.4000000000000000;
47,50, 0.4000000000000000;
48,31, 0.4000000000000000;
48,38, -0.2286264000000000;
48,40, -1.0500000000000000;
48,42, 0.1161859000000000;
48,63, 0.6500000000000000;
48,65, 0.1131608000000000;
49,32, -1.0500000000000000;
49,38, -0.2024528000000000;
49,53, 0.1192061000000000;
49,55, 0.6000000000000000;
49,58, 0.4500000000000000;
50,32, -1.0500000000000000;
50,38, -0.2024528000000000;
50,42, 0.0787541100000000;
50,44, 0.4000000000000000;
50,61, 0.6500000000000000;
50,65, 0.1150555000000000;
51,32, 1.0000000000000000;
51,33, 1.0000000000000000;
51,37, 1.0000000000000000;
51,39, 1.0000000000000000;
51,40, 1.0000000000000000;
52,33, -1.0500000000000000;
52,38, -0.2232997000000000;
52,53, 0.1284235000000000;
52,54, 0.6000000000000000;
52,66, 0.4500000000000000;
53,33, -1.0500000000000000;
53,38, -0.2232997000000000;
53,42, 0.1340093000000000;
53,43, 0.4000000000000000;
53,64, 0.6500000000000000;
53,65, 0.0978901500000000;
54,34, -1.0000000000000000;
54,52, 0.3333333000000000;
54,53, -0.1656874000000000;
55,34, 0.6000000000000000;
55,37, -1.0500000000000000;
55,38, -0.1385226000000000;
55,53, 0.0994124600000000;
55,57, 0.4500000000000000;
56,34, 1.0000000000000000;
56,35, 1.0000000000000000;
56,36, 1.0000000000000000;
56,54, 1.0000000000000000;
56,55, 1.0000000000000000;
57,35, -1.0000000000000000;
57,48, 0.3333333000000000;
57,53, -0.2071759000000000;
58,35, 0.6000000000000000;
58,38, -0.2070986000000000;
58,39, -1.0500000000000000;
58,53, 0.1243055000000000;
58,59, 0.4500000000000000;
59,36, -1.0000000000000000;
59,50, 0.3333333000000000;
59,53, -0.2144206000000000;
60,36, 0.6000000000000000;
60,38, -0.2286264000000000;
60,40, -1.0500000000000000;
60,53, 0.1286524000000000;
60,60, 0.4500000000000000;
61,37, -1.0500000000000000;
61,38, -0.1385226000000000;
61,42, 0.0316298900000000;
61,45, 0.4000000000000000;
61,62, 0.6500000000000000;
61,65, 0.0924190900000000;
62,48, 1.0000000000000000;
62,49, 1.0000000000000000;
62,50, 1.0000000000000000;
62,51, 1.0000000000000000;
62,52, 1.0000000000000000;
63,49, 0.3333333000000000;
63,53, -0.2140392000000000;
63,54, -1.0000000000000000;
64,51, 0.3333333000000000;
64,53, -0.1986768000000000;
64,55, -1.0000000000000000;
65,56, 1.0000000000000000;
65,61, 1.0000000000000000;
65,62, 1.0000000000000000;
65,63, 1.0000000000000000;
65,64, 1.0000000000000000;
66,57, 1.0000000000000000;
66,58, 1.0000000000000000;
66,59, 1.0000000000000000;
66,60, 1.0000000000000000;
66,66, 1.0000000000000000;
];
I = InMatrix(:,1);
J = InMatrix(:,2);
X = InMatrix(:,3);
S = sparse(I,J,X);
