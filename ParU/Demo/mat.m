%%%%%%%% START %%%%%%%%%%%%%%%%%%% 
clear all
% A  is  37 x 37 
LU = zeros(37,37);
npivots =[]; 
S = zeros(37,37);
% nf=11
% anz = 233  rjsize=144
S(1,1)= 0.1800664010750090;
S(1,2)= 0.1800664010750090;
S(1,23)= 0.6366555998208320;
S(1,24)= 0.1466223992833270;
S(1,25)= 0.0733111996416637;
S(1,27)= 0.1466223992833270;
S(1,29)= 0.0900332005375044;
S(1,30)= 0.0733111996416637;
S(1,37)= 0.1466223992833270;
S(2,1)= 0.5099667994624950;
S(2,2)= 0.0733111996416637;
S(2,23)= 0.0366555998208319;
S(2,26)= 0.0366555998208319;
S(2,29)= 0.0733111996416637;
S(3,1)= 0.0900332005375044;
S(3,25)= 0.0733111996416637;
S(3,26)= 0.6366555998208320;
S(3,29)= 0.0900332005375044;
S(3,30)= 0.0733111996416637;
S(3,31)= 0.1466223992833270;
S(4,1)= 0.1466223992833270;
S(4,2)= 0.0733111996416637;
S(4,19)= 0.0300110668458348;
S(4,20)= 0.0300110668458348;
S(4,23)= 0.0366555998208319;
S(4,26)= 0.0733111996416637;
S(4,27)= 0.0733111996416637;
S(4,29)= 0.4132890659499940;
S(4,36)= 0.1466223992833270;
S(5,1)= 0.0733111996416637;
S(5,2)= 0.5466223992833270;
S(5,23)= 0.0366555998208319;
S(5,29)= 0.0366555998208319;
S(5,31)= 0.0733111996416637;
S(5,36)= 0.0733111996416637;
S(6,2)= 0.0900332005375044;
S(6,25)= 0.0366555998208319;
S(6,26)= 0.0733111996416637;
S(6,30)= 0.0366555998208319;
S(6,31)= 0.5633444001791680;
S(6,36)= 0.0900332005375044;
S(7,2)= 0.0366555998208319;
S(7,27)= 0.0366555998208319;
S(7,29)= 0.0366555998208319;
S(7,31)= 0.0366555998208319;
S(7,33)= 0.0300110668458348;
S(7,34)= 0.0600221336916696;
S(7,36)= 0.1699889331541650;
S(8,3)= 0.1200442673833390;
S(8,4)= 0.7699889331541649;
S(8,5)= 0.2199335989249910;
S(8,17)= 0.1099667994624960;
S(8,21)= 0.0600221336916696;
S(9,3)= 0.8199335989249910;
S(9,4)= 0.1099667994624960;
S(9,21)= 0.1099667994624960;
S(10,3)= 0.0600221336916696;
S(10,4)= 0.0300110668458348;
S(10,17)= 0.1099667994624960;
S(10,21)= 0.6132890659499940;
S(10,22)= 0.0733111996416637;
S(10,30)= 0.0733111996416637;
S(10,32)= 0.2199335989249910;
S(11,4)= 0.0600221336916696;
S(11,5)= 0.7200442673833390;
S(11,17)= 0.0600221336916696;
S(12,4)= 0.0300110668458348;
S(12,5)= 0.0600221336916696;
S(12,17)= 0.5633444001791680;
S(12,18)= 0.0733111996416637;
S(12,21)= 0.0600221336916696;
S(12,28)= 0.0733111996416637;
S(12,32)= 0.1200442673833390;
S(13,6)= 0.1099667994624960;
S(13,23)= 0.0600221336916696;
S(13,24)= 0.6300110668458350;
S(13,27)= 0.0600221336916696;
S(13,28)= 0.1099667994624960;
S(13,37)= 0.1200442673833390;
S(14,6)= 0.2199335989249910;
S(14,7)= 0.2199335989249910;
S(14,23)= 0.0600221336916696;
S(14,25)= 0.6300110668458350;
S(14,26)= 0.1200442673833390;
S(14,27)= 0.0600221336916696;
S(14,28)= 0.1099667994624960;
S(14,30)= 0.1200442673833390;
S(14,31)= 0.1200442673833390;
S(15,6)= 0.4900332005375040;
S(15,7)= 0.0600221336916696;
S(15,24)= 0.0300110668458348;
S(15,25)= 0.0300110668458348;
S(15,28)= 0.0600221336916696;
S(16,6)= 0.1200442673833390;
S(16,7)= 0.0600221336916696;
S(16,17)= 0.0366555998208319;
S(16,18)= 0.0366555998208319;
S(16,24)= 0.0600221336916696;
S(16,25)= 0.0300110668458348;
S(16,28)= 0.3867109340500060;
S(16,30)= 0.0600221336916696;
S(16,35)= 0.1200442673833390;
S(17,6)= 0.0600221336916696;
S(17,7)= 0.5200442673833390;
S(17,25)= 0.0300110668458348;
S(17,28)= 0.0300110668458348;
S(17,35)= 0.0600221336916696;
S(17,37)= 0.0600221336916696;
S(18,7)= 0.1099667994624960;
S(18,23)= 0.0300110668458348;
S(18,24)= 0.0600221336916696;
S(18,27)= 0.0300110668458348;
S(18,35)= 0.1099667994624960;
S(18,37)= 0.5699889331541650;
S(19,7)= 0.0300110668458348;
S(19,28)= 0.0300110668458348;
S(19,30)= 0.0300110668458348;
S(19,32)= 0.0733111996416637;
S(19,33)= 0.0366555998208319;
S(19,35)= 0.1633444001791680;
S(19,37)= 0.0300110668458348;
S(20,8)= 0.1800664010750090;
S(20,9)= 0.7633444001791680;
S(20,10)= 0.1466223992833270;
S(20,15)= 0.0733111996416637;
S(20,20)= 0.0900332005375044;
S(21,8)= 0.7466223992833270;
S(21,9)= 0.0733111996416637;
S(21,20)= 0.0733111996416637;
S(22,8)= 0.0733111996416637;
S(22,9)= 0.0366555998208319;
S(22,15)= 0.0733111996416637;
S(22,19)= 0.0600221336916696;
S(22,20)= 0.5699889331541650;
S(22,29)= 0.0600221336916696;
S(22,34)= 0.1466223992833270;
S(23,9)= 0.0900332005375044;
S(23,10)= 0.7800664010750090;
S(23,15)= 0.0900332005375044;
S(24,9)= 0.0366555998208319;
S(24,10)= 0.0733111996416637;
S(24,15)= 0.5867109340500060;
S(24,16)= 0.0600221336916696;
S(24,20)= 0.0900332005375044;
S(24,27)= 0.0600221336916696;
S(24,34)= 0.1800664010750090;
S(25,11)= 0.0900332005375044;
S(25,13)= 0.8000000000000000;
S(25,14)= 0.1099667994624960;
S(25,16)= 0.1099667994624960;
S(25,22)= 0.0900332005375044;
S(26,11)= 0.0600221336916696;
S(26,12)= 0.7333333333333329;
S(26,14)= 0.0733111996416637;
S(26,18)= 0.0733111996416637;
S(26,19)= 0.0600221336916696;
S(27,11)= 0.7832779991041590;
S(27,12)= 0.1099667994624960;
S(27,13)= 0.0733111996416637;
S(27,19)= 0.1099667994624960;
S(27,22)= 0.0733111996416637;
S(28,11)= 0.0300110668458348;
S(28,12)= 0.0300110668458348;
S(28,16)= 0.0733111996416637;
S(28,19)= 0.6132890659499940;
S(28,20)= 0.1099667994624960;
S(28,29)= 0.1099667994624960;
S(28,33)= 0.0733111996416637;
S(29,11)= 0.0366555998208319;
S(29,13)= 0.0366555998208319;
S(29,18)= 0.1099667994624960;
S(29,21)= 0.0900332005375044;
S(29,22)= 0.6366555998208320;
S(29,30)= 0.0900332005375044;
S(29,33)= 0.1099667994624960;
S(30,12)= 0.0900332005375044;
S(30,13)= 0.0600221336916696;
S(30,14)= 0.7500553342291740;
S(30,16)= 0.0600221336916696;
S(30,18)= 0.0900332005375044;
S(31,12)= 0.0366555998208319;
S(31,14)= 0.0366555998208319;
S(31,17)= 0.0900332005375044;
S(31,18)= 0.5867109340500060;
S(31,22)= 0.0600221336916696;
S(31,28)= 0.0900332005375044;
S(31,33)= 0.0600221336916696;
S(32,13)= 0.0300110668458348;
S(32,14)= 0.0300110668458348;
S(32,15)= 0.1099667994624960;
S(32,16)= 0.6300110668458350;
S(32,19)= 0.0900332005375044;
S(32,27)= 0.1099667994624960;
S(32,33)= 0.0900332005375044;
S(33,15)= 0.0300110668458348;
S(33,16)= 0.0300110668458348;
S(33,23)= 0.0733111996416637;
S(33,24)= 0.0733111996416637;
S(33,25)= 0.0366555998208319;
S(33,27)= 0.3933554670250030;
S(33,29)= 0.0900332005375044;
S(33,30)= 0.0366555998208319;
S(33,36)= 0.1800664010750090;
S(33,37)= 0.0733111996416637;
S(34,15)= 0.0366555998208319;
S(34,20)= 0.0366555998208319;
S(34,33)= 0.0600221336916696;
S(34,34)= 0.3933554670250030;
S(34,36)= 0.1200442673833390;
S(35,16)= 0.0366555998208319;
S(35,18)= 0.0300110668458348;
S(35,19)= 0.0366555998208319;
S(35,22)= 0.0300110668458348;
S(35,32)= 0.1800664010750090;
S(35,33)= 0.4666666666666670;
S(35,34)= 0.2199335989249910;
S(35,35)= 0.1800664010750090;
S(35,36)= 0.2199335989249910;
S(36,17)= 0.0300110668458348;
S(36,21)= 0.0300110668458348;
S(36,32)= 0.4066445329749970;
S(36,33)= 0.0733111996416637;
S(36,35)= 0.1466223992833270;
S(37,21)= 0.0366555998208319;
S(37,22)= 0.0366555998208319;
S(37,23)= 0.0300110668458348;
S(37,25)= 0.0600221336916696;
S(37,26)= 0.0600221336916696;
S(37,27)= 0.0300110668458348;
S(37,28)= 0.1099667994624960;
S(37,30)= 0.4066445329749970;
S(37,31)= 0.0600221336916696;
S(37,35)= 0.2199335989249910;
%~~~~~~~  Assemble Front 0 start ~~~~
% fp=1 pivotal columns:clo1=0...col2=0
rowMark=0;
%L part:
cols{1} = [1 ];
rows{1} = [2 1 3 4 5 ];
Luf{1}= [  0.5099667994624950 ;
     0.3530943607795625 ;
     0.1765471803897810 ;
     0.2875136174313053 ;
     0.1437568087156531 ;
   ];
Us{1} =[];
Ucols{1}=[];
Urows{1}=[];
% rowCount=5;
% U part After TRSM: 1 x 4
Ucols{1} = [2 23 26 29 ];
Urows{1} = [2 ];
Us{1} = [ 0.0733111996416637  0.0366555998208319  0.0366555998208319  0.0733111996416637 ;
    ];

%After DGEMM:
% Element 5 is 4 x 4:
	% 1	% 22	% 25	% 28	
% 0	-0.0259	-0.0129	-0.0129	-0.0259	
% 2	-0.0129	-0.0065	-0.0065	-0.0129	
% 3	-0.0211	-0.0105	-0.0105	-0.0211	
% 4	-0.0105	-0.0053	-0.0053	-0.0105	

%||||  Start FourPass 0 ||||

%||||  Finish FourPass 0 ||||
fp =1;
%~~~~~~~  Assemble Front 1 start ~~~~
% fp=1 pivotal columns:clo1=1...col2=1
rowMark=5;
%L part:
cols{2} = [2 ];
rows{2} = [5 4 1 6 7 3 ];
Luf{2}= [  0.5360834151797254 ;
     0.0974348953081282 ;
     0.2876056701882166 ;
     0.1679462523706692 ;
     0.0683766719560665 ;
     -0.0241434172765611 ;
   ];
Us{2} =[];
Ucols{2}=[];
Urows{2}=[];
% rowCount=6;
% U part After TRSM: 1 x 5
Ucols{2} = [23 29 31 36 26 ];
Urows{2} = [5 ];
Us{2} = [ 0.0313861077690310  0.0261166157172302  0.0733111996416637  0.0733111996416637  -0.0052694920518009 ;
    ];
%Row Tuple reallocating space
%Row Tuple reallocating space
%Row Tuple reallocating space

%After DGEMM:
% Element 8 is 5 x 5:
	% 22	% 28	% 30	% 35	% 25	
% 3	-0.0031	-0.0025	-0.0071	-0.0071	0.0005	
% 0	-0.0090	-0.0075	-0.0211	-0.0211	0.0015	
% 5	-0.0053	-0.0044	-0.0123	-0.0123	0.0009	
% 6	-0.0021	-0.0018	-0.0050	-0.0050	0.0004	
% 2	0.0008	0.0006	0.0018	0.0018	-0.0001	

%||||  Start FourPass 1 ||||

%||||  Finish FourPass 1 ||||
fp =1;
%~~~~~~~  Assemble Front 2 start ~~~~
% fp=1 pivotal columns:clo1=2...col2=2
rowMark=11;
%L part:
cols{3} = [3 ];
rows{3} = [9 8 10 ];
Luf{3}= [  0.8199335989249910 ;
     0.1464073036410853 ;
     0.0732036518205428 ;
   ];
Us{3} =[];
Ucols{3}=[];
Urows{3}=[];
% rowCount=3;
% U part After TRSM: 1 x 2
Ucols{3} = [4 21 ];
Urows{3} = [9 ];
Us{3} = [ 0.1099667994624960  0.1099667994624960 ;
    ];

%After DGEMM:
% Element 12 is 2 x 2:
	% 3	% 20	
% 7	-0.0161	-0.0161	
% 9	-0.0080	-0.0080	

%||||  Start FourPass 2 ||||

%||||  Finish FourPass 2 ||||
fp =1;
%~~~~~~~  Assemble Front 3 start ~~~~
% fp=2 pivotal columns:clo1=3...col2=4
rowMark=14;
%L part:
cols{4} = [4 5 ];
rows{4} = [8 11 10 12 ];
Luf{4}= [  0.7538889905548209  0.2199335989249910 ;
     0.0796166736000437  0.7025338858240451 ;
     0.0291304102080077  -0.0091194974142686 ;
     0.0398083368000218  0.0729743346854913 ;
   ];
Us{4} =[];
Ucols{4}=[];
Urows{4}=[];
% rowCount=4;
% U part After TRSM: 2 x 2
Ucols{4} = [17 21 ];
Urows{4} = [8 11 ];
Us{4} = [ 0.1099667994624960  0.0439221910923256 ;
     0.0512669429120226  -0.0034969387519964 ;
    ];
%Row Tuple reallocating space

%After DGEMM:
% Element 15 is 2 x 2:
	% 16	% 20	
% 9	-0.0027	-0.0013	
% 11	-0.0081	-0.0015	

%||||  Start FourPass 3 ||||

%||||  Finish FourPass 3 ||||
fp =2;
%~~~~~~~  Assemble Front 4 start ~~~~
% fp=1 pivotal columns:clo1=5...col2=5
rowMark=18;
%L part:
cols{5} = [6 ];
rows{5} = [15 14 13 16 17 ];
Luf{5}= [  0.4900332005375040 ;
     0.4488136695304560 ;
     0.2244068347652290 ;
     0.2449717024309082 ;
     0.1224858512154543 ;
   ];
Us{5} =[];
Ucols{5}=[];
Urows{5}=[];
% rowCount=5;
% U part After TRSM: 1 x 4
Ucols{5} = [7 24 25 28 ];
Urows{5} = [15 ];
Us{5} = [ 0.0600221336916696  0.0300110668458348  0.0300110668458348  0.0600221336916696 ;
    ];

%After DGEMM:
% Element 21 is 4 x 4:
	% 6	% 23	% 24	% 27	
% 13	-0.0269	-0.0135	-0.0135	-0.0269	
% 12	-0.0135	-0.0067	-0.0067	-0.0135	
% 15	-0.0147	-0.0074	-0.0074	-0.0147	
% 16	-0.0074	-0.0037	-0.0037	-0.0074	

%||||  Start FourPass 4 ||||

%||||  Finish FourPass 4 ||||
fp =1;
%~~~~~~~  Assemble Front 5 start ~~~~
% fp=1 pivotal columns:clo1=6...col2=6
rowMark=23;
%L part:
cols{6} = [7 ];
rows{6} = [17 16 14 18 19 13 ];
Luf{6}= [  0.5126924052463471 ;
     0.0883929797944059 ;
     0.3764339843439884 ;
     0.2144888403596642 ;
     0.0585362032648301 ;
     -0.0262718481876692 ;
   ];
Us{6} =[];
Ucols{6}=[];
Urows{6}=[];
% rowCount=6;
% U part After TRSM: 1 x 5
Ucols{6} = [25 28 35 37 24 ];
Urows{6} = [17 ];
Us{6} = [ 0.0263351357773388  0.0226592047088429  0.0600221336916696  0.0600221336916696  -0.0036759310684960 ;
    ];
%Row Tuple reallocating space
%Row Tuple reallocating space
%Row Tuple reallocating space

%After DGEMM:
% Element 24 is 5 x 5:
	% 24	% 27	% 34	% 36	% 23	
% 15	-0.0023	-0.0020	-0.0053	-0.0053	0.0003	
% 13	-0.0099	-0.0085	-0.0226	-0.0226	0.0014	
% 17	-0.0056	-0.0049	-0.0129	-0.0129	0.0008	
% 18	-0.0015	-0.0013	-0.0035	-0.0035	0.0002	
% 12	0.0007	0.0006	0.0016	0.0016	-0.0001	

%||||  Start FourPass 5 ||||

%||||  Finish FourPass 5 ||||
fp =1;
%~~~~~~~  Assemble Front 6 start ~~~~
% fp=1 pivotal columns:clo1=7...col2=7
rowMark=29;
%L part:
cols{7} = [8 ];
rows{7} = [21 20 22 ];
Luf{7}= [  0.7466223992833270 ;
     0.2411746570258974 ;
     0.0981904637632546 ;
   ];
Us{7} =[];
Ucols{7}=[];
Urows{7}=[];
% rowCount=3;
% U part After TRSM: 1 x 2
Ucols{7} = [9 20 ];
Urows{7} = [21 ];
Us{7} = [ 0.0733111996416637  0.0733111996416637 ;
    ];

%After DGEMM:
% Element 28 is 2 x 2:
	% 8	% 19	
% 19	-0.0177	-0.0177	
% 21	-0.0072	-0.0072	

%||||  Start FourPass 6 ||||

%||||  Finish FourPass 6 ||||
fp =1;
%~~~~~~~  Assemble Front 7 start ~~~~
% fp=2 pivotal columns:clo1=8...col2=9
rowMark=32;
%L part:
cols{8} = [9 10 ];
rows{8} = [20 23 22 24 ];
Luf{8}= [  0.7456635967494327  0.1466223992833270 ;
     0.1207423842735325  0.7623628629976342 ;
     0.0395045959832138  -0.0075977712542318 ;
     0.0491583603928426  0.0867086870366556 ;
   ];
Us{8} =[];
Ucols{8}=[];
Urows{8}=[];
% rowCount=4;
% U part After TRSM: 2 x 2
Ucols{8} = [15 20 ];
Urows{8} = [20 23 ];
Us{8} = [ 0.0733111996416637  0.0723523971077691 ;
     0.0811814314988170  -0.0087360009346975 ;
    ];
%Row Tuple reallocating space

%After DGEMM:
% Element 31 is 2 x 2:
	% 14	% 19	
% 21	-0.0023	-0.0029	
% 23	-0.0106	-0.0028	

%||||  Start FourPass 7 ||||

%||||  Finish FourPass 7 ||||
fp =2;
%~~~~~~~  Assemble Front 8 start ~~~~
% fp=1 pivotal columns:clo1=10...col2=10
rowMark=36;
%L part:
cols{9} = [11 ];
rows{9} = [27 26 25 28 29 ];
Luf{9}= [  0.7832779991041590 ;
     0.0766294135164238 ;
     0.1149441202746356 ;
     0.0383147067582119 ;
     0.0467976885125782 ;
   ];
Us{9} =[];
Ucols{9}=[];
Urows{9}=[];
% rowCount=5;
% U part After TRSM: 1 x 4
Ucols{9} = [12 13 19 22 ];
Urows{9} = [27 ];
Us{9} = [ 0.1099667994624960  0.0733111996416637  0.1099667994624960  0.0733111996416637 ;
    ];

%After DGEMM:
% Element 37 is 4 x 4:
	% 11	% 12	% 18	% 21	
% 25	-0.0084	-0.0056	-0.0084	-0.0056	
% 24	-0.0126	-0.0084	-0.0126	-0.0084	
% 27	-0.0042	-0.0028	-0.0042	-0.0028	
% 28	-0.0051	-0.0034	-0.0051	-0.0034	

%||||  Start FourPass 8 ||||

%||||  Finish FourPass 8 ||||
fp =1;
%~~~~~~~  Assemble Front 9 start ~~~~
% fp=3 pivotal columns:clo1=11...col2=13
rowMark=41;
%L part:
cols{10} = [12 13 14 ];
rows{10} = [26 25 30 31 28 29 32 ];
Luf{10}= [  0.7249066419842437  -0.0056177942327261  0.0733111996416637 ;
     -0.0174367791541198  0.7914753524135415  0.1112451106601713 ;
     0.1241997180368797  0.0767173127580069  0.7324156779574186 ;
     0.0505659593909868  0.0003589109302933  0.0449316181795369 ;
     0.0355876462942533  -0.0032963414897243  -0.0030614734718816 ;
     -0.0070991100507573  0.0419279307875128  -0.0056577612286456 ;
     0.0000000000000000  0.0379178792546330  0.0352161879495617 ;
   ];
Us{10} =[];
Ucols{10}=[];
Urows{10}=[];
% rowCount=7;
% U part After TRSM: 3 x 4
Ucols{10} = [18 19 22 16 ];
Urows{10} = [26 25 30 ];
Us{10} = [ 0.0733111996416637  0.0515954423425803  -0.0056177942327261  0.0000000000000000 ;
     0.0012783111976753  -0.0117403786901472  0.0815085529510459  0.1099667994624960 ;
     0.0808299016131103  -0.0055074490870671  -0.0055553886895042  0.0515857763443083 ;
    ];
%Row Tuple reallocating space
%Row Tuple reallocating space

%After DGEMM:
% Element 41 is 4 x 4:
	% 17	% 18	% 21	% 15	
% 30	-0.0073	-0.0024	0.0005	-0.0024	
% 27	-0.0024	-0.0019	0.0005	0.0005	
% 28	0.0009	0.0008	-0.0035	-0.0043	
% 31	-0.0029	0.0006	-0.0029	-0.0060	

%||||  Start FourPass 9 ||||

%||||  Finish FourPass 9 ||||
fp =3;
%~~~~~~~  Assemble Front 10 start ~~~~
% fp=23 pivotal columns:clo1=14...col2=36
rowMark=48;
%L part:
cols{11} = [15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 ];
rows{11} = [24 32 12 31 28 22 10 29 1 13 14 3 33 16 4 37 6 36 35 34 19 7 18 ];
Luf{11}= [  0.5760679403401708  0.0600221336916696  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0872339624961932  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0600221336916696  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.1800664010750090  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.1908920663030824  0.6125669555022097  0.0000000000000000  -0.0028949918567952  0.0906723221612539  -0.0166522713527039  0.0000000000000000  -0.0028949918567952  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0985090503381733  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0900332005375044  -0.0343732473729680  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  0.0000000000000000  0.5552256337389806  0.0733111996416637  0.0000000000000000  0.0000000000000000  0.0585288510945342  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0733111996416637  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.1200442673833390  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  -0.0038482335217728  0.1621560588462136  0.5674726100215014  -0.0020083724228560  -0.0000640818288331  -0.0094908078222866  0.0605154205348692  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0003790858297094  0.0781453453343243  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0194659052859633  0.0603686024720505  -0.0001322762827928  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  0.1205282387302433  0.0000000000000000  -0.0035391530575897  0.5962483154872286  0.1119736416040951  -0.0000335894215232  -0.0017941990872387  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0118717806934756  0.0002765683378764  0.1099667994624960  0.0000000000000000  0.0000000000000000  -0.0000688928182116  0.0626733102776529  0.0041424788192921  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.1233046717084264  -0.0120819600594076  0.0000000000000000  -0.0000616367651379  0.1025034498034492  0.5374306104233001  0.0000028580488903  0.0001526643948528  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0049939052363991  -0.0000235325824417  0.0487501573829197  0.0000000000000000  0.0000000000000000  0.0000058619361011  -0.0053327320557744  0.1235794480883855  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  0.0000000000000000  0.1931304027087233  -0.0249503170017627  -0.0000840413755583  0.0000145350055835  0.5923872307724832  0.0748209255615665  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000085331772828  -0.0122088467863207  0.0000085331772828  0.0733111996416637  0.0000000000000000  0.1962637148397115  0.0015115604311278  -0.0000047484235373  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  -0.0070503645563066  0.0000000000000000  0.1953760417169161  -0.0055130714592905  0.0009534889679117  0.1551135546535952  0.6062765568537196  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0005584484748768  -0.0133724234477066  0.0005584484748768  0.0786616397651666  0.0000000000000000  -0.0266403763354205  0.0989231313812743  -0.0003107576222281  0.0000000000000000  0.0000000000000000  0.0000000000000000 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.6146858916735921  0.1466223992833270  0.0733111996416637  -0.0114273497946184  0.1466223992833270  0.0000000000000000  0.0566361425956462  0.0733111996416637  -0.0210847167052428  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0210847167052428  0.1466223992833270 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0976468380106279  0.6088625911525050  -0.0132014226654319  0.0011158445742859  0.0457049200201210  0.0970927216110571  -0.0055303402415839  -0.0071586068357743  0.0020588559166168  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0015768923842475  0.0020588559166168  0.1073039460960379 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0976468380106279  -0.0433642101100123  0.5988971736173393  0.1212084996761943  0.0476868777749348  0.0787086998567692  -0.0057701590777998  0.1125752332166432  0.1221924039605102  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0225259902417113  0.0021481365771712  -0.0322584337417893 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0092952758062377  0.0022384289337452  0.1235975075009720  0.6149671469267598  -0.0046333908081622  -0.0099455342782490  0.0789728644204273  0.0600946532627679  0.1330891087425242  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0027806264863284  0.0013038816975581  0.0051097653894348 ;
     0.0520964017336447  0.0438876589973454  0.0000000000000000  0.0002238952385848  -0.0066733040364841  -0.0057058488794164  0.0000032362308891  0.0001685056199225  0.1192661172718002  0.0916859011955345  0.0486267815312545  -0.0075343712770952  0.3617659506276415  -0.0128178099611241  0.0856729858763584  0.0235335129837281  -0.0026131563388993  0.0000077859265409  -0.0035937263503887  -0.0071393985379765  0.0009717378792517  0.1822976922708458  0.0475929790540987 ;
     0.0000000000000000  0.0000000000000000  0.0660192858423828  0.0560655196645269  0.0001888482376262  -0.0000326614143342  -0.0056245632385556  -0.0049014754429353  0.0000000000000000  0.0870396684332785  0.0358666094650250  -0.0072271513782604  -0.0158621436959339  0.3490996873570900  0.0026016315116112  0.0582130312267452  -0.0036414244944460  -0.0058604134781389  -0.0029602400251165  -0.0001041259259076  0.1154449206916652  0.0026448070286197  -0.0126963787747857 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0503331683567291  0.0453548387805100  0.0000026351608894  0.0001372088192926  0.0375126774584725  -0.0090335633235269  -0.0047910649778648  0.1045667506033531  0.1928341073699821  0.0136204322253108  0.3549040962617332  -0.0139012513509669  -0.0191112865288007  0.0000846598071647  -0.0021929454227498  -0.0044352495929135  -0.0021442335266286  0.1049735862137911  -0.0142243450622062 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0618777683189294  0.0528238269507302  0.0488234190053139  -0.0117573438404615  0.0939854437254136  0.0800064539506951  0.0532074116918685  0.3055013566412603  -0.0394176552073434  0.3593168556660338  0.0394416506761943  -0.0089438269566795  -0.0043099009051387  0.0002535616183685  0.1864419819103600  -0.0056223457563196  -0.0024882170655625 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0085754028971756  0.0020650737389286  0.0623004022218633  0.1082083078808558  -0.0036116670659104  -0.0116705676444293  -0.0330662326362585  0.0670370588947265  0.5255050882957542  0.0005340010035713  0.0001688836127405  -0.0001906553592150  -0.0101193715300410  0.0817981285826288  0.0022127232158701 ;
     0.0000000000000000  0.0000000000000000  0.0540520196154046  -0.0069829245166739  -0.0000235209268788  0.0000040679582080  0.0452089249242877  -0.0048823293290572  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000130718798568  -0.0083933072502103  0.0000716964018205  -0.0067933901564042  0.0004543881715095  0.3909068326756130  0.0740948865690882  -0.0000012849758665  0.1488626762519327  -0.0000630734143370  -0.0001240757679431 ;
     0.0000000000000000  0.0598393359151736  0.0000000000000000  0.0531907632244693  0.0525563835148440  -0.0090896575973619  0.0008552085571274  0.0445294087530661  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0148195531303397  -0.0107574054841489  -0.0114498186911799  -0.0076524162793178  0.0000097146061809  0.4655674827482490  0.4156800878352041  0.2227617650049469  0.1134193470394420  0.2238510948162243  0.0004445636895373 ;
     0.0636306887676940  -0.0062348281666228  0.0000000000000000  -0.0000318073092025  0.0009480319734900  0.0574862343582599  -0.0000007331887981  -0.0000381760254977  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0080347353961499  -0.0002862533978091  -0.0062484181332741  0.0003393805064464  -0.0002946487674300  0.0000008767723115  0.1462527077565903  0.3419105863941690  -0.0166268174798754  0.0894528381851417  0.0002263608501866 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0003534049411801  -0.0025662058205849  0.0005051515579974  0.0003000906059219  0.0826733000946045  -0.0008271004921458  0.0708035465764230  -0.0043024699743094  0.1904068003412297  0.0555651102347250  -0.0362352869464476  0.1017314926947076  -0.0086173466334150  0.0275907383183825 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  -0.0034913402503150  0.0008407622535105  0.0004459089335228  0.0004316124066331  0.1025796299880008  0.0034443170669199  0.0739451638901269  -0.0038984971932937  0.0633735037021230  -0.0001421892196196  0.0734581597337162  0.1308301577698772  -0.0501143953545952  0.1046434301203381  -0.0021821722985481 ;
     0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0000000000000000  0.0488234190053139  0.0881183600826292  -0.0134657647267512  0.0034014191540561  0.0538550346081265  -0.0333194013142490  -0.0201502375333086  -0.0034645815448976  0.0034475035617883  -0.0005802059055587  0.0001881122610336  0.0007349366299555  0.9942848469082044  0.0135406936231982  0.5093561786892320 ;
   ];
Us{11} =[];
Ucols{11}=[];
Urows{11}=[];
err = 1e-9; [m n] =size(S);
format long g
oldR=[]; c=[];
%Finalizing the permutation
for f=1:11
	npivots(f) = length(cols{f});
	oldR=[oldR, rows{f}(1:npivots(f))]; c=[c, cols{f}];
end
oldR = [oldR setdiff(1:m,oldR)];
newR(oldR)=1:length(oldR);
for f=1:11
	LU(newR(rows{f}),cols{f})=Luf{f};
	LU(newR(Urows{f}),Ucols{f})=Us{f};
end
L=tril(LU,-1)+eye(size(LU));
U=triu(LU); U=U(1:n,:);
spparms('spumoni',3);
[l,u,p]=lu(S);
matlabErr = norm(p*S-l*u)
fprintf('Paru\n');
myErr = norm(S(oldR,c)-L*U)
if(myErr <= 100*matlabErr || myErr<err)
	fprintf('Pass\n')
else
	fprintf('Fail\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
