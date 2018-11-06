%%%%%%%% START %%%%%%%%%%%%%%%%%%% 
clear all
% A  is  66 x 66 
LU = zeros(66,66);
npivots =[]; 
S = zeros(66,66);
% nf=66
% anz = 66  rjsize=66
S(1,1)= 0.0921385805100000;
S(2,2)= 0.0921385805100000;
S(3,3)= 0.0921385805100000;
S(4,4)= 0.1379957379830000;
S(5,5)= 0.1379957379830000;
S(6,6)= 0.1379957379830000;
S(7,7)= 0.1379957379830000;
S(8,8)= 0.1379957379830000;
S(9,9)= 0.1379957379830000;
S(10,10)= 0.0921385805100000;
S(11,11)= 0.0921385805100000;
S(12,12)= 0.0921385805100000;
S(13,13)= 0.1728285734550000;
S(14,14)= 0.1728285734550000;
S(15,15)= 0.1728285734550000;
S(16,16)= 0.0852383576022000;
S(17,17)= 0.0852383576022000;
S(18,18)= 0.0852383576022000;
S(19,19)= 0.0852383576022000;
S(20,20)= 0.0852383576022000;
S(21,21)= 0.0852383576022000;
S(22,22)= 0.1728285734550000;
S(23,23)= 0.1728285734550000;
S(24,24)= 0.1728285734550000;
S(25,25)= 0.0617332189107000;
S(26,26)= 0.0617332189107000;
S(27,27)= 0.0617332189107000;
S(28,28)= 0.1413083414760000;
S(29,29)= 0.1413083414760000;
S(30,30)= 0.1413083414760000;
S(31,31)= 0.1413083414760000;
S(32,32)= 0.1413083414760000;
S(33,33)= 0.1413083414760000;
S(34,34)= 0.0617332189107000;
S(35,35)= 0.0617332189107000;
S(36,36)= 0.0617332189107000;
S(37,37)= 0.1254266384520000;
S(38,38)= 0.1254266384520000;
S(39,39)= 0.1254266384520000;
S(40,40)= 0.0533208927371000;
S(41,41)= 0.0533208927371000;
S(42,42)= 0.0533208927371000;
S(43,43)= 0.0533208927371000;
S(44,44)= 0.0533208927371000;
S(45,45)= 0.0533208927371000;
S(46,46)= 0.1254266384520000;
S(47,47)= 0.1254266384520000;
S(48,48)= 0.1254266384520000;
S(49,49)= 0.0231706100487000;
S(50,50)= 0.0231706100487000;
S(51,51)= 0.0231706100487000;
S(52,52)= 0.0305931884075000;
S(53,53)= 0.0305931884075000;
S(54,54)= 0.0305931884075000;
S(55,55)= 0.0648568699369000;
S(56,56)= 0.0648568699369000;
S(57,57)= 0.0648568699369000;
S(58,58)= 0.0648568699369000;
S(59,59)= 0.0648568699369000;
S(60,60)= 0.0648568699369000;
S(61,61)= 0.0305931884075000;
S(62,62)= 0.0305931884075000;
S(63,63)= 0.0305931884075000;
S(64,64)= 0.0197469386650000;
S(65,65)= 0.0197469386650000;
S(66,66)= 0.0197469386650000;
%~~~~~~~  Assemble Front 0 start ~~~~
% fp=1 pivotal columns:clo1=0...col2=0
%L part:
cols{1} = [1 ];
rows{1} = [1 ];
Luf{1}= [  0.0921385805100000 ;
   ];
U{1} =[];
Ucols{1}=[];
Urows{1}=[];
%~~~~~~~  Assemble Front 1 start ~~~~
% fp=1 pivotal columns:clo1=1...col2=1
%L part:
cols{2} = [2 ];
rows{2} = [2 ];
Luf{2}= [  0.0921385805100000 ;
   ];
U{2} =[];
Ucols{2}=[];
Urows{2}=[];
%~~~~~~~  Assemble Front 2 start ~~~~
% fp=1 pivotal columns:clo1=2...col2=2
%L part:
cols{3} = [3 ];
rows{3} = [3 ];
Luf{3}= [  0.0921385805100000 ;
   ];
U{3} =[];
Ucols{3}=[];
Urows{3}=[];
%~~~~~~~  Assemble Front 3 start ~~~~
% fp=1 pivotal columns:clo1=3...col2=3
%L part:
cols{4} = [4 ];
rows{4} = [4 ];
Luf{4}= [  0.1379957379830000 ;
   ];
U{4} =[];
Ucols{4}=[];
Urows{4}=[];
%~~~~~~~  Assemble Front 4 start ~~~~
% fp=1 pivotal columns:clo1=4...col2=4
%L part:
cols{5} = [5 ];
rows{5} = [5 ];
Luf{5}= [  0.1379957379830000 ;
   ];
U{5} =[];
Ucols{5}=[];
Urows{5}=[];
%~~~~~~~  Assemble Front 5 start ~~~~
% fp=1 pivotal columns:clo1=5...col2=5
%L part:
cols{6} = [6 ];
rows{6} = [6 ];
Luf{6}= [  0.1379957379830000 ;
   ];
U{6} =[];
Ucols{6}=[];
Urows{6}=[];
%~~~~~~~  Assemble Front 6 start ~~~~
% fp=1 pivotal columns:clo1=6...col2=6
%L part:
cols{7} = [7 ];
rows{7} = [7 ];
Luf{7}= [  0.1379957379830000 ;
   ];
U{7} =[];
Ucols{7}=[];
Urows{7}=[];
%~~~~~~~  Assemble Front 7 start ~~~~
% fp=1 pivotal columns:clo1=7...col2=7
%L part:
cols{8} = [8 ];
rows{8} = [8 ];
Luf{8}= [  0.1379957379830000 ;
   ];
U{8} =[];
Ucols{8}=[];
Urows{8}=[];
%~~~~~~~  Assemble Front 8 start ~~~~
% fp=1 pivotal columns:clo1=8...col2=8
%L part:
cols{9} = [9 ];
rows{9} = [9 ];
Luf{9}= [  0.1379957379830000 ;
   ];
U{9} =[];
Ucols{9}=[];
Urows{9}=[];
%~~~~~~~  Assemble Front 9 start ~~~~
% fp=1 pivotal columns:clo1=9...col2=9
%L part:
cols{10} = [10 ];
rows{10} = [10 ];
Luf{10}= [  0.0921385805100000 ;
   ];
U{10} =[];
Ucols{10}=[];
Urows{10}=[];
%~~~~~~~  Assemble Front 10 start ~~~~
% fp=1 pivotal columns:clo1=10...col2=10
%L part:
cols{11} = [11 ];
rows{11} = [11 ];
Luf{11}= [  0.0921385805100000 ;
   ];
U{11} =[];
Ucols{11}=[];
Urows{11}=[];
%~~~~~~~  Assemble Front 11 start ~~~~
% fp=1 pivotal columns:clo1=11...col2=11
%L part:
cols{12} = [12 ];
rows{12} = [12 ];
Luf{12}= [  0.0921385805100000 ;
   ];
U{12} =[];
Ucols{12}=[];
Urows{12}=[];
%~~~~~~~  Assemble Front 12 start ~~~~
% fp=1 pivotal columns:clo1=12...col2=12
%L part:
cols{13} = [13 ];
rows{13} = [13 ];
Luf{13}= [  0.1728285734550000 ;
   ];
U{13} =[];
Ucols{13}=[];
Urows{13}=[];
%~~~~~~~  Assemble Front 13 start ~~~~
% fp=1 pivotal columns:clo1=13...col2=13
%L part:
cols{14} = [14 ];
rows{14} = [14 ];
Luf{14}= [  0.1728285734550000 ;
   ];
U{14} =[];
Ucols{14}=[];
Urows{14}=[];
%~~~~~~~  Assemble Front 14 start ~~~~
% fp=1 pivotal columns:clo1=14...col2=14
%L part:
cols{15} = [15 ];
rows{15} = [15 ];
Luf{15}= [  0.1728285734550000 ;
   ];
U{15} =[];
Ucols{15}=[];
Urows{15}=[];
%~~~~~~~  Assemble Front 15 start ~~~~
% fp=1 pivotal columns:clo1=15...col2=15
%L part:
cols{16} = [16 ];
rows{16} = [16 ];
Luf{16}= [  0.0852383576022000 ;
   ];
U{16} =[];
Ucols{16}=[];
Urows{16}=[];
%~~~~~~~  Assemble Front 16 start ~~~~
% fp=1 pivotal columns:clo1=16...col2=16
%L part:
cols{17} = [17 ];
rows{17} = [17 ];
Luf{17}= [  0.0852383576022000 ;
   ];
U{17} =[];
Ucols{17}=[];
Urows{17}=[];
%~~~~~~~  Assemble Front 17 start ~~~~
% fp=1 pivotal columns:clo1=17...col2=17
%L part:
cols{18} = [18 ];
rows{18} = [18 ];
Luf{18}= [  0.0852383576022000 ;
   ];
U{18} =[];
Ucols{18}=[];
Urows{18}=[];
%~~~~~~~  Assemble Front 18 start ~~~~
% fp=1 pivotal columns:clo1=18...col2=18
%L part:
cols{19} = [19 ];
rows{19} = [19 ];
Luf{19}= [  0.0852383576022000 ;
   ];
U{19} =[];
Ucols{19}=[];
Urows{19}=[];
%~~~~~~~  Assemble Front 19 start ~~~~
% fp=1 pivotal columns:clo1=19...col2=19
%L part:
cols{20} = [20 ];
rows{20} = [20 ];
Luf{20}= [  0.0852383576022000 ;
   ];
U{20} =[];
Ucols{20}=[];
Urows{20}=[];
%~~~~~~~  Assemble Front 20 start ~~~~
% fp=1 pivotal columns:clo1=20...col2=20
%L part:
cols{21} = [21 ];
rows{21} = [21 ];
Luf{21}= [  0.0852383576022000 ;
   ];
U{21} =[];
Ucols{21}=[];
Urows{21}=[];
%~~~~~~~  Assemble Front 21 start ~~~~
% fp=1 pivotal columns:clo1=21...col2=21
%L part:
cols{22} = [22 ];
rows{22} = [22 ];
Luf{22}= [  0.1728285734550000 ;
   ];
U{22} =[];
Ucols{22}=[];
Urows{22}=[];
%~~~~~~~  Assemble Front 22 start ~~~~
% fp=1 pivotal columns:clo1=22...col2=22
%L part:
cols{23} = [23 ];
rows{23} = [23 ];
Luf{23}= [  0.1728285734550000 ;
   ];
U{23} =[];
Ucols{23}=[];
Urows{23}=[];
%~~~~~~~  Assemble Front 23 start ~~~~
% fp=1 pivotal columns:clo1=23...col2=23
%L part:
cols{24} = [24 ];
rows{24} = [24 ];
Luf{24}= [  0.1728285734550000 ;
   ];
U{24} =[];
Ucols{24}=[];
Urows{24}=[];
%~~~~~~~  Assemble Front 24 start ~~~~
% fp=1 pivotal columns:clo1=24...col2=24
%L part:
cols{25} = [25 ];
rows{25} = [25 ];
Luf{25}= [  0.0617332189107000 ;
   ];
U{25} =[];
Ucols{25}=[];
Urows{25}=[];
%~~~~~~~  Assemble Front 25 start ~~~~
% fp=1 pivotal columns:clo1=25...col2=25
%L part:
cols{26} = [26 ];
rows{26} = [26 ];
Luf{26}= [  0.0617332189107000 ;
   ];
U{26} =[];
Ucols{26}=[];
Urows{26}=[];
%~~~~~~~  Assemble Front 26 start ~~~~
% fp=1 pivotal columns:clo1=26...col2=26
%L part:
cols{27} = [27 ];
rows{27} = [27 ];
Luf{27}= [  0.0617332189107000 ;
   ];
U{27} =[];
Ucols{27}=[];
Urows{27}=[];
%~~~~~~~  Assemble Front 27 start ~~~~
% fp=1 pivotal columns:clo1=27...col2=27
%L part:
cols{28} = [28 ];
rows{28} = [28 ];
Luf{28}= [  0.1413083414760000 ;
   ];
U{28} =[];
Ucols{28}=[];
Urows{28}=[];
%~~~~~~~  Assemble Front 28 start ~~~~
% fp=1 pivotal columns:clo1=28...col2=28
%L part:
cols{29} = [29 ];
rows{29} = [29 ];
Luf{29}= [  0.1413083414760000 ;
   ];
U{29} =[];
Ucols{29}=[];
Urows{29}=[];
%~~~~~~~  Assemble Front 29 start ~~~~
% fp=1 pivotal columns:clo1=29...col2=29
%L part:
cols{30} = [30 ];
rows{30} = [30 ];
Luf{30}= [  0.1413083414760000 ;
   ];
U{30} =[];
Ucols{30}=[];
Urows{30}=[];
%~~~~~~~  Assemble Front 30 start ~~~~
% fp=1 pivotal columns:clo1=30...col2=30
%L part:
cols{31} = [31 ];
rows{31} = [31 ];
Luf{31}= [  0.1413083414760000 ;
   ];
U{31} =[];
Ucols{31}=[];
Urows{31}=[];
%~~~~~~~  Assemble Front 31 start ~~~~
% fp=1 pivotal columns:clo1=31...col2=31
%L part:
cols{32} = [32 ];
rows{32} = [32 ];
Luf{32}= [  0.1413083414760000 ;
   ];
U{32} =[];
Ucols{32}=[];
Urows{32}=[];
%~~~~~~~  Assemble Front 32 start ~~~~
% fp=1 pivotal columns:clo1=32...col2=32
%L part:
cols{33} = [33 ];
rows{33} = [33 ];
Luf{33}= [  0.1413083414760000 ;
   ];
U{33} =[];
Ucols{33}=[];
Urows{33}=[];
%~~~~~~~  Assemble Front 33 start ~~~~
% fp=1 pivotal columns:clo1=33...col2=33
%L part:
cols{34} = [34 ];
rows{34} = [34 ];
Luf{34}= [  0.0617332189107000 ;
   ];
U{34} =[];
Ucols{34}=[];
Urows{34}=[];
%~~~~~~~  Assemble Front 34 start ~~~~
% fp=1 pivotal columns:clo1=34...col2=34
%L part:
cols{35} = [35 ];
rows{35} = [35 ];
Luf{35}= [  0.0617332189107000 ;
   ];
U{35} =[];
Ucols{35}=[];
Urows{35}=[];
%~~~~~~~  Assemble Front 35 start ~~~~
% fp=1 pivotal columns:clo1=35...col2=35
%L part:
cols{36} = [36 ];
rows{36} = [36 ];
Luf{36}= [  0.0617332189107000 ;
   ];
U{36} =[];
Ucols{36}=[];
Urows{36}=[];
%~~~~~~~  Assemble Front 36 start ~~~~
% fp=1 pivotal columns:clo1=36...col2=36
%L part:
cols{37} = [37 ];
rows{37} = [37 ];
Luf{37}= [  0.1254266384520000 ;
   ];
U{37} =[];
Ucols{37}=[];
Urows{37}=[];
%~~~~~~~  Assemble Front 37 start ~~~~
% fp=1 pivotal columns:clo1=37...col2=37
%L part:
cols{38} = [38 ];
rows{38} = [38 ];
Luf{38}= [  0.1254266384520000 ;
   ];
U{38} =[];
Ucols{38}=[];
Urows{38}=[];
%~~~~~~~  Assemble Front 38 start ~~~~
% fp=1 pivotal columns:clo1=38...col2=38
%L part:
cols{39} = [39 ];
rows{39} = [39 ];
Luf{39}= [  0.1254266384520000 ;
   ];
U{39} =[];
Ucols{39}=[];
Urows{39}=[];
%~~~~~~~  Assemble Front 39 start ~~~~
% fp=1 pivotal columns:clo1=39...col2=39
%L part:
cols{40} = [40 ];
rows{40} = [40 ];
Luf{40}= [  0.0533208927371000 ;
   ];
U{40} =[];
Ucols{40}=[];
Urows{40}=[];
%~~~~~~~  Assemble Front 40 start ~~~~
% fp=1 pivotal columns:clo1=40...col2=40
%L part:
cols{41} = [41 ];
rows{41} = [41 ];
Luf{41}= [  0.0533208927371000 ;
   ];
U{41} =[];
Ucols{41}=[];
Urows{41}=[];
%~~~~~~~  Assemble Front 41 start ~~~~
% fp=1 pivotal columns:clo1=41...col2=41
%L part:
cols{42} = [42 ];
rows{42} = [42 ];
Luf{42}= [  0.0533208927371000 ;
   ];
U{42} =[];
Ucols{42}=[];
Urows{42}=[];
%~~~~~~~  Assemble Front 42 start ~~~~
% fp=1 pivotal columns:clo1=42...col2=42
%L part:
cols{43} = [43 ];
rows{43} = [43 ];
Luf{43}= [  0.0533208927371000 ;
   ];
U{43} =[];
Ucols{43}=[];
Urows{43}=[];
%~~~~~~~  Assemble Front 43 start ~~~~
% fp=1 pivotal columns:clo1=43...col2=43
%L part:
cols{44} = [44 ];
rows{44} = [44 ];
Luf{44}= [  0.0533208927371000 ;
   ];
U{44} =[];
Ucols{44}=[];
Urows{44}=[];
%~~~~~~~  Assemble Front 44 start ~~~~
% fp=1 pivotal columns:clo1=44...col2=44
%L part:
cols{45} = [45 ];
rows{45} = [45 ];
Luf{45}= [  0.0533208927371000 ;
   ];
U{45} =[];
Ucols{45}=[];
Urows{45}=[];
%~~~~~~~  Assemble Front 45 start ~~~~
% fp=1 pivotal columns:clo1=45...col2=45
%L part:
cols{46} = [46 ];
rows{46} = [46 ];
Luf{46}= [  0.1254266384520000 ;
   ];
U{46} =[];
Ucols{46}=[];
Urows{46}=[];
%~~~~~~~  Assemble Front 46 start ~~~~
% fp=1 pivotal columns:clo1=46...col2=46
%L part:
cols{47} = [47 ];
rows{47} = [47 ];
Luf{47}= [  0.1254266384520000 ;
   ];
U{47} =[];
Ucols{47}=[];
Urows{47}=[];
%~~~~~~~  Assemble Front 47 start ~~~~
% fp=1 pivotal columns:clo1=47...col2=47
%L part:
cols{48} = [48 ];
rows{48} = [48 ];
Luf{48}= [  0.1254266384520000 ;
   ];
U{48} =[];
Ucols{48}=[];
Urows{48}=[];
%~~~~~~~  Assemble Front 48 start ~~~~
% fp=1 pivotal columns:clo1=48...col2=48
%L part:
cols{49} = [49 ];
rows{49} = [49 ];
Luf{49}= [  0.0231706100487000 ;
   ];
U{49} =[];
Ucols{49}=[];
Urows{49}=[];
%~~~~~~~  Assemble Front 49 start ~~~~
% fp=1 pivotal columns:clo1=49...col2=49
%L part:
cols{50} = [50 ];
rows{50} = [50 ];
Luf{50}= [  0.0231706100487000 ;
   ];
U{50} =[];
Ucols{50}=[];
Urows{50}=[];
%~~~~~~~  Assemble Front 50 start ~~~~
% fp=1 pivotal columns:clo1=50...col2=50
%L part:
cols{51} = [51 ];
rows{51} = [51 ];
Luf{51}= [  0.0231706100487000 ;
   ];
U{51} =[];
Ucols{51}=[];
Urows{51}=[];
%~~~~~~~  Assemble Front 51 start ~~~~
% fp=1 pivotal columns:clo1=51...col2=51
%L part:
cols{52} = [52 ];
rows{52} = [52 ];
Luf{52}= [  0.0305931884075000 ;
   ];
U{52} =[];
Ucols{52}=[];
Urows{52}=[];
%~~~~~~~  Assemble Front 52 start ~~~~
% fp=1 pivotal columns:clo1=52...col2=52
%L part:
cols{53} = [53 ];
rows{53} = [53 ];
Luf{53}= [  0.0305931884075000 ;
   ];
U{53} =[];
Ucols{53}=[];
Urows{53}=[];
%~~~~~~~  Assemble Front 53 start ~~~~
% fp=1 pivotal columns:clo1=53...col2=53
%L part:
cols{54} = [54 ];
rows{54} = [54 ];
Luf{54}= [  0.0305931884075000 ;
   ];
U{54} =[];
Ucols{54}=[];
Urows{54}=[];
%~~~~~~~  Assemble Front 54 start ~~~~
% fp=1 pivotal columns:clo1=54...col2=54
%L part:
cols{55} = [55 ];
rows{55} = [55 ];
Luf{55}= [  0.0648568699369000 ;
   ];
U{55} =[];
Ucols{55}=[];
Urows{55}=[];
%~~~~~~~  Assemble Front 55 start ~~~~
% fp=1 pivotal columns:clo1=55...col2=55
%L part:
cols{56} = [56 ];
rows{56} = [56 ];
Luf{56}= [  0.0648568699369000 ;
   ];
U{56} =[];
Ucols{56}=[];
Urows{56}=[];
%~~~~~~~  Assemble Front 56 start ~~~~
% fp=1 pivotal columns:clo1=56...col2=56
%L part:
cols{57} = [57 ];
rows{57} = [57 ];
Luf{57}= [  0.0648568699369000 ;
   ];
U{57} =[];
Ucols{57}=[];
Urows{57}=[];
%~~~~~~~  Assemble Front 57 start ~~~~
% fp=1 pivotal columns:clo1=57...col2=57
%L part:
cols{58} = [58 ];
rows{58} = [58 ];
Luf{58}= [  0.0648568699369000 ;
   ];
U{58} =[];
Ucols{58}=[];
Urows{58}=[];
%~~~~~~~  Assemble Front 58 start ~~~~
% fp=1 pivotal columns:clo1=58...col2=58
%L part:
cols{59} = [59 ];
rows{59} = [59 ];
Luf{59}= [  0.0648568699369000 ;
   ];
U{59} =[];
Ucols{59}=[];
Urows{59}=[];
%~~~~~~~  Assemble Front 59 start ~~~~
% fp=1 pivotal columns:clo1=59...col2=59
%L part:
cols{60} = [60 ];
rows{60} = [60 ];
Luf{60}= [  0.0648568699369000 ;
   ];
U{60} =[];
Ucols{60}=[];
Urows{60}=[];
%~~~~~~~  Assemble Front 60 start ~~~~
% fp=1 pivotal columns:clo1=60...col2=60
%L part:
cols{61} = [61 ];
rows{61} = [61 ];
Luf{61}= [  0.0305931884075000 ;
   ];
U{61} =[];
Ucols{61}=[];
Urows{61}=[];
%~~~~~~~  Assemble Front 61 start ~~~~
% fp=1 pivotal columns:clo1=61...col2=61
%L part:
cols{62} = [62 ];
rows{62} = [62 ];
Luf{62}= [  0.0305931884075000 ;
   ];
U{62} =[];
Ucols{62}=[];
Urows{62}=[];
%~~~~~~~  Assemble Front 62 start ~~~~
% fp=1 pivotal columns:clo1=62...col2=62
%L part:
cols{63} = [63 ];
rows{63} = [63 ];
Luf{63}= [  0.0305931884075000 ;
   ];
U{63} =[];
Ucols{63}=[];
Urows{63}=[];
%~~~~~~~  Assemble Front 63 start ~~~~
% fp=1 pivotal columns:clo1=63...col2=63
%L part:
cols{64} = [64 ];
rows{64} = [64 ];
Luf{64}= [  0.0197469386650000 ;
   ];
U{64} =[];
Ucols{64}=[];
Urows{64}=[];
%~~~~~~~  Assemble Front 64 start ~~~~
% fp=1 pivotal columns:clo1=64...col2=64
%L part:
cols{65} = [65 ];
rows{65} = [65 ];
Luf{65}= [  0.0197469386650000 ;
   ];
U{65} =[];
Ucols{65}=[];
Urows{65}=[];
%~~~~~~~  Assemble Front 65 start ~~~~
% fp=1 pivotal columns:clo1=65...col2=65
%L part:
cols{66} = [66 ];
rows{66} = [66 ];
Luf{66}= [  0.0197469386650000 ;
   ];
U{66} =[];
Ucols{66}=[];
Urows{66}=[];
err = 1e-12; [m n] =size(S);
oldR=[]; c=[];
%Finalizing the permutation
for f=1:66
	npivots(f) = length(cols{f});
	oldR=[oldR, rows{f}(1:npivots(f))]; c=[c, cols{f}];
end
oldR = [oldR setdiff(1:m,oldR)];
newR(oldR)=1:length(oldR);
for f=1:66
	LU(newR(rows{f}),cols{f})=Luf{f};
	LU(newR(Urows{f}),Ucols{f})=U{f};
end
L=tril(LU,-1)+eye(size(LU));
U=triu(LU); U=U(1:n,:);
norm(S(oldR,c)-L*U)
if( (norm(S(oldR,c)-L*U)) < err )
	fprintf('Pass\n')
else
	fprintf('Fail\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*

