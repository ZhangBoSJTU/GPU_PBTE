V34 :0x24 disp
8 Disp.f90 S624 0
01/01/2022  20:58:38
use crystal private
enduse
D 64 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 67 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 70 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 73 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 76 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 79 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 133 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 136 23 10 1 11 111 0 0 1 0 0
 0 110 11 11 111 111
D 139 23 10 1 11 111 0 0 1 0 0
 0 110 11 11 111 111
D 142 23 14 2 112 113 0 0 1 0 0
 0 110 11 11 111 111
 0 110 111 11 111 111
D 145 23 14 2 112 113 0 0 1 0 0
 0 110 11 11 111 111
 0 110 111 11 111 111
D 148 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 151 23 14 2 114 118 0 0 1 0 0
 0 116 11 11 117 117
 0 116 117 11 117 117
D 154 23 10 1 11 46 0 0 0 0 0
 0 46 11 11 46 46
D 157 23 14 3 119 124 0 0 1 0 0
 0 121 11 11 122 122
 0 121 122 11 122 122
 0 46 123 11 46 46
S 624 24 0 0 0 9 1 0 5013 10005 0 A 0 0 0 0 B 0 25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 25 0 0 0 0 0 0 disp
S 625 23 5 0 0 0 626 624 5018 0 0 A 0 0 0 0 B 0 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cal_disp
S 626 14 5 0 0 0 1 625 5018 0 400000 A 0 0 0 0 B 0 28 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 28 0 624 0 0 0 0 cal_disp cal_disp 
F 626 0
S 634 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
R 646 6 2 crystal nband
R 667 7 23 crystal e1$ac
R 669 7 25 crystal e2$ac
R 671 7 27 crystal e3$ac
S 718 23 5 0 0 0 724 624 5497 0 0 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cal_disp_point
S 719 7 3 0 0 133 1 718 5512 800004 3000 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 qq
S 720 7 3 0 0 136 1 718 5515 800204 3000 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ff
S 721 7 3 0 0 142 1 718 5518 800204 3000 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 eigenv
S 722 7 3 0 0 139 1 718 5525 800204 3000 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ww
S 723 7 3 0 0 145 1 718 5528 800204 3000 A 0 0 0 0 B 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dm
S 724 14 5 0 0 0 1 718 5497 200 400000 A 0 0 0 0 B 0 45 0 0 0 0 0 3 5 0 0 0 0 0 0 0 0 0 0 0 0 45 0 624 0 0 0 0 cal_disp_point cal_disp_point 
F 724 5 719 720 721 722 723
S 725 6 1 0 0 7 1 718 5531 40800006 3000 A 0 0 0 0 B 0 50 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_402
S 726 6 1 0 0 7 1 718 5539 40800006 3000 A 0 0 0 0 B 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_404
S 727 6 1 0 0 7 1 718 5547 40800006 3000 A 0 0 0 0 B 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_406
S 728 23 5 0 0 0 732 624 5555 0 0 A 0 0 0 0 B 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cal_dm
S 729 7 3 0 0 148 1 728 5512 800004 3000 A 0 0 0 0 B 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 qq
S 730 7 3 0 0 151 1 728 5528 800204 3000 A 0 0 0 0 B 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dm
S 731 6 3 0 0 6 1 728 5562 800004 3000 A 0 0 0 0 B 0 60 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ndim
S 732 14 5 0 0 0 1 728 5555 200 400000 A 0 0 0 0 B 0 60 0 0 0 0 0 9 3 0 0 0 0 0 0 0 0 0 0 0 0 60 0 624 0 0 0 0 cal_dm cal_dm 
F 732 3 729 730 731
S 733 6 1 0 0 7 1 728 5567 40800006 3000 A 0 0 0 0 B 0 66 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_407
S 734 6 1 0 0 7 1 728 5575 40800006 3000 A 0 0 0 0 B 0 66 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_409
S 735 6 1 0 0 7 1 728 5583 40800006 3000 A 0 0 0 0 B 0 66 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_411
S 736 23 5 0 0 0 740 624 5591 0 0 A 0 0 0 0 B 0 86 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cal_ddm
S 737 7 3 0 0 154 1 736 5512 800004 3000 A 0 0 0 0 B 0 86 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 qq
S 738 7 3 0 0 157 1 736 5599 800204 3000 A 0 0 0 0 B 0 86 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ddm
S 739 6 3 0 0 6 1 736 5562 800004 3000 A 0 0 0 0 B 0 86 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ndim
S 740 14 5 0 0 0 1 736 5591 200 400000 A 0 0 0 0 B 0 86 0 0 0 0 0 13 3 0 0 0 0 0 0 0 0 0 0 0 0 86 0 624 0 0 0 0 cal_ddm cal_ddm 
F 740 3 737 738 739
S 741 6 1 0 0 7 1 736 5603 40800006 3000 A 0 0 0 0 B 0 92 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_412
S 742 6 1 0 0 7 1 736 5611 40800006 3000 A 0 0 0 0 B 0 92 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_414
S 743 6 1 0 0 7 1 736 5619 40800006 3000 A 0 0 0 0 B 0 92 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_416
S 744 6 1 0 0 7 1 736 5627 40800006 3000 A 0 0 0 0 B 0 92 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_419
A 12 2 0 0 0 10 618 0 0 0 12 0 0 0 0 0 0 0 0 0 0 0
A 13 2 0 0 0 10 617 0 0 0 13 0 0 0 0 0 0 0 0 0 0 0
A 46 2 0 0 0 7 634 0 0 0 46 0 0 0 0 0 0 0 0 0 0 0
A 92 1 0 5 0 64 667 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 97 1 0 5 0 70 669 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 102 1 0 5 0 76 671 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 109 1 0 0 0 6 646 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 110 7 0 0 0 7 109 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 111 1 0 0 0 7 725 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 112 1 0 0 0 7 727 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 113 1 0 0 0 7 726 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 114 1 0 0 0 7 735 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 115 1 0 0 0 6 731 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 116 7 0 0 0 7 115 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 117 1 0 0 0 7 733 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 118 1 0 0 0 7 734 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 119 1 0 0 44 7 744 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 120 1 0 0 0 6 739 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 121 7 0 0 0 7 120 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 122 1 0 0 36 7 741 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 123 1 0 0 38 7 742 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 124 1 0 0 42 7 743 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 30 1 1
V 92 64 7 0
R 0 67 0 0
A 0 10 0 0 1 12 1
A 0 10 0 0 1 13 1
A 0 10 0 0 1 13 0
J 31 1 1
V 97 70 7 0
R 0 73 0 0
A 0 10 0 0 1 13 1
A 0 10 0 0 1 12 1
A 0 10 0 0 1 13 0
J 32 1 1
V 102 76 7 0
R 0 79 0 0
A 0 10 0 0 1 13 1
A 0 10 0 0 1 13 1
A 0 10 0 0 1 12 0
Z
