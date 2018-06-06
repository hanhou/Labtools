function DrawMapping(monkey_hemi)
% DrawMapping.m
% First version by HH @ Gu Lab 2013

monkey_hemis = {'Polo_L','Polo_R','Messi_L','Messi_R'};

if nargin == 0
    % monkey_hemi = 'Polo_R'; % AP 0:5
    % monkey_hemi = 'Polo_L'; % AP -3:3.2
    % monkey_hemi = 'Messi_L'; % AP -4:0.8
    monkey_hemi = 3; % 
end;

%% Parameters
clc; clear global
global maxX maxY maxZ movDis GMTypes data;
global hemisphere monkey; global gridRange;
global MRI_path MRI_offset AP0; global MRI_offset_new;
global linWid start_end_markers overlapping Polo_right_AP0;
global overlapTransparent; 

linWid = 1.3;
% overlapping = [-10 10; -10 10]; start_end_markers = false; % First row for area annotation; second for unit annotation
% overlapping = [-1 0; -1 0]; start_end_markers = false; % First row for area annotation; second for unit annotation
overlapping = [0 0; 0 0]; start_end_markers =0; % First row for area annotation; second for unit annotation
overlapTransparent = 1;

maxX = 30; % Grid size
maxY = 30;
maxZ = 200;   % Drawing range (in 100 um)
movDis = 150; % Maximum travel distance of the microdriver (in 100 um)

V = ones(30,25,maxZ + 30,3);  % Matrix to be rendered. The z-resolution = 100 um.
Polo_right_AP0 = 9; % For MRI alignment. This controls MRI-AP alignment so will affect all monkeys, whereas "AP0" of each monkey below only affects AP-Grid alignment.  HH20150918

%% Define our area types here
GMTypes = { % 'Type', ColorCode (RGB);
    'GM',[1 1 1];    % Gray matter without visual modulation
    'VM',[0 0 0];          % Visual modulation but not sure to which area it belongs
    'LIP',[0.6 0.8 0.4];
    'VIP',[1 1 0];
    'MST',[0 0 1];
    'MT',[0.8 0 0];
    'MIP',[0.6 0.6 0.6];  % [0 0.6 0]
    'AUD',[1 0 1];
    };
for i = 1:length(GMTypes)
    eval([GMTypes{i,1} '= -' num2str(i) ';']);  % This is a trick.
end

%% Input our mapping data here

switch monkey_hemis{monkey_hemi}
    
    case 'Messi_R'
        
        %  Messi_right
        %  %{
        
        % Header
        toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
        toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
        gridRange = [1 20; 1 25];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
        
        monkey = 10;
        hemisphere = 2;
        AP0 = 10 % 6;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , vode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {8,[11,10],[1.9 0.0],[GM 8 17; GM 58 72; VIP 95 118]}
            {40,[13,18],[2.2 0.0-0.3],[GM 31 60; MST 75 103;  MT 120 147]}
            {40,[13,8],[2.05 0.0],[GM 37 59; VIP 67 86]}
            {41,[13,13],[2.05 0.2],[GM 7 18; LIP 31 43], 43}
            {42,[13,12],[1.9 0],[GM 13 56; LIP 74 84],84}
            {43,[13,14],[2.05 0],[GM 8 28; LIP 38 70; MST 92 107; MST 117 135]}
            {44,[13,15],[2.2 0],[GM 7 13; LIP 23 45; MST 69 90; MT 100 130]}
            {45,[13,11],[1.9 0],[GM 5 17; GM 22 63; LIP 73 96], 150}
            {45,[13,2],[1.9 0],[GM 9 49]}
            {46,[11,14],[2.1 0],[GM 3 37; LIP 51 63],63}
            {47,[9,15],[2.2 0+0.1],[GM 10 20; LIP 32 63; MST 114 130]}
            {48,[7,17],[2.2 0],[GM 12 20; LIP 32 63; GM 82 112; GM 132 146]}
            {49,[9,15],[2.2 0],[GM 4 23; LIP 35 40],40}
            {50,[10,14],[2.2 0],[GM 0 22; LIP 46 60],60}
            {53,[13,12],[1.9 0],[GM 13 58; LIP 73 85],85}
            {54,[12,13],[2.05 0],[GM 4 34; LIP 52 67],67}
            {55,[12,12],[2.05 0],[GM 14 51; LIP 69 75],75}
            {56,[12,11],[2.1 0],[GM 10 57; LIP 65 78],78}
            {57,[13,14],[2.05 0.1],[GM 0 22; LIP 33 60],60}
            {58,[13,15],[2.2 0],[LIP 20 39; MST 63 70],70}
            {72,[14,12],[2.05 0],[GM 2 30; LIP 51 75]}
            {73,[14,13],[1.9 0],[GM 13 34; LIP 50 67],67}
            {74,[14,11],[2.05 0],[GM 5 47; LIP 62 70],70}
            {75,[14,14],[1.9+0.15 0],[LIP 33 66; MST 92 95],95}
            {76,[15,10],[1.9 0],[GM 27 62; LIP 75 93],93}
            {77,[15,11],[2.05-0.2 0],[GM 27 59; LIP 72 83],83}
            {78,[15,12],[2.05+0.1 0],[GM 0 20; LIP 34 57],57}
            {79,[16,10],[2.05 0],[GM 2 32; LIP 43 75]}
            }';
        
        MRI_path = 'Z:\Data\MOOG\Polo\Mapping\MRI\PoloOutput\forDrawMapping\';
%         MRI_offset = {[-68 79]*1.1-9.5 ,[-240 500]*1.1+23 , [0 -2.5]};   % [x1 x2],[y1 y2], [dx/dxSelect slope, dy/dxSelect slope]
        
        MRI_offset = {[-70.9547, 67.0547], [-233.02, 565.02],[0, -2.5]}; % Manual adjust HH20180604
        
        % Update version with ratio fixed {[x_center, y_center, x_range], [slopes]}. HH20180606
        MRI_offset_new = {[-1.9500,  166,  138.0094 ], [0, -2.5]};

        %}

    case 'Messi_L'
        %  Messi_left
        %  %{
        
        % Header
        toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
        toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
        gridRange = [0 20; 0 21];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
        
        monkey = 10;
        hemisphere = 1;
        AP0 = 10 % 6;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {1,[15,15],[2.2 0.0],[GM 4 37; GM 56 74; MST 88 117]}
            {2,[15,13],[2.1 0.0],[GM 5 18; GM 31 52; MST 74 96; MST 106 124]}
            {3,[15,17],[2.4 0.0],[GM 26 56; MT 69 122]}
            {3,[15,11],[2.1 0.0 + 0.1],[GM 0 20; GM 29 47; MST 75 96;]}
            {3,[15,9],[2.1 0.0 + 0.1],[GM 6 36; GM 43 67; GM 116 149]}
            {4,[15,19],[2.25 0.0],[GM 7 79],79}
            {4,[13,20],[2.25 0.0],[GM 9 39; MST 89 135; MT 142 150]}
            {4,[13,13],[2.1 0.0],[GM 7 37; GM 44 72; MST 92 122; MT 134 150]}
            {5,[13,17],[2.25 0.0],[GM 5 50; MST 66 90; MT 109 140]}
            {6,[15,7],[2.05 0.0],[GM 33 59; GM 71 90; LIP 108 140]}
            {6,[11,5],[1.8 0.0],[GM 1 26; GM 39 58; GM 66 84; VIP 104 121; VIP 129 143]}
            {7,[11,9],[1.9 0.0],[GM 2 30; GM 65 97; LIP 105 125],125}
            {8,[11,10],[1.9 0.0],[GM 4 33; GM 55 85; LIP 97 123]}
            {9,[9,12],[2.0 0-0.2],[GM 18 39; GM 69 98; LIP 108 128]}
            {10,[9,11],[2.0 0],[GM 7 35;GM 67 95;LIP 106 121],121}
            {11,[13,9],[1.9 0+0.1],[GM 44 73; LIP 82 104]}
            {12,[13,8],[1.9 0],[GM 58 86; LIP 93 102],102}
            {13,[13,7],[1.8 0],[GM 0 13; GM 29 50; GM 70 96; LIP 104 116],116}
            {14,[13,6],[1.7 0],[GM 10 24; GM 37 49; GM 55 70; GM 89 111; LIP 117 129],129}
            {15,[12,8],[1.8 0.0 + 0.1],[GM 6 19; GM 64 93; LIP 102 116],116}
            {16,[12,7],[1.8 0.0 + 0.1],[GM 2 9; GM 31 46; GM 74 97; LIP 105 114; VIP 114 128]}
            {17,[12,9],[1.9 0.0 + 0.1],[GM 0 7; GM 49 79; LIP 86 90],90}
            {18,[12,10],[2.0 0 + 0.1],[GM 3 12; GM 37 64; LIP 73 80],80}
            {19,[12,11],[2.0 0 + 0.1],[GM 3 10; GM 26 55; LIP 64 88],96}
            {20,[11,11],[2.0 0 - 0.05],[GM 4 19; GM 44 79; LIP 88 105],105}
            {21,[14,7],[1.8 0 + 0.1],[GM 20 34; GM 56 81; LIP 89 101],101}
            {22,[14,8],[1.9 0 + 0.1],[GM 42 69; LIP 80 101],120}
            {23,[14,9],[1.9 0 + 0.1],[GM 36 64; LIP 72 79],79}
            {24,[10,12],[2.0 0-0.1],[GM 9 26; GM 48 82;LIP 94 115]}
            {24,[12,10],[2.0 0+0.1],[GM 3 12; GM 36 62; LIP 72 77],77}
            {25,[12,10],[2.0 0+0.1],[GM 36 62; LIP 72 86],86}
            {26,[12,11],[2.0 0+0.1],[GM 30 58; LIP 67 75],75}
            {27,[18,5],[1.7 0.4+0.05],[GM 11 36; LIP 47 65; GM 77 100; GM 109 130]}
            {28,[18,4],[1.7 0.4],[GM 18 42; LIP 53 61; GM 71 91; GM 103 125]}
            {29,[12,12],[2.0 0],[GM 9 18; GM 29 58; LIP 66 74],74}
            {30,[12,13],[2.1 0],[GM 19 39; LIP 50 60],60}
            {31,[11,13],[2.1 0],[GM 0 11; GM 22 55; LIP 65 80],80}
            {32,[11,12],[2.1 0],[GM 0 7; GM 26 57; LIP 69 82],82}
            {33,[11,11],[2.0 0.1],[GM 0 10; GM 34 61; LIP 76 97]}
            {34,[13,11],[2.0 0.1],[GM 26 52; LIP 62 80]}
            {35,[13,10],[2.0 0],[GM 35 65; LIP 75 96],96}
            {36,[14,10],[2.0 0+0.1],[GM 20 46; LIP 55 80],90}
            {37,[13,12],[2.0 0+0.1],[GM 3 8; GM 20 47; LIP 56 82]}
            {38,[14,11],[2.0 0+0.1],[GM 0 38; LIP 48 65]}
            {39,[15,10],[2.0 0],[GM 22 48; LIP 59 81; GM 120 150]}
            {59,[13,11],[2.05 0.0 + 0.1],[GM 13 39; LIP 56 74],74}
            {60,[14,11],[2.0 0],[GM 3 52;LIP 57 60],60}
            {61,[14,10],[2.0 0.0 + 0.1],[GM 21 45; LIP 57 75],75}
            {62,[13,12],[2.0 0.0 + 0.2],[GM 2 36;LIP 45 65],65}
            {64,[11,13],[2.1 0],[GM 15 33; LIP 65 78],78}
            {65,[14,9],[1.9 0.0 + 0.2],[GM 22 51; LIP 63 84],84}
            {66,[14,7],[1.8 0],[GM 29 41; GM 65 89; LIP 99 120],120}
            {66,[17,7],[1.8 0.0 + 0.2],[GM 2 6; GM 35 59; LIP 70 90; LIP 113 150],150}
            {67,[15,8],[1.9 0 + 0.1],[GM 34 62; LIP 74 84],84}
            {68,[15,9],[1.9 0],[GM 47 58; LIP 75 83],83}
            {69,[12,9],[1.9 0 + 0.15],[GM 45 73; LIP 86 97],97}
            {70,[16,7],[1.9 0],[GM 41 69; LIP 77 101; VIP 138 150],150}
            {71,[16,6],[1.9 0 + 0.1],[GM 7 19; GM 39 67; LIP 76 80],80}
            {71,[13,7],[1.9 0],[GM 20 38; GM 62 89; LIP 96 110],110}
            }';
        
        MRI_path = 'Z:\Data\MOOG\Polo\Mapping\MRI\PoloOutput\forDrawMapping\';
%         MRI_offset = {[-68 79]*0.9-4  ,[-240 500]*1.1+45 , [0 -2.5]};   % [x1 x2],[y1 y2], [dx/dxSelect slope, dy/dxSelect slope]
        
        MRI_offset = {[-69.2, 63.1], [-175.525, 591.525],[0, -4.33333]};  % Manual 20180604
        
        % Update version with ratio fixed {[x_center, y_center, x_range], [slopes]}. HH20180606
        MRI_offset_new = {[-3.0500, 208, 132.3000], [0, -4.33333]};
        %}
        
       
    case 'Polo_L'
        
        %  Polo_left
        % %{
        
        % Header
        toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
        toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
        gridRange = [0 20; 3 25];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
        
        monkey = 5;
        hemisphere = 1;
        AP0 = 7 % 4;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {62,[6,18],[2.0 0.0],[GM 30 43; LIP 64 72; VIP 72 84; MST 136 150]}
            {63,[6,19],[2.0 0],[GM 4 19; GM 22 39; LIP 56 63],63}
            {64,[6,15],[2.0 0],[GM 15 27; GM 48 69; VIP 81 84],84}
            {65,[6,16],[1.9 0],[GM 52 75; LIP 82 97],97}
            {66,[6,20],[2.1 0],[GM 2 19;LIP 37 60],60}
            {67,[6,21],[2.1 + 0.1 0],[LIP 20 43],43}
            {68,[6,17],[2.0 + 0.1 0],[GM 0 13; GM 25 48; LIP 55 64; LIP 64 75],75}
            {69,[5,19],[2.0  0],[GM 0 41; LIP 53 59],59}
            {70,[5,18],[2.0 0],[GM 0 3; GM 28 53; LIP 62 72],72}
            {71,[5,17],[2.0 0.1],[GM 0 11; GM 27 48; LIP 57 70],70}
            {72,[6,19],[2.1 0],[LIP 46 65]}
            {73,[6,18],[2.1 0],[GM 0 12; GM 26 38; LIP 47 59],59}
            {74,[7,18],[2.0 0],[GM 6 49; LIP 59 73],73}
            {75,[7,17],[2.0 0],[GM 7 54; LIP 62 75],75}
            {76,[7,16],[2.0 0],[GM 2 27; GM 40 61; LIP 73 80; VIP 80 86],86}
            {77,[7,19],[2.1 0],[GM 0 11; LIP 37 53],53}
            {78,[8,17],[2.0 0],[GM 9 50; LIP 66 78],78}
            {79,[8,18],[2.0 0],[GM 0 44; LIP 56 57],57}
            {80,[8,19],[2.1 0],[GM 0 24; LIP 48 60]}
            {81,[7,17],[2.1 0],[GM 0 14; GM 26 46; LIP 59 66]}
            {82,[6,18],[2.1 0],[GM 0 13; GM 27 41; LIP 56 71]}
            {83,[6,12],[2.0 0],[VIP 63 71],71}
            {83,[6,10],[2.0 - 0.1 0],[GM 17 73; VIP 102 123]}
            {84,[6,14],[2.0 0],[GM 5 20; GM 52 81; VIP 92 113]}
            {85,[7,18],[2.1-0.1 0],[GM 8 48; LIP 58 73],73}
            {86,[8,18],[2.1 0],[GM 0 35; LIP 46 65; MST 99 100],100}
            {87,[4,18],[2.0 0],[GM 0 10; GM 31 51; LIP 62 87]}
            {88,[4,17],[2.0 0],[GM 35 58; LIP 71 94 ]}
            {89,[4,19],[2.1-0.1 0],[GM 20 45; LIP 61 87]}
            {90,[3,17],[2.0 0],[GM 35 70; LIP 83 97; VIP 97 110]}
            {91,[3,18],[2.0 0],[GM 0 10; GM 37 60; LIP 73 92; VIP 92 98]}
            {91,[3,19],[2.0 0],[GM 30 45; LIP 63 80]}
            {92,[10,17],[2.0 0],[GM 3 45; LIP 52 69],69}
            {93,[12,15],[2.1 0],[GM 3 27; LIP 37 61]}
            {94,[12,13],[2.0 0],[GM 21 49; LIP 56 72],72}
            {95,[14,12],[1.9 -0.1],[GM 10 25; GM 35 56; LIP 67 89]}
            {96,[16,10],[1.9 -0.1],[GM 9 31; GM 42 70; LIP 100 120; LIP 133 148]}
            {97,[10,14],[2.0 0],[GM 3 10; GM 25 61; LIP 71 81],81}
            {98,[8,15],[2.0 -0.05],[GM 3 22; GM 42 65; LIP 78 89],89}
            {99,[9,15],[2.0 0],[GM 5 62; LIP 73 86],86}
            {100,[11,14],[2.0 0],[GM 3 56; LIP 66 81],81}
            {101,[13,13],[2.0 0],[GM 11 42; LIP 52 54],54}
            {102,[8,11],[2.0 0],[GM 2 17; GM 35 42; VIP 65 89; VIP 95 110]}
            {103,[8,10],[1.9 0],[GM 11 26; GM 36 52; VIP 73 92; VIP 101 119]}
            {104,[8,9],[1.9 0],[GM 14 30; GM 40 59; VIP 82 96 ; VIP 107 115]}
            {105,[8,12],[1.9 0],[VIP 67 91; VIP 103 126]}
            {106,[9,11],[1.8 0-0.1],[GM 30 58; GM 80 102; VIP 117 132],132}
            {107,[9,10],[1.8 0],[GM 16 52; GM 74 91; VIP 107 125]}
            {108,[10,10],[1.8 0],[GM 15 58; GM 68 83; VIP 98 125]}
            {109,[10,9],[1.8 0],[GM 4 25; GM 34 51; GM 72 85; VIP 101 117],117}
            {111,[10,15],[2.1 0+0.1],[GM 0 35; LIP 44 59],59}
            }';
        
        MRI_path = 'Z:\Data\MOOG\Polo\Mapping\MRI\PoloOutput\forDrawMapping\';
%         MRI_offset = {[-68 79] - 5,[-240 500]+15, [0 -1.5]};   % [x1 x2],[y1 y2], [dx/dxSelect slope, dy/dxSelect slope]
        
        MRI_offset = {[-61.4862, 66.4862], [-220, 520],[0, -1.5]}; % Manual adjust HH20180604
        
        % Update version with ratio fixed {[x_center, y_center, x_range], [slopes]}. HH20180606
        MRI_offset_new = {[ 2.5000,   150,  127.9724 ], [0, -1.5]};

        %}
    case 'Polo_R'
        
        %  Polo_right
        % %{
        
        % Header
        toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
        toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
        gridRange = [0 20; 0 25];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]
        monkey = 5;
        hemisphere = 2;  % L = 1, R = 2
        AP0 =  Polo_right_AP0 % 6;
        
        data = {
            %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
            % When you are not sure about one area, use "AreaType-100" instead
            {1,[5,15],[2.0 0.0],[GM 0 20; LIP 44 69; MST 147 150]}
            {2,[5,17],[2.0 0.0],[GM 0 17; LIP 32 64],64}
            {3,[5,13],[1.9 0.0],[GM 21 37; LIP 71 96]}
            {3,[5,11],[1.9 0.0],[GM 31 58; LIP 74 97; VIP 97 103]}
            {4,[6,7],[1.7 0.0],[GM 0 21; VIP 78 100; VIP 111 129],129}
            {5,[5,16],[1.7 0.0],[GM 21 51; LIP 76 82],82}
            {6,[5,14],[1.8 0.0],[GM 0 15; GM 30 43; LIP 79 98]}
            {7,[5,12],[1.8 0.0],[GM 38 50; LIP 79 107],107}
            {8,[9,10],[1.8 0.0],[GM 43 73; LIP 91 111],111}
            {9,[9,9],[1.9 0.0],[GM 0 12; GM 35 59; LIP 75 94],94}
            {10,[9,11],[1.9 0.0],[GM 28 56; LIP  75 90],90}
            {11,[9,12],[1.9 0.0],[GM 24 49; LIP 68 69], 69}
            {12,[9,13],[1.9 0.0],[GM 17 35; LIP 53 62; MST 101 105 ], 105}
            {13,[7,11],[1.9 0.0],[GM 0 9; GM 25 47; LIP 61 86; MST 127 149]}
            {14,[7,12],[1.9 0.0],[GM 0 11; GM 26 47; LIP 58 76],76}
            {15,[7,13],[1.9 0.0],[GM 0 11; GM 20 41; LIP 50 68],68}
            {16,[2,10],[1.9 0.0],[VIP 65 86],86}
            {17,[2,11],[1.9 0.0],[GM 49 59 ; VIP 59 77; VIP 90 103]}
            {18,[2,12],[1.9 0.0],[GM 47 75; VIP 88 102]}
            {19,[8,11],[1.9 0.0],[GM 25 53; LIP 61 84; MST 118 138]}
            {20,[3,10],[1.9 0.0],[VIP 57 73],73}
            {21,[3,9],[1.9 0.0],[VIP 68 92],92}
            {22,[3,11],[1.9 0.0], [GM 47 73; LIP 92 99; VIP 99 110]}
            {23,[2,9],[1.9 0.0],[VIP 77 105]}
            {25,[4,9],[1.9 0.0],[GM 68 80; VIP 80 85],85}
            {26,[4,10],[1.9 0],[GM 49 81; VIP 95 107],107}
            {27,[6,11],[1.9 0],[GM 36 66; LIP 79 100],100}
            {28,[6,12],[1.9 0],[GM 29 55; LIP 71 90; MST 142 150]}
            {29,[6,13],[2.1 0],[GM 13 36; LIP 44 60],60}
            {30,[8,10],[1.9 0],[GM 9 25; GM 34 54; LIP 63 73],73}
            {31,[8,9],[1.9 0],[GM 11 22; GM 39 63; LIP 77 84],84}
            {32,[7,10],[1.9 0],[GM 2 17; GM 39 64; LIP 74   98],98}
            {33,[8,12],[2.0 0],[GM 4 34; LIP 45 65; MST 108 110],110}
            {34,[4,12],[2.0 0],[GM 24 43; LIP 67 85],85}
            {35,[4,13],[2.0 0],[GM 19 38; LIP 60 82],82}
            {36,[4,14],[2.0 0],[GM 0 38; LIP 55 78],78}
            {37,[3,13],[2.0 0],[GM 22 43; LIP 72 86]}
            {39,[8,10],[2.0 0],[GM 10 47; LIP 57 81],81}
            {40,[7,11],[2.0 0],[GM 24 45; LIP 58 70; VIP 70 80; MST 125 129],129}
            {38,[3,14],[2.0 0],[GM 18 46; LIP 66 84]}
            {41,[7,12],[2.0 0],[GM 14 37; LIP 48 66], 66}
            {42,[6,14],[2.0 0],[GM 5 30; LIP 45 58],58}
            {43,[6,15],[2.0 0],[GM 6 27; LIP 43 70; MST 144 150]}
            {44,[6,16],[2.0 0],[LIP 30 67; MST 106 135]}
            {45,[3,15],[2.0 0],[GM 11 35; LIP 61 63],63}
            {46,[3,16],[2.0 0],[GM 13 35; LIP 52 56],56}
            {47,[3,18],[2.0 0],[GM 10 38; LIP 53 57],57}
            {48,[3,17],[2.1 0],[GM 0 25; LIP 42 56],56}
            {49,[4,15],[2.0 0],[GM 2 35; LIP 55 71],71}
            {50,[2,14],[2.0 0],[GM 0 2; GM 30 55; LIP 69 73; VIP 73 86];}
            {51,[2,15],[2 0],[GM 27 53; LIP 68 73], 73}
            {52,[2,16],[2.1 0],[GM 0 43; LIP 54 65],65}
            {53,[10,10],[1.9 0],[GM 0 46; VIP 60 68],68}
            {54,[10,12],[2.0 0],[GM 6 15; LIP 33 59; MST 85 90],90}
            {55,[10,13],[2.0 0],[LIP 30 47; MST 92 93],93}
            {56,[8,13],[2.0 0],[GM 0 22; LIP 41 60],60}
            {57,[7,14],[2.0 0],[GM 5 28; LIP 36 58; MST 120 144]}
            {58,[6,13],[2.0 0],[GM 19 33; LIP 46 48],48}
            {59,[6,12],[1.9 0],[GM 40 50; LIP 71 86; MST 148 150]}
            {60,[5,13],[2.0 0],[GM 12 27; LIP 57 62],62}
            {61,[5,14],[2.0 0],[LIP 58 67],67}
            }';
        
        MRI_path = 'Z:\Data\MOOG\Polo\Mapping\MRI\PoloOutput\forDrawMapping\';
%          MRI_offset = {[-68 79] - 10,[-240 500], [0 -1.5]};   % [x1 x2],[y1 y2], [dx/dxSelect slope, dy/dxSelect slope]
        
        MRI_offset = {[-71.9026, 63.9026], [-244.948, 524.948],[0, -1.66667]};  % Manual adjust 20180604
        
        % Update version with ratio fixed {[x_center, y_center, x_range], [slopes]}. HH20180606
        MRI_offset_new = {[ -4.0000,   140 ,  135.8052 ], [0, -1.66667]};

        %}
end

%  Hetao_left
%{

toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
gridRange = [5 25; 1 22];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]

monkey = 2;
hemisphere = 1;
data = {
    %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] }
    % When you are not sure about one area, use "AreaType-100" instead
    {83,[11,16],[1.9 0.0],[GM 4 35]}
    {83,[11,18],[2.1 0.0],[GM 9 25; GM 68 81; VIP 94 107]}
    {84,[11,17],[2.1 0.0],[GM 8 30; GM 72 89; VIP 93 101]}
    {84,[11,19],[2.1 0.0],[GM 45 70; VIP 77 90]}
    {85,[11,20],[2.1 0.0],[GM 40 66; VIP 73 95]}
    {86,[11,21],[2.15 0],[GM 37 56; AUD 106 120]}
    {86,[11,22],[2.15 0],[GM 18 52; LIP 62 82; AUD 105 140]}
    {87,[12,17],[1.9 0.0],[GM 24 38; GM 67 90; VIP 105 123]}
    {88,[12,18],[2.05 0],[GM 0 40; GM 49 63; VIP 87 99]}
    {89,[12,19],[2.05 0],[GM 34 59; VIP 67 90]}
    {90,[12,16],[1.9 0.0],[GM 8 27; VIP 85 120]}
    {90,[12,15],[1.9 0.0],[GM 12 30; VIP 99 119]}
    };
%}

%  Hetao_right
%{

toPlotTypes3D = [ MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?
toPlotTypes_zview = [MST VIP MT LIP AUD];    % Which area types do we want to plot ?
gridRange = [1 25; 1 22];  % Grid plotting ragne = [xLow xHigh; yLow yHigh]

monkey = 2;
hemisphere = 2;
data = {
    %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] }
    % When you are not sure about one area, use "AreaType-100" instead
    {1,[20,15],[1.8 0],[GM,43 82; MT,122 128; MT,131 150]};
    {2,[20,13],[1.8 0],[GM,40 89; VM,105 114; GM,124 150]};
    {3,[22,15],[2.0 0],[GM,40 60; MST, 82 130]};
    {4,[18,15],[2.0 0],[GM,40 90; MST, 103 122; MT, 135 140]};
    {5,[18,13],[1.8 0.2],[GM,38 60; MST, 80 150]};
    {6,[18,11],[1.7 0.1],[GM,60 100; MST, 135 150]};                % u-stim
    {7,[18,10],[1.7 0.3],[GM,50 86; MST, 117 150]};
    {8,[18,9],[1.7 0.5],[GM 43 54; LIP 55 76; MST 110 120]};
    {9,[18,8],[1.7 0.1],[GM, 38 63; LIP 91 103]};
    {9,[18,17],[2.0 0.2],[GM 39 50; MST 95 137; MT 153 189]};
    {10,[18,7],[1.7 0.2],[GM 41 57; VM 57 69; LIP 91 100]};
    {10,[18,6],[1.7 0],[GM 69 85; LIP 108 125]}
    {11,[18,5],[1.6 0.3],[GM 60 83; LIP 90 109]}
    {11,[18,4],[1.6,0.5],[GM 40 75; LIP 82 108]}
    {12,[16,15],[1.8 0.5],[MST 123 144]}
    {[12,15],[16,13],[1.8 0],[GM 60 76; LIP 77 100; MST 130 177]}    % u-stim
    {13,[16,10],[1.9 0.2],[GM 43 67; LIP 68 80]}
    {13,[16,8],[1.9 0.2],[GM 0 74; LIP 81 106]}
    {14,[16,17],[2.1 0.4],[MST 42 72; MST 85 115]}
    {16,[14,17],[2.0 0.5],[GM 12 40; VM 57 77; MST 91 130]}     % u-stim
    {17,[14,15],[2.1 0.3],[GM 43 74; MST 95 124]}               % u-stim
    {17,[14,19],[2.1 0.5],[VM 30 56; VM 75 97; VM 110 130]}
    {18,[12,15],[2.1 0.3],[GM 20 59; GM 76 84; VM 103 130]}    % u-stim
    {18,[12,17],[2.1 0.5],[GM 27 50; GM 63 79; VM 93 137]}
    {19,[12,19],[2.1 0.5],[GM 44 77; GM 96 118]}  % Sound!!
    {19,[20,17],[2.1 0.5],[GM 28 61; MST 62 99; MT 117 136]}  % u-stim
    {20,[16,19],[2.3 0.2],[GM 30 56; MST 85 109; MT 116 139]}   % u-stim at MT
    {21,[18,19],[2.3 0.3],[GM 23 36; VM 37 45]}
    {22,[15,19],[2.3 0.5],[GM 48 74; VM 96 114; MT 123 137]}
    {22,[20,19],[2.3 0.5],[GM 21 32; VM 33 64; MT 82 130]}
    {23,[22,17],[2.3 0.5],[MST 42 62; MT 63 103]}   % u-stim at MT
    {24,[22,13],[2.3 0.2],[VM 63 100]}
    {24,[22,19],[2.3 0.5],[VM 33 45]}
    {25,[20,11],[2.1 0.3],[GM 25 51; MST 63 103]}   % u-stim, sigma, miu
    {26,[20,9],[2.1 0.3],[GM 27 37; MST 42 84]}    % u-stim
    {27,[22,11],[2.1 0.3],[MST 56 87; VM 148 150]}      % u-stim
    {28,[20,12],[2.1 0.3],[GM 27 54; MST 79 117]}
    {[29,2],[20,13],[2.1 0.3],[MST 50 82; VM 145 150]}
    {30,[20,10],[2.0 0.2],[VM 53 82; MST 100 130]}  % LIP?
    {30,[19,11],[2.0 0.3],[MST 52 67; MST 103 126]}
    {[1,31],[20,15],[2.0 0.2],[VM 25 75; MST 80 108; MT 123 135]}  % Oscillation?   % u-stim
    {32,[20,7],[1.8 0.1],[VM 0 45; LIP 59 83]}
    {33,[20,8],[2.1 0.1],[LIP 0 42]}
    {34,[18,3],[2.1 0.1],[MIP 35 72; VIP 92 106]}
    {35,[18,2],[2.1 0.0],[GM 27 35; MIP 56 87; VIP 88 109]}
    {36,[18,1],[2.0 0.1],[GM 29 36; MIP 74 100; VIP 101 127]}
    {37,[14,5],[2.0 0.1],[GM 33 45; GM 51 57; GM 90 102; VIP 107 118]}
    {38,[14,7],[1.9 0.0],[GM 89 115; LIP 125 138]}
    {38,[14,6],[1.9 0.2],[GM 78 102; LIP 118 125]}
    {39,[14,2],[1.9 0.1],[GM 34 49; GM 78 95]}
    {39,[14,3],[1.9 0.1],[GM 34 40; GM 72 87]}
    {40,[14,4],[1.9 0.1],[GM 45 59; LIP 101 110; VIP 115 130]}
    {41,[16,4],[1.7 0.2],[GM 35 43; GM 52 67; VIP 94 114; VIP 137 149]}
    {42,[16,6],[1.7 0.3],[GM 27 36; VM 62 80; LIP 98 122]}
    {42,[16,2],[1.7 0.4],[GM 28 39; VIP 83 118]}
    {43,[12,8],[2.0 0.0],[VIP 94 115]}
    {45,[12,6],[2.0 0.0],[VIP 84 112]}
    {45,[12,20],[2.3 0.2],[GM 62 94; AUD 110 130]}
    {46,[12,11],[1.9 0.0],[GM 46 76; LIP 95 113]}
    {46,[12,4],[1.8 0.0],[GM 45 63; GM 87 102]}
    {47,[10,10],[1.9 0.0],[VIP 65 89; VIP 95 106]}
    {48,[10,12],[1.9 0.0],[GM 45 70; VIP 81 102]}
    {48,[10,8],[1.9 0.0],[GM 0 1]}
    {49,[10,14],[1.9 0.0],[GM 9 32; GM 40 61; LIP 78 91; VIP 93 98]}
    {50,[8,12],[1.9 0.0],[VIP 64 90; AUD 127 140]}
    {50,[12,12],[1.9 0.0],[GM 27 47; LIP 67 89]}
    {51,[8,10],[1.8 0.0],[GM 30 40]}
    {51,[8,8],[1.8 0.0],[GM 30 60]}
    {51,[8,14],[1.9 0.0],[GM 14 25; GM 65 89; LIP 98 110; GM 133 145]}
    {52,[6,12],[1.9 0.0],[GM 7 60]}
    {52,[6,14],[1.9 0.0],[GM 12 30; VIP 88 105; GM 138 150]}
    {52,[6,16],[1.9 0.3],[GM 56 77; GM 114 134]}
    {53,[18,12],[2.1 0.0],[VM 7 56; MST 110 128]}
    {54,[12,8],[1.9 0.1],[GM 66 74; VIP 74 117]}
    {55,[12,10],[1.9 0.0],[GM 4 12; GM 54 80; VIP 93 109]}
    {56,[11,11],[1.9 0.0],[GM 60 87; VIP 93 106]}
    {57,[12,9],[1.7 0.0],[GM 20 40; VIP 96 125]}
    {58,[12,7],[1.8 0.0],[GM 10 23; VIP 100 118]}
    {59,[11,10],[1.8 0.0],[GM 44 90; VIP 111 146]}
    {60,[11,9],[1.8 0.0],[GM 3 15; VIP 87 113]}
    {61,[12,11],[1.9 0.0],[GM 35 65; LIP 79 100]}
    {63,[13,9],[1.9 0.0],[GM 12 29; LIP 86 115]}
    {65,[10,11],[1.75 0.0],[GM 9 22; LIP 92 102; VIP 104 124]}
    {66,[11,12],[1.9 0.0],[GM 12 25; GM 67 93; VIP 106 109]}
    {66,[11,13],[1.9 0.0],[GM 12 23; GM 53 78; VIP 96 111]}
    {67,[10,9],[1.75 0.0],[GM 0 1]}
    {68,[9,11],[1.75 0.0],[GM 10 17]}
    {69,[9,12],[1.75 0.0],[GM 9 14; VIP 89 113]}
    {70,[9,14],[1.85 0.0],[GM 9 26; LIP 70 91]}
    {71,[20,12],[2.1 0.3],[GM 38 65; MT 116 133]}
    {73,[20,11],[2.0 0.3],[GM 22 60; MST 91 96]}
    {74,[19,13],[2.0 0.0],[GM 8 23; MST 86 103]}
    {75,[19,12],[2.0 0.0],[MST 62 95; MST 124 150]}
    {76,[18,14],[2.1 0.2],[MST 95 122; MST 134 139]}
    {77,[18,16],[2.1 0.2],[GM 40 57; MST 102 120]}
    {78,[12,9],[1.7 0.0],[GM 31 47; GM 100 108; VIP 112 131]}
    {78,[12,10],[1.7 0.2],[GM 3 16; GM 70 80; VIP 82 95]}
    {79,[9,13],[2.0 0.0],[GM 10 74; LIP 84 99]}
    {80,[13,8],[1.9 0.0],[GM 70 80; VIP 82 97]}
    {80,[13,7],[1.9 0.0],[GM 26 68; VIP 95 119]}
    {81,[10,13],[1.8 0.0],[GM 7 32; GM 49 78; LIP 90 107]}
    {81,[11,12],[2.1 0.0],[GM 0 20; GM 42 70; VIP 80 100]}
    {82,[13,6],[1.9 0.0],[GM 11 24; VIP 120 141]}
    {82,[13,5],[1.9 0.1],[GM 0 17; GM 54 83]}
    {82,[12,5],[1.9 0.0],[GM 4 18]}
    
    };
%}

%% ========== 3-D Visualization ========== %%

% toPlotTypes3D = [GM MST VIP MT LIP AUD];    % Which area types do we want to plot in 3D plot?

% toPlotTypes3D = VIP ;    % Which area types do we want to plot in 3D plot?


for channel = 1:length(data)
    gridLoc = data{channel}{2};
    GMData = data{channel}{4};
    if isempty(GMData); continue; end;
    
    GuideTubeAndOffset = data{channel}{3};
    offSet = round((GuideTubeAndOffset(2) + GuideTubeAndOffset(1) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
    
    for GMType = -size(GMTypes,1):-1  % For each area type
        GMThisType = find(GMData(:,1) == GMType);   % Read out ranges for this type
        if isempty(GMThisType) || isempty(intersect(toPlotTypes3D,GMType)); continue; end       % If absent in this channel or absent in areas we want to plot, next type
        
        GMPos = [];    % Clear cache
        for i = 1:length(GMThisType)    % For each appearance
            % Attach each appearance to position cache
            GMPos = [GMPos (maxZ - offSet - GMData(GMThisType(i),3)):(maxZ - offSet - GMData(GMThisType(i),2))];
        end
        
        % Add color to positions in the cache
        try
            if hemisphere == 2
                V(gridLoc(1), maxY - gridLoc(2), GMPos, :) = repmat(GMTypes{-GMType,2},length(GMPos),1);
            else
                V(gridLoc(1), gridLoc(2), GMPos, :) = repmat(GMTypes{-GMType,2},length(GMPos),1);
            end
        catch
            fprintf('Warning: Out of Range (Session %g, Channel [%g,%g], GMType %s)\n', channel, data{channel}{2}(1), data{channel}{2}(2), GMTypes{-GMType,1});
            keyboard;
        end
    end
    
    % Add start and end markers
    if hemisphere == 2
        V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - 2):(maxZ - offSet - 1),:) = repmat([0 0 0],2,1);
        if length(data{channel}) <5
            V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - movDis - 1):(maxZ - offSet - movDis),:) = repmat([0 0 0],2,1);
        else
            V(gridLoc(1), (maxY - gridLoc(2)),(maxZ - offSet - data{channel}{5} - 1):(maxZ - offSet - data{channel}{5}),:) = repmat([0 0 0],2,1);
        end
    else
        V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - 2):(maxZ - offSet - 1),:) = repmat([0 0 0],2,1);
        if length(data{channel})<5
            V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - movDis - 1):(maxZ - offSet - movDis),:) = repmat([0 0 0],2,1);
        else
            V(gridLoc(1), (gridLoc(2)),(maxZ - offSet - data{channel}{5} - 1):(maxZ - offSet - data{channel}{5}),:) = repmat([0 0 0],2,1);
        end
    end
    
end

% Render the mapping results using "vol3d" function
% close all;
set(figure(801),'Position',[10 100 600 600]);  clf;
set(0, 'DefaultAxesXTickMode', 'auto', 'DefaultAxesYTickMode', 'auto', 'DefaultAxesZTickMode', 'auto');
vol3d('cdata',V);
view(-150,30);

axis tight;
daspect([1,1,10]);
alphamap('rampup');
alphamap(.5 .* alphamap);  % Transparency

grid on;  grid minor;
set(gca,'GridLineStyle','-');

set(gca,'xtick',-0.5:5:maxY-0.5);
if hemisphere == 2
    %     set(gca,'xticklabel','25|20|15|10|5|0');
    set(gca,'xticklabel',25:-5:0);
else
    %     set(gca,'xticklabel','0|5|10|15|20|25');
    set(gca,'xticklabel',0:5:25);
end
xlim([-1 maxY]);

set(gca,'ytick',-0.5:5:maxX-0.5);
% set(gca,'yticklabel','0|5|10|15|20|25|30');
set(gca,'yticklabel',0:5:30);
ylim([-1 maxX]);

set(gca,'ztick',0:50:maxZ);
% set(gca,'zticklabel','-250|-200|-150|-100|-50|0');
set(gca,'zticklabel',-250:50:0);

ylabel('Posterior (X)')
xlabel('Lateral (Y)')

for i = 1:length(toPlotTypes3D)
    text(0, 0, maxZ-i*40, GMTypes{-toPlotTypes3D(i),1},'color',GMTypes{-toPlotTypes3D(i),2},'FontSize',20);
end;

set(gcf,'color','w');
% set(findall(gcf,'fontsize',10),'fontsize',20);
SetFigure(20);

% k=1;
% for i = 1:150
%     view(-170+i*.5,20+i/10);
%     drawnow;
%     mov(k) = getframe(gcf);
%     k=k+1;
% end
%
% for i = 150:-1:1
%     view(-170+i*.5,20+i/10);
%     drawnow;
%     mov(k) = getframe(gcf);
%     k=k+1;
% end
%
% movie2avi(mov,'Test.avi');

%% ============ 2-D Visualization (Grid view from the top) =============== %

radius = 0.42;  % Radius of each hole (interval = 1)

% Plot grid outline
set(figure(802),'Position',[10 100 600 600]); clf
hold on; axis equal ij;  

global h_grid_axis; h_grid_axis = gca;

x1 = gridRange(1,1); x2 = gridRange(1,2);
y1 = gridRange(2,1); y2 = gridRange(2,2);
xlim([y1-2,y2+2]);

% Quick select monkey_hemi. HH20160122
h_monkey_hemi = uicontrol('style','popupmenu','unit','norm','position',[0.01 0.94 0.131 0.035],...
        'string','Polo_L|Polo_R|Messi_L|Messi_R','callback',@change_monkey_hemi);
set(h_monkey_hemi,'value',monkey_hemi);

% Frame
interval = 5;
xLoc = intersect(x1:x2,0:interval:100);
yLoc = intersect(y1:y2,0:interval:100);
xLines = line(repmat([y1-1;y2+1],1,length(xLoc)),repmat(xLoc,2,1));  % xLines
yLines = line(repmat(yLoc+0.25,2,1), repmat([x1-1;x2+1],1,length(yLoc)));  % yLines
set(xLines,'LineWidth',5,'Color','g');
set(yLines,'LineWidth',5,'Color','g');
set(gca,'xtick', yLoc);
set(gca,'ytick', xLoc);

% Parallel drawing (to speed up)
xOffsets = repmat(x1:x2, y2-y1+1,1);
xOffsets = xOffsets(:)';
yOffsets = repmat([(y1:y2)+0.5*~mod(x1,2) (y1:y2)+0.5*mod(x1,2)], 1, fix((x2-x1 + 1)/2));
if mod(x2-x1+1,2) == 1
    yOffsets = [yOffsets y1:y2];
end
t = linspace(0,2*pi,100)';
xGrid = radius * repmat(sin(t),1,length(xOffsets)) + repmat(xOffsets,length(t),1);
yGrid = radius * repmat(cos(t),1,length(yOffsets)) + repmat(yOffsets,length(t),1);
set(fill(yGrid,xGrid,[1 1 1]),'LineWidth',1.5,'ButtonDownFcn',@SelectChannel);    % White color

% Plot mapping result
for channel = 1:length(data)
    xCenter = data{channel}{2}(1);
    yCenter =  data{channel}{2}(2) + 0.5 * ~mod(data{channel}{2}(1),2);
    
    % Paint
    haveTypes = intersect(toPlotTypes_zview,data{channel}{4}(:,1));    % Find which areas this channel has.
    if isempty(haveTypes)    % If there is no types of interest, paint using the color of GM
        t = linspace(0,2*pi,100)';
        xMap = radius * sin(t) + xCenter ;
        yMap = radius * cos(t) + yCenter;
        set(fill(yMap,xMap,GMTypes{-GM,2}),'ButtonDownFcn',@SelectChannel);
        c = 'k';
    else                     % Else, we paint colors for each type
        for iType = 1:length(haveTypes)
            t = linspace(2*pi/length(haveTypes)*(iType-1), 2*pi/length(haveTypes)*iType, 100)';
            xMap = 0.99 * radius * [0; sin(t)] + xCenter;
            yMap = 0.99 * radius * [0; cos(t)] + yCenter;
            set(fill(yMap,xMap,GMTypes{-haveTypes(iType),2}),'LineStyle','none','ButtonDownFcn',@SelectChannel);
        end
        
        if [0.21 0.72 0.07]*(GMTypes{-haveTypes(1),2})' > 0.5
            c = 'k';
        else
            c = 'y';
        end
    end
    
    %------ Denotations -----%
    % Session No.
    text(yCenter, xCenter-0.1, num2str(data{channel}{1}),'FontSize',7,'color',c,'HorizontalAlignment','center','ButtonDownFcn',@SelectChannel);
    
end

% Legend
for i = 1:length(toPlotTypes_zview)
    text(y1+(i-1)*4, x1-2, GMTypes{-toPlotTypes_zview(i),1},'color',GMTypes{-toPlotTypes_zview(i),2},'FontSize',15);
end;

set(gcf,'color','w')
set(gca,'fontsize',15)

if hemisphere == 1
    set(gca,'xdir','rev');
end

%% Load xls file
global num txt raw xls;
if isempty(num)
    XlsData = ReadXls('Z:\Labtools\HH_Tools\DataHub\DataHub.xlsm',2,3);
    num = XlsData.num;
    xls = XlsData.header;
    txt = XlsData.txt;
    
    disp('Xls Loaded');
end

%% Call back for 2-D Grid View
function SelectChannel(source,event)

global maxX maxY maxZ movDis GMTypes data Polo_right_AP0
global MRI_path MRI_offset AP0; global MRI_offset_new;
global monkey hemisphere ;
global gridRange linWid start_end_markers overlapping;
global handles;
global num txt raw xls;
global overlapTransparent; 
global h_grid_axis;

hemisphere_text = {'Left','Right'};

% Which channel do we select?

pos = get(h_grid_axis,'CurrentPoint');
pos = pos(1,1:2);
xSelect = round(pos(2));
ySelect = round(pos(1)- 0.5 * ~mod(xSelect,2));

figure(802);

delete(handles);

% handles(100) = rectangle('Position',[gridRange(2,2)+1.4,xSelect-0.1,0.2,0.2],'FaceColor','k');
handles(1) = plot(xlim,repmat(xSelect - 0.5 + overlapping(1,1),1,2),'r--','LineWid',1);
handles(2) = plot(xlim,repmat(xSelect + 0.5 + overlapping(1,2),1,2),'r--','LineWid',1);
handles(3) = plot(repmat(ySelect + 0.25,1,2),[xSelect - 0.5 + overlapping(1,1) xSelect + 0.5 + overlapping(1,2)],'r','LineWid',1);

% Show corresponding sessions
sessions_match = [];
for channel = 1:length(data)
    if data{channel}{2}(1) == xSelect && data{channel}{2}(2) == ySelect
        sessions_match = [sessions_match data{channel}{1}];
    end
end

title(sprintf('Monkey %g, session(s) for %s [%g,%g]: %s',monkey,hemisphere_text{hemisphere},xSelect,ySelect,num2str(sessions_match)));


figurePosition = [ 700  100  459  600 ];
set(figure(803),'Position',figurePosition,'color','w'); clf


% 2-D Visualization (Coronal)
h_coronal = axes('Position',[0.15 0.1 0.8 0.8]);

axis ij; hold on;

% Overlapping MRI data
try
    %     if hemisphere == 1 && monkey == 5
    fileNo = (xSelect - AP0) + Polo_right_AP0 ;  % Because I use Polo_right MRI as standard images
    %     elseif monkey == 10
    %         fileNo = xSelect + Polo_right_AP0 - AP0;
    %     else
    %         fileNo = xSelect + Polo_right_AP0 - AP0;
    %     end
    
    MRI = imread([MRI_path num2str(fileNo) '.bmp']);
    
    % h_MRI = image(MRI_offset{1}+xSelect*MRI_offset{3}(1), MRI_offset{2} + xSelect*MRI_offset{3}(2),MRI);
    
    % MRI_offset_new = {[x_center, y_center, x_range], [slopes]}
    ratioyx = size(MRI,1)/size(MRI,2) * 10 * 0.8; % Auto keep ratio !!!   HH20180606
    x_lims = [MRI_offset_new{1}(1) - MRI_offset_new{1}(3)/2, MRI_offset_new{1}(1) + MRI_offset_new{1}(3)/2];
    y_lims = [MRI_offset_new{1}(2) - ratioyx * MRI_offset_new{1}(3)/2, MRI_offset_new{1}(2) + ratioyx * MRI_offset_new{1}(3)/2];
    
    h_MRI = image(x_lims + xSelect * MRI_offset_new{2}(1), y_lims + xSelect * MRI_offset_new{2}(2), MRI);
    
    set(h_MRI,'AlphaData',1);
    
    % 20180604
    uicontrol('style','pushbutton','unit','norm','pos', [0.0630    0.9510    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 1});
    uicontrol('style','pushbutton','unit','norm','pos', [0.0630    0.9210    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 2});
    uicontrol('style','pushbutton','unit','norm','pos', [0.0330    0.9350    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 3});
    uicontrol('style','pushbutton','unit','norm','pos', [0.0930    0.9350    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 4});

    uicontrol('style','pushbutton','unit','norm','pos', [0.1630    0.9510    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 5});
    uicontrol('style','pushbutton','unit','norm','pos', [0.1630    0.9210    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 6});
                 
    uicontrol('style','pushbutton','unit','norm','pos', [0.2430    0.9510    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 7});
    uicontrol('style','pushbutton','unit','norm','pos', [0.2430    0.9210    0.0300    0.0210],...
                     'callback',{@ manual_adjust_mri, 8});

catch
    disp('No MRI data found...');
end

% Frame
xlim([gridRange(2,1) gridRange(2,2)+1]);
ylim([-30 maxZ]);

grid minor;
set(h_coronal,'XMinorGrid','on','XMinorTick','on');
% title(sprintf('Monkey %g, %s [%g, %g]',monkey, hemisphere,xSelect,ySelect));
title(sprintf('Monkey %g, %s[%g], AP \\approx %g',monkey, hemisphere_text{hemisphere},xSelect,(AP0-xSelect)*0.8));

% Keep scale
aspectRatio = (range(ylim) * 100) / (range(xlim) * 800);  % grid interval = 0.8 mm
set(figure(803),'Position',[figurePosition(1:2) figurePosition(4)/aspectRatio figurePosition(4)]);
daspect([1/0.8 10 1]); % This is betther man

for channel = 1:length(data)
    %     if data{channel}{2}(1) == xSelect  % Only plot the line we select
    if data{channel}{2}(1) >= xSelect + overlapping(1,1) && data{channel}{2}(1) <= xSelect + overlapping(1,2)   % Overlapping neighboring slices
        yLoc = data{channel}{2}(2)-0.5;
        GMData = data{channel}{4};
        if isempty(GMData); continue; end;
        
        GuideTubeAndOffset = data{channel}{3};
        offSet = round((GuideTubeAndOffset(2) + GuideTubeAndOffset(1) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
        
        for GMType = -size(GMTypes,1):-1  % For each area type
            GMThisType = find(GMData(:,1) == GMType);   % Read out ranges for this type
            if isempty(GMThisType); continue; end       % If absent, next type
            
            for i = 1:length(GMThisType)  % For each appearance
                zBegin = GMData(GMThisType(i),2) + offSet;
                zEnd = GMData(GMThisType(i),3) + offSet;
                
                if overlapping(1,1)~=0 || overlapping(1,2)~=0  % If we overlap neighboring slices
                    p = patch([yLoc yLoc yLoc+1 yLoc+1],[zEnd,zBegin,zBegin,zEnd],GMTypes{-GMType,2},'EdgeColor','none','FaceAlpha',overlapTransparent);
                else
                    rectangle('Position',[yLoc,zBegin,1,zEnd-zBegin],'LineWidth',linWid,...
                        'EdgeColor',GMTypes{-GMType,2});
                end
            end
            
        end
        
        for GMType = -size(GMTypes,1):-1  % For each area type (that we are NOT SURE!!)
            GMThisType = find(GMData(:,1) == GMType - 100);   % Read out ranges for this type (that we are NOT SURE!!)
            if isempty(GMThisType); continue; end       % If absent, next type
            
            for i = 1:length(GMThisType)  % For each appearance
                zBegin = GMData(GMThisType(i),2) + offSet;
                zEnd = GMData(GMThisType(i),3) + offSet;
                
                % Note we use dotted line here to mark areas that we are not sure
                rectangle('Position',[yLoc,zBegin,1,zEnd-zBegin],'LineWidth',linWid,'EdgeColor',GMTypes{-GMType,2},'LineStyle',':');
            end
            
        end
        
        % Add start and end markers
        if start_end_markers
            if yLoc+0.5 == ySelect
                col_temp = 'c';
                rectangle('Position',[yLoc,offSet,1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                offSet_selected = offSet;
                if length(data{channel})<5
                    rectangle('Position',[yLoc,offSet + movDis,1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                else
                    rectangle('Position',[yLoc,offSet + data{channel}{5},1,1],'LineWidth',linWid,'FaceColor',col_temp,'EdgeColor',col_temp);
                end
            else
                rectangle('Position',[yLoc,offSet,1,1],'LineWidth',linWid,'FaceColor','k');
                if length(data{channel})<5
                    rectangle('Position',[yLoc,offSet + movDis,1,1],'LineWidth',linWid,'FaceColor','k');
                else
                    rectangle('Position',[yLoc,offSet + data{channel}{5},1,1],'LineWidth',linWid,'FaceColor','k');
                end
            end
            
            
        end
        
    end
end

if exist('h_MRI')
    if exist('offSet_selected')
        set(h_MRI,'ButtonDownFcn',{@ShowDepth,offSet_selected});
    else
        set(h_MRI,'ButtonDownFcn','');
    end
else
    set(gca,'color',[0.2 0.2 0.2]);
end

if start_end_markers
    % Channel indicator
    plot([ySelect ySelect],ylim,'r--','LineW',0.5);
end

xlabel('Grid Y No. (x 0.8 mm)');


% Set y scale to mm
ytick_temp = 0:50:200;
set(gca,'ytick',ytick_temp);
set(gca,'yticklabel',ytick_temp/10);
ylabel('Depth (mm)');

% rectangle('Position',[ySelect+0.2,maxZ*0.95,0.6,10],'FaceColor','r','EdgeColor','r');

if hemisphere == 1
    set(gca,'xdir','rev');
end

for oldsize = 5:100
    set(findall(gcf,'fontsize',oldsize),'fontsize',13);
end
set(findall(gcf,'tickdir','i'),'tickdir','o');

drawnow;



%% Plot Tuning.  HH20140624
figure(803);
% --------------- Tuning Properties
% Mask: monkey & hemishpere & xLoc & significant visual tuning
% mask_tuning = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ...
%     & (num(:,xls.p_vis) < 0.05);
% to_plot = num(mask_tuning,:);
%
% for i =  1:size(to_plot,1)
%     offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
%     xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
%     yy = offSet + round(to_plot(i,xls.Depth)/100);
%
%     if to_plot(i,xls.Pref_vis) > 0
%         plot(xx,yy,'r>','linewid',1.2);
%     else
%         plot(xx,yy,'r<','linewid',1.2);
%     end
% %     plot(xx,yy,'ro' );
% end

% --------------- Memsac Properties
% % Mask: monkey & hemishpere & xLoc & significant visual tuning
% % mask_memsac = (num(:,xls.Monkey) == monkey) & (strcmp(txt(:,xls.Hemisphere),hemisphere)) & (num(:,xls.Xloc) == xSelect) ...
% %     & (num(:,xls.p_M) < 0.05);
% mask_memsac = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ...
%     & (1|num(:,xls.p_M) < 0.05) & (strcmp(txt(:,xls.Area),'LIP')) & (num(:,xls.Chan1)>0);
% to_plot = num(mask_memsac,:);
%
% for i =  1:size(to_plot,1)
%     offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
%     xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.pref_M) > 0) * 0.2 + (to_plot(i,xls.pref_M) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
%     yy = offSet + round(to_plot(i,xls.Depth)/100);
%     %
%     %     if to_plot(i,xls.pref_M) > 0
%     %         plot(xx,yy,'k>','linewid',1.2);
%     %     else
%     %         plot(xx,yy,'k<','linewid',1.2);
%     %     end
%     if to_plot(i,xls.p_M)<0.05
%         if to_plot(i,xls.pref_M) > 0
%             plot(xx,yy,'k>','linewid',1,'markerfacecol','g','markersize',6);
%         else
%             plot(xx,yy,'k<','linewid',1,'markerfacecol','g','markersize',6);
%             %            plot(xx,yy,'ko','linewid',1,'markerfacecol','k');
%         end
%
%     else
%         plot(to_plot(i,xls.Yloc),yy,'go','markersize',5);
%     end
% end

% --------------- MemSac Properties (MU & SU)
% Mask: monkey & hemishpere & xLoc & significant visual tuning
%{
mask_tuning = (num(:,xls.Monkey) == monkey) & num(:,xls.Hemisphere)==hemisphere & ...
    (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) ;
to_plot_num = num(mask_tuning,:);
to_plot_txt = txt(mask_tuning,:);

for i =  1:size(to_plot_num,1)
    
    if strcmp(to_plot_txt(i,xls.Protocol),'MemSac')
        offSet = round((to_plot_num(i,xls.guidetube)  + to_plot_num(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
        xx = to_plot_num(i,xls.Yloc); % + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
        yy = offSet + round(to_plot_num(i,xls.Depth)/100);
        
        if to_plot_num(i,xls.p_M) < 0.01 %to_plot_num(i,xls.HD_MemSac) >= 0.8  % T site
            plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
            pref_M = headingToAzi(to_plot_num(i,xls.pref_M))*pi/180;
            aspect = daspect;
            plot([xx xx+0.5*cos(pref_M)*sign((hemisphere==2)-0.5)],[yy yy+(-1)*0.5*sin(pref_M)*aspect(2)/aspect(1)],'r-','linew',1.5);
        else  % non-T site
            plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',5);
        end
    end
    
end
%}

% --------------- LIP SUs -------
% % Mask: monkey & hemishpere & xLoc & significant visual tuning
% %{

mask_tuning = (strcmp(txt(:,xls.Protocol),'HD') | strcmp(txt(:,xls.Protocol),'HD_dt')) & (strcmpi(txt(:,xls.Area),'LIP') | strcmpi(txt(:,xls.Area),'LIP-V')) ...
    & (num(:,xls.HD_rep) >= 8) & (num(:,xls.Units_RealSU) == 1) & num(:,xls.Monkey) == monkey & num(:,xls.Hemisphere)==hemisphere...
    & (num(:,xls.Xloc) >= xSelect + overlapping(2,1) & num(:,xls.Xloc) <= xSelect + overlapping(2,2)) & (num(:,xls.HD_TargFirst) == 1);

to_plot = num(mask_tuning,:);
text(max(xlim)*0.6, max(ylim),sprintf('No. cells = %g',sum(mask_tuning)));

for i =  1:size(to_plot,1)
    offSet = round((to_plot(i,xls.guidetube)  + to_plot(i,xls.offset) - 2.0) * 100);  % Related to Guide Tube 2.0 cm!!
    %    xx = to_plot(i,xls.Yloc) + ((to_plot(i,xls.Pref_vis) > 0) * 0.2 + (to_plot(i,xls.Pref_vis) <= 0) * -0.2)*sign((hemisphere==2)-0.5);  % Left and Right Prefer
    xx = to_plot(i,xls.Yloc) + (rand(1) - 0.5)*0;  % Left and Right Prefer
    yy = offSet + (to_plot(i,xls.Depth)/100);

    % Plot according to the max choice divergence (if all non-significant, black circles)
    max_color = {'b','r','g'};
    choice_pref_this = to_plot(i,xls.HD_vest_ChoicePref : xls.HD_comb_ChoicePref_p);
    [~, max_ind] = max(choice_pref_this([1 3 5]));
    any_sig = any(choice_pref_this([2 4 6]) < 0.05);
    
    if any_sig
        plot(xx,yy,['ko'],'linewid',1,'markerfacecol',max_color{max_ind},'markersize',7);
    else
        plot(xx,yy,'ko','linewid',1,'markerfacecol','none','markersize',7);
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',7);
    end

        
%     if to_plot(i,xls.HD_MemSac) >= 0.9 % T cell
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
%     else % non-T cell
%         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','r','markersize',5);
% %         plot(xx,yy,'ro','linewid',0.4,'markerfacecol','none','markersize',7);
%     end
end
%}

function ShowDepth(~,~,offset)
persistent depthLine;
pos = get(gca,'CurrentPo');
depth = pos(1,2);

try
    set(depthLine(1),'ydata',[depth depth]);
    set(depthLine(2),'position',[max(xlim)*0.8 depth],'string',sprintf('%0.0f\n',(depth-offset)*100));
catch
    depthLine(1) = plot(xlim,[depth depth],'b--','ButtonDownFcn',{@ShowDepth,offset});
    depthLine(2) = text(max(xlim)*0.8,depth,sprintf('%0.0f\n',(depth-offset)*100),'color','b','fontsize',15);
end


% function ReadXls()
%
% %% Read xls for plotting tuning. HH20140624
% global num txt raw xls;
% [num,txt,raw] = xlsread('Z:\Data\MOOG\Results\Result.xlsm',2);
%
% % Get Header infomation
% HEADS_N = 3;
%
% header_all = txt(HEADS_N-1,:);
% header_all(strcmp(header_all,'')) = txt(HEADS_N-2,strcmp(header_all,''));
%
% for i = 1:length(header_all)
%     try
%         if i == num(1,i)
%             eval(['xls.' header_all{i} '=' num2str(i) ';']);
%         else
%             disp('Header info error...');
%             keyboard;
%         end
%     catch
%     end
% end
%
% % Delete headers
% end_line = find(~isnan(num(:,1)),1,'last');
% num = num(HEADS_N+1:end_line,:);
% txt = txt(HEADS_N : end_line - 1,:); % Note here
% raw = raw(HEADS_N+1:end_line,:);

function change_monkey_hemi(~,~)
DrawMapping(get(gcbo,'value'));


% 20180604 Add manual align buttons
function manual_adjust_mri(source,event,choice)  
global MRI_offset_new hemisphere;

manual_step_size = 0.5;
manual_zoom_size = 1.02;

switch choice
    case {1,2}
        MRI_offset_new{1}(2) = MRI_offset_new{1}(2) + sign(choice - 1.5) * manual_step_size * 10;
    case {3,4}
        MRI_offset_new{1}(1) = MRI_offset_new{1}(1) + sign(3.5 - choice) * sign(1.5 - hemisphere) * manual_step_size;
    case {5,6}
        MRI_offset_new{1}(3) = MRI_offset_new{1}(3) * manual_zoom_size^(sign(5.5 - choice));
    case {7,8}
        MRI_offset_new{2}(2) = MRI_offset_new{2}(2) + sign(choice - 7.5) * manual_step_size/3;
end

fprintf('MRI_offset_new = {[%g, %g, %g],[%g, %g]};\n',...
        MRI_offset_new{1}, MRI_offset_new{2});
SelectChannel(source,event);
