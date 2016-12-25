% Raw data from ORIGIN
%Meanx Meany minx miny maxx maxy
%Left hemisphere raw
raw= [18.72727 -8.81818 13 -17 30	-1,
8.35714 -11.85714 3.5 -22 20 -7,
10.6 -14.4	7	-18	13	-11,
18.03571 -2.35714	5	-14	26	6,
3.3	-17.6 -2	-24	7	-1,
13	-32	12 -33	14	-31,
8.1	2.35 1	-16	11	14.5,
26.13636 -4.13636	19	-7	30	0.5,
5.91667	19 3.5	17.5	9	21.5,
21	2.5	20.5 1	21.5	4,
14	9	10	7 18	11];

figure;
hold on;
line([-1,1],[0,0]);
line([0,0],[-1,1]);
for i=1:11
    patch([raw(i,3),raw(i,1),raw(i,5),raw(i,1)],[raw(i,2),raw(i,6),raw(i,2),raw(i,4)],'y')
end
axis equal
% xlim([-5 35]);
% ylim([-35 35]);