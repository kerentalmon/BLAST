%% TOA localization scheme - Solving the Ambiguity problem of 2 receivers - SHARKS DATA
% with Foward Backward algorithm
clear
close all
warning off
Vs=1548;
t=readtable('TruePositionsWithTI.xlsx');
ID=table2array(t(:,1));
AccurateTime=table2array(t(:,2));
TOA1=table2array(t(:,3));
TOA2=table2array(t(:,5));
TOA3=table2array(t(:,7));
Anchor1=table2array(t(:,4));
Anchor2=table2array(t(:,6));
Anchor3=table2array(t(:,8));
Pos=table2array(t(:,9:10));
TI=table2array(t(:,11));
TI1=table2array(t(:,12));
TI2=table2array(t(:,13));
TI3=table2array(t(:,14));
%coordinates of receivers [North, East]
P1=[32+28/60+0.47/3600,  34+52/60+41.47/3600]; %555
P3=[32+27/60+48.37/3600, 34+52/60+36.67/3600]; %557
P4=[32+27/60+38.11/3600, 34+52/60+47.36/3600]; %558
P5=[32+27/60+35.95/3600, 34+52/60+32.05/3600]; %559
Pc = [P1; P3; P4; P5];
[P_x, P_y, Buoy_utmzone] = deg2utm(Pc(:,1), Pc(:,2));
%UTM of TBRs
P = [P_x, P_y];
Anchors=P';
L=8;   % number of Lobsters' positions
Lobsters=Pos(1:L,:)';
R1=((Lobsters(1,:)-Anchors(1,1)).^2+(Lobsters(2,:)-Anchors(2,1)).^2).^0.5; %555
R2=((Lobsters(1,:)-Anchors(1,2)).^2+(Lobsters(2,:)-Anchors(2,2)).^2).^0.5; %557
R3=((Lobsters(1,:)-Anchors(1,3)).^2+(Lobsters(2,:)-Anchors(2,3)).^2).^0.5; %558
R4=((Lobsters(1,:)-Anchors(1,4)).^2+(Lobsters(2,:)-Anchors(2,4)).^2).^0.5; %559
STAT=zeros(9,L); %register the simulation resulsts
rangeerror=zeros(9,1); %register the simulation resulsts - range error
VLmax=1;
%% Here we assume that Lobster #1 is localized thus its transmission time is 
% known from now on. The next Lobster positions are calcolated based on TOA
% scheme
syms xl yl       % Lobster's coordinates at times i+1, i+2,.....
S12=zeros(4,L);  % Lobster's position - two possible solutions
S=zeros(2,L);    % Lobster's position - after resolving ambiguity
S12(1,1)=Lobsters(1,1);
S12(2,1)=Lobsters(1,1);
S12(3,1)=Lobsters(2,1);
S12(4,1)=Lobsters(2,1);

for n=2:L
   eqns=([%sqrt((xl-Anchors(1,2))^2+(yl-Anchors(2,2))^2)/Vs==TOA1(n)-TOA1(n-1)-TI(n)+R2(n-1)/Vs,...
          sqrt((xl-Anchors(1,3))^2+(yl-Anchors(2,3))^2)/Vs==TOA2(n)-TOA2(n-1)-TI3(n)+R3(n-1)/Vs,...
          sqrt((xl-Anchors(1,1))^2+(yl-Anchors(2,1))^2)/Vs==TOA3(n)-TOA3(n-1)-TI3(n)+R1(n-1)/Vs...
             ]);
    [tempx,tempy]=solve(eqns,xl,yl);
   S12(:,n)=eval([tempx;tempy]);
end
%Rearrange S12 so that S1 is actual position and S2 is the imaginery
%position (this is for convinience)
 rearrS12=sqrt((S12(1,:)-Lobsters(1,1:L)).^2+(S12(3,:)-Lobsters(2,1:L)).^2)>...
     sqrt((S12(2,:)-Lobsters(1,1:L)).^2+(S12(4,:)-Lobsters(2,1:L)).^2); 
 for i=1:L
     if rearrS12(i)
         S12(:,i)=S12([2 1 4 3],i);
     end
 end
         
%% Prepare observation data for the forward backword algorithm

delta=zeros(4,length(S12)-1);
delta(1,:)=sqrt(sum((S12([1 3],2:end)-S12([1 3],1:end-1)).^2));
delta(2,:)=sqrt(sum((S12([1 3],2:end)-S12([2 4],1:end-1)).^2));
delta(3,:)=sqrt(sum((S12([2 4],2:end)-S12([2 4],1:end-1)).^2));
delta(4,:)=sqrt(sum((S12([2 4],2:end)-S12([1 3],1:end-1)).^2));
D=1;
p=forback1(D,delta,VLmax,TI(2:L));
pS1=double(p(1,2:end)>=p(2,2:end));
pS2=double(p(1,2:end)<p(2,2:end));
S(1,2:end)=S12(1,2:end).*pS1+S12(2,2:end).*pS2;
S(2,2:end)=S12(3,2:end).*pS1+S12(4,2:end).*pS2;
S(:,1)=Lobsters(:,1);
%% plot the results
minx=min([Anchors(1,:) S12(1,:) S12(2,:)])-100;
maxx=max([Anchors(1,:) S12(1,:) S12(2,:)]);
miny=min([Anchors(2,:) S12(3,:) S12(4,:)])-100;
maxy=max([Anchors(1,:) S12(3,:) S12(4,:)]);
anc1='  TBR 555';
anc2='  TBR 557';
anc3='  TBR 558';
anc4='  TBR 559';

figure(1)
hold on
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^g','MarkerSize',10,'MarkerFaceColor','g')
text(Anchors(1,1)-minx,Anchors(2,1)-miny,anc1,'FontSize',10)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2,'FontSize',10)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3,'FontSize',10)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4,'FontSize',10)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'.k','MarkerSize',15)
axis([0 800 0 900])
% title('Sharks Reception by 3 Anchors (TBR 555, 557, 558)')
legend('Anchors (TBR)', 'True Positions','Location', 'northwest')
xlabel('Position X-axis [m]', 'fontsize', 16)
ylabel('Position Y-axis [m]', 'fontsize', 16)

figure(2)
hold on
plot(S12(1,1)-minx,S12(3,1)-miny,'*k','MarkerSize',14)
plot(S12(1,2:L)-minx,S12(3,2:L)-miny,'sk','MarkerSize',14,'MarkerFaceColor','k')
plot(S12(2,2:L)-minx,S12(4,2:L)-miny,'dk','MarkerSize',14,'MarkerFaceColor','k')
% plot(S(1,2:L)-minx,S(2,2:L)-miny,'om','MarkerSize',16)
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^k','MarkerSize',14,'MarkerFaceColor','k')
text(Anchors(1,1)-minx-110,Anchors(2,1)-miny,anc1,'FontSize',14)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2,'FontSize',14)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3,'FontSize',14)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4,'FontSize',14)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'ok','MarkerSize',16)
axis([0 800 0 900])
% title('Bilatertion - Sharks Reception by 2 Anchors (TBR 555, 558)')
xlabel('Position X-axis [m]', 'fontsize', 18)
ylabel('Position Y-axis [m]', 'fontsize', 18)
set(gca,'FontSize',18, 'FontWeight','bold')
legend('First Shark', 'Sharks Solution 1','Sharks Solution 2', 'Receivers (TBR)',...
    'True Positions','Location', 'northeast', 'FontSize',16, 'FontWeight','bold')

figure(3)
hold on
plot(S12(1,1)-minx,S12(3,1)-miny,'*k','MarkerSize',14)
% plot(S12(1,2:L)-minx,S12(3,2:L)-miny,'sr','MarkerSize',10,'MarkerFaceColor','r')
plot(S12(2,2:L)-minx,S12(4,2:L)-miny,'dk','MarkerSize',14,'MarkerFaceColor','k')
plot(S(1,2:L)-minx,S(2,2:L)-miny,'ok','MarkerSize',14)
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^k','MarkerSize',14,'MarkerFaceColor','k')
text(Anchors(1,1)-minx-110,Anchors(2,1)-miny,anc1,'FontSize',14)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2,'FontSize',14)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3,'FontSize',14)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4,'FontSize',14)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'.k','MarkerSize',16)
axis([0 800 0 900])
% title('Bilatertion - Sharks Reception by 2 Anchors (TBR 555, 558)')
xlabel('Position X-axis [m]', 'fontsize', 18)
ylabel('Position Y-axis [m]', 'fontsize', 18)
set(gca,'FontSize',18, 'FontWeight','bold')
legend('First Shark','Ambiguous solution', 'BLAST - Estimated position','Receivers (TBR)',...
    'True Positions','Location', 'northeast','FontSize',16, 'FontWeight','bold')
   


figure(4)
hold on
plot(S12(1,1)-minx,S12(3,1)-miny,'or','MarkerSize',10,'MarkerFaceColor','b')
plot(S12(1,2:L)-minx,S12(3,2:L)-miny,'sr','MarkerSize',10,'MarkerFaceColor','r')
plot(S12(2,2:L)-minx,S12(4,2:L)-miny,'db','MarkerSize',10,'MarkerFaceColor','b')
plot(S(1,2:L)-minx,S(2,2:L)-miny,'om','MarkerSize',12)
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^g','MarkerSize',10,'MarkerFaceColor','g')
text(Anchors(1,1)-minx,Anchors(2,1)-miny,anc1,'FontSize',10)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2,'FontSize',10)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3,'FontSize',10)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4,'FontSize',10)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'.k','MarkerSize',15)
axis([0 800 0 900])
% title('Bilatertion - Sharks Reception by 2 Anchors (TBR 555, 558)')
legend('First Shark', 'Sharks Solution 1','Ambiguous solution', 'BLAS - Estimated position','Anchors (TBR)', 'True Positions','Location', 'northwest')
xlabel('Position X-axis [m]', 'fontsize', 16)
ylabel('Position Y-axis [m]', 'fontsize', 16)


%% workshop figures
figure('DefaultAxesFontSize',18)
set(0, 'defaultTextFontSize',18);
hold on
plot(S12(1,1)-minx,S12(3,1)-miny,'*k','MarkerSize',10)
plot(S12(1,2:L)-minx,S12(3,2:L)-miny,'sb','MarkerSize',10,'MarkerFaceColor','b')
plot(S12(2,2:L)-minx,S12(4,2:L)-miny,'dr','MarkerSize',10,'MarkerFaceColor','r')
% plot(S(1,2:L)-minx,S(2,2:L)-miny,'om','MarkerSize',16)
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^g','MarkerSize',10,'MarkerFaceColor','g')
text(Anchors(1,1)-minx,Anchors(2,1)-miny,anc1)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'ok','MarkerSize',13)
axis([0 800 0 900])
% title('Bilatertion - Sharks Reception by 2 Anchors (TBR 555, 558)')
legend('First Shark', 'Sharks Solution 1','Sharks Solution 2 (Ambiguity)', 'Anchors (TBR)', 'True Positions','Location', 'northwest')
xlabel('Position X-axis [m]')
ylabel('Position Y-axis [m]')

figure('DefaultAxesFontSize',18)
set(0, 'defaultTextFontSize',18);
hold on
plot(S12(1,1)-minx,S12(3,1)-miny,'*k','MarkerSize',10)
% plot(S12(1,2:L)-minx,S12(3,2:L)-miny,'sr','MarkerSize',10,'MarkerFaceColor','r')
plot(S12(2,2:L)-minx,S12(4,2:L)-miny,'dr','MarkerSize',10,'MarkerFaceColor','r')
plot(S(1,2:L)-minx,S(2,2:L)-miny,'ob','MarkerSize',10)
plot(Anchors(1,:)-minx,Anchors(2,:)-miny,'^g','MarkerSize',10,'MarkerFaceColor','g')
text(Anchors(1,1)-minx,Anchors(2,1)-miny,anc1)
text(Anchors(1,2)-minx,Anchors(2,2)-miny,anc2)
text(Anchors(1,3)-minx,Anchors(2,3)-miny,anc3)
text(Anchors(1,4)-minx,Anchors(2,4)-miny,anc4)
plot(Lobsters(1,:)-minx,Lobsters(2,:)-miny,'.k','MarkerSize',13)
axis([0 800 0 900])
% title('Bilatertion - Sharks Reception by 2 Anchors (TBR 555, 558)')
legend('First Shark','Sharks Solution 2 (Ambiguity) ', 'BLAS - Estimated position','Anchors (TBR)', 'True Positions','Location', 'northwest')
xlabel('Position X-axis [m]')
ylabel('Position Y-axis [m]')