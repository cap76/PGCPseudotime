D1 = importdata('/AllProcesses_TSCAN_data.xlsx')
D1.data = D1.data(:,2:end);
D1.data(3,:) = double(int64(D1.data(3,:)));
%D1 = importdata('./PGC2Data.txt')
%t  = D1.data(1,:);
%D1.data = D1.data(2:end,:);
%D1.data(3,:) = double(int64(D1.data(3,:)));
%l1 = importdata('Pseudotime_time_ED_ESC.csv');
%l2 = importdata('Pseudotime_time_ED_PGC.txt');
%l2 = importdata('Pseudotime_time_ED_ESC_PGC_reversed.csv');
%l3 = importdata('Pseudotime_time_ED_Soma.csv');
%l2 = importdata('Wishbone_pgs_soma_branches.csv');
%l2 = importdata('./TSCAN_Pseudotime_time_PGC.csv')
l2 = importdata('./Wishbone_pgc_pseudotimes.csv')
l3 = importdata('./Wishbone_soma_pseudotimes.csv')

%for i = 1:length(l1.data)
%    inds1(i,1) = find(strcmp(l1.textdata{i+1,1}, D1.textdata(6,2:end))==1);
%end

for i = 1:length(l2.data)
    inds2(i,1) = find(strcmp(l2.textdata{i,1}, D1.textdata(1,2:end))==1);
end

for i = 1:length(l3.data)
    inds3(i,1) = find(strcmp(l3.textdata{i,1}, D1.textdata(1,2:end))==1);
end


T1 = zeros(1,size(D1.data,2));
T2 = zeros(1,size(D1.data,2));

T1(inds2) = l2.data(:,1); 
T2(inds3) = l3.data(:,1); 
T1(find(T1==0)) = T2(find(T1==0));
Vt = ones(1,length(T1));

uT = unique(D1.data(3,:));%unique(Output.Data.orig.origt)
for i = 1:length(T1)
    Ind(i,1) = find( D1.data(3,i) == uT);
end

subplot(2,2,1);
inds = find((D1.data(1,:)~=2 ) );
%inds = find((D1.data(1,:)==0 | D1.data(1,:)==-1 ) );
vs = violinplot(T1(inds)+randn(1,length(inds))*0.01,Ind(inds));
set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14],'XTickLabel',{'Ooc.','Zy.','2C','4C','8C','Mor.','Blast','4W','7W','8W','10W','11W','17W','19W'},'YTick',[])
ylim([-0.2 1.2])
ax = gca;
ax.XTickLabelRotation = 30; 
ax.FontSize=11;
ylabel('Pseudotime')

subplot(2,2,2);
inds = find((D1.data(1,:)~=1 ));
vs = violinplot(T1(inds)+randn(1,length(inds))*0.01,Ind(inds));
ylim([-0.2 1.2])
set(gca,'XTick',[1,2,3,4,5,6,7,8,9,10,11,12,13,14],'XTickLabel',{'Ooc.','Zy.','2C','4C','8C','Mor.','Blast','4W','7W','8W','10W','11W','17W','19W'},'YTick',[])
ax = gca;
ax.XTickLabelRotation = 30; 
ax.FontSize=11;
ylabel('Pseudotime')



tnew = D1.data(3,:);
tnew = (tnew - min(tnew))./max(tnew - min(tnew))
T1 = (T1-min(T1))./max((T1-min(T1)));

uT2 = unique(tnew);
DELTA = zeros(length(uT2),3);
for i = 1:length(uT2)
    inds1 = find((D1.data(1,:)==0 | D1.data(1,:)==1 | D1.data(1,:)==-1) & tnew==uT2(i));
    inds2 = find((D1.data(1,:)==0 | D1.data(1,:)==2 | D1.data(1,:)==-1) & tnew==uT2(i));
    DELTA(i,1)= uT2(i);
    DELTA(i,2)= mean(T1(inds1));
    DELTA(i,3)= mean(T1(inds2));
end

uT = unique(D1.data(3,:));%unique(Output.Data.orig.origt)
for i = 1:length(D1.data(3,:))
    Ind(i,1) = find( D1.data(3,i) == uT);
end
T = T1;%Output.Param.Store{257}.update.x(:,1);

inds1 = find((D1.data(1,:)==0 | D1.data(1,:)==1) & Vt<=1);
inds2 = find((D1.data(1,:)==0 | D1.data(1,:)==2) & Vt<=1);

Summary(1,1)=sum((DELTA(:,2)-DELTA(:,3)).^2);
Summary(2,1)=corr(Ind(inds1),T(inds1)');
Summary(3,1)=corr(Ind(inds2),T(inds2)');
Summary(4,1)=corr(Ind(inds1),T(inds1)','Type','Spearman');
Summary(5,1)=corr(Ind(inds2),T(inds2)','Type','Spearman');





