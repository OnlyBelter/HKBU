

%% Load gut samples from HMP, genus taxonomoc level, relative abundance, unassigned reads are excluded
A=load('HMP_genus.txt');

%% normalize to relative abundances
% [NumSpecies,NumSamples]=size(A);
% A=A./repmat(sum(A),NumSpecies,1);

%% make null model
A_null=func_make_null(A,1);

%% get Overlap and Dissimilarity for the real samples
[Overlap,RootJSD]=func_Cal_Overlap_rJSD_from_relative_abundance(A);

%% get Overlap and Dissimilarity for the null model
[Overlap_null,RootJSD_null]=func_Cal_Overlap_rJSD_from_relative_abundance(A_null);

%% make Bootstrap for the real samples
N_times=10;
[BS, xs]=func_cal_rlowess_bootstrap(Overlap,RootJSD,N_times);

%% make Bootstrap for the null model
[BS_null, xs_null]=func_cal_rlowess_bootstrap(Overlap_null,RootJSD_null,N_times);

%% plot figure
c14_7=[  0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

c14_light_02=1-0.5*(1-c14_7);

figure('position',[138         156        1032         480])

subplot(1,2,1)
hold on

% real
mean_BS=mean(BS,1);
Prctile_BS = prctile(BS,[3  97],1)'; 
h_area=area(xs',[Prctile_BS(:,1), (Prctile_BS(:,2)-Prctile_BS(:,1)) ]);
h_area(1).FaceColor = 'none';
h_area(2).FaceColor = c14_light_02(6,:);
h_area(1).EdgeColor = 'none';
h_area(2).EdgeColor = 'none';

h_point=plot(Overlap(:),RootJSD(:),'.','color',0.5*[1 1 1],...
    'MarkerSize',14);

plot(xs,mean_BS','linewidth',3,'color',c14_7(1,:))
plot(xs,Prctile_BS(:,1),'linewidth',2,'color',c14_7(6,:))
plot(xs,Prctile_BS(:,2),'linewidth',2,'color',c14_7(6,:))

xc=detect_negative_slope_Lowess_03(xs,mean_BS');
plot(xc*[1 1],[0 1],'-g')


xlim([0.5 1])
ylim([0.05 0.6])


box on
axis square

set(gca,'fontsize',18)
set(gca,'xtick',0:0.1:1)
set(gca,'ytick',0.1:0.1:0.5)

xlabel('Overlap','fontsize',22)
ylabel('Dissimilarity','fontsize',22)
title({'HMP, stool, genus, '; [num2str(size(BS,1)) ' bootstrap realization(s)']},...
    'fontsize',22)

subplot(1,2,2)
hold on

% null
mean_BS_null=mean(BS_null,1);
Prctile_BS_null = prctile(BS_null,[3  97],1)'; 
h_area_null=area(xs_null',[Prctile_BS_null(:,1), (Prctile_BS_null(:,2)-Prctile_BS_null(:,1)) ]);
h_area_null(1).FaceColor = 'none';
h_area_null(2).FaceColor = c14_light_02(3,:);
h_area_null(1).EdgeColor = 'none';
h_area_null(2).EdgeColor = 'none';

h_point_null=plot(Overlap_null(:),RootJSD_null(:),'.','color',0.5*[1 1 1],...
    'MarkerSize',14);


plot(xs_null,mean_BS_null','linewidth',3,'color',c14_7(2,:))
plot(xs_null,Prctile_BS_null(:,1),'linewidth',2,'color',c14_7(3,:))
plot(xs_null,Prctile_BS_null(:,2),'linewidth',2,'color',c14_7(3,:))


xlim([0.5 1])
ylim([0.05 0.6])

box on
axis square

set(gca,'fontsize',18)
set(gca,'xtick',0:0.1:1)
set(gca,'ytick',0.1:0.1:0.5)

xlabel('Overlap','fontsize',22)
ylabel('Dissimilarity','fontsize',22)

title({'Randomized samples'; [num2str(size(BS,1)) ' bootstrap realization(s)']},...
    'fontsize',22)


pause(1)
h_marker=h_point_null.MarkerHandle;  h_marker.EdgeColorData(4)=100;
h_marker=h_point.MarkerHandle;  h_marker.EdgeColorData(4)=100;



%%
% Fraction of data with negative slope
Fns=sum(Overlap(:)>xc)/sum(~isnan(Overlap(:)));

% P-value for the slope to the right of the changing point
slope_BS_real=nan(N_times,1);
slope_BS_null=nan(N_times,1);
for k=1:N_times
    P_real=polyfit(xs(xs>xc),BS(k,xs>xc),1);
    slope_BS_real(k)=P_real(1);

    P_null=polyfit(xs_null(xs>xc),BS_null(k,xs>xc),1);
    slope_BS_null(k)=P_null(1);
end
[~,P_value_slopes]=ttest2(slope_BS_real,slope_BS_null);

% Display the universality scores
disp(['Fns=' num2str(Fns)])
disp(['p-value=' num2str(P_value_slopes)])



