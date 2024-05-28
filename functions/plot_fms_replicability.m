
%%% Replicability plot of a model 

function plot_fms_replicability(FMS_cp_real)

I=length(FMS_cp_real);
FMS_cp_real_allR = [];
g = [];
for R=1:I
    FMS_cp_real_allR = [FMS_cp_real_allR;FMS_cp_real{R}];
    g = [g;repmat(['R=',num2str(R)],size(FMS_cp_real{R}))];
end

%% plot
figure
boxplot(FMS_cp_real_allR,g); 
ylabel('FMS-CP-Real'); xlabel('Number of Components')
% set(gca,'XTick', 1:1:I, 'XTickLabel',{'$\textbf{R=1}$', '$\textbf{R=2}$','$\textbf{R=3}$','$\textbf{R=4}$'},'TickLabelInterpreter','latex')
set(gca,'Fontsize',18)
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.5);
for R=1:I
    index = round(length(FMS_cp_real{R})*0.95);
    aaa   = sort(FMS_cp_real{R},'descend');
    hold on;
    plot(R-0.25:0.005:R+0.25, ones(101,1)*aaa(index),'g-','Linewidth',2);
end
yticks([0.2 0.4 0.6 0.7 0.8 0.9 1.0])
ylim([0 1.0])
