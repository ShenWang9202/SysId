% =================================================================
%
% Draw Fig. 2 in our paper
% Y. Zheng, N. Li, Non-asymptotic  Identification  of  Partially  Observable  
%                     Linear Time-invariant  Systems  using  Multiple  Trajectories
%
% =================================================================

clc;clear;close all

Markersize = 6;
linewidth  = 1.2;
Fontsize   = 10;

% First plot all
load data_marginally_stable_varyingT 
fig = figure;
for idx = 1:4 
    subplot(1,4,idx);
    shadedErrorBar(Tind, Err1{idx}, {@mean,@std}, 'lineprops', '-r'); hold on
    shadedErrorBar(Tind, Err2{idx}, {@mean,@std}, 'lineprops', '-b');
    shadedErrorBar(Tind, Err3{idx}, {@mean,@std}, 'lineprops', '-g');
    h1 = semilogy(Tind,mean(Err1{idx}),'r','linewidth',linewidth,'marker','o','markersize',Markersize,'markerfacecolor','r');
    h2 = semilogy(Tind,mean(Err2{idx}),'b','linewidth',linewidth,'marker','d','markersize',Markersize,'markerfacecolor','b');
    h3 = semilogy(Tind,mean(Err3{idx}),'g','linewidth',linewidth,'marker','s','markersize',Markersize,'markerfacecolor','g');

    set(gca,'Position',[0.08+(idx-1)*0.16+(idx-1)*0.07 0.27 0.16 0.7],...
        'TickLabelInterpreter','latex','fontsize',Fontsize)
    
    if idx == 1
       ylabel('$\|\hat{G} - G\|/\|G\|$','Interpreter','latex','FontSize',Fontsize);
    end
end

annotation('textbox',[.4 0 .5 .2],'String','Length of each experiment $T$',...
  'Interpreter','latex','FontSize',Fontsize,'EdgeColor','none')

% Set legend
h = legend([h1,h2,h3],'Multi-rollout (all data)','Multi-rollout (Sun et al.)',...
    'Single-rollout (Simchowitz et al.)');
set(h,'orientation','horizontal','FontSize',Fontsize,...
    'Position',[0.17 0.04 0.7 0.03],'Interpreter','latex')
set(h,'box','off');

% Set paper position for printing & print
set(gcf,'Position',[250 150 1000 200]);
print(gcf,'mstable_T','-painters','-dpng','-r600')
