
%% Draw figures

clc;close all;

%% Figure 2a
load data_marginally_stable

Msize  = 6;
lwidth = 1.2;
Fsize  = 12;

% line
figure;
     shadedErrorBar(Num*T, Err1, {@mean,@std}, 'lineprops', '-r'); hold on
     shadedErrorBar(Num*T, Err2, {@mean,@std}, 'lineprops', '-b');
     shadedErrorBar(Num*T, Err3, {@mean,@std}, 'lineprops', '-g');
h1 = semilogy(Num*T,mean(Err1),'r','linewidth',lwidth,'marker','o','markersize',Msize,'markerfacecolor','r');
h2 = semilogy(Num*T,mean(Err2),'b','linewidth',lwidth,'marker','d','markersize',Msize,'markerfacecolor','b');
h3 = semilogy(Num*T,mean(Err3),'g','linewidth',lwidth,'marker','s','markersize',Msize,'markerfacecolor','g');

% limits     
%ylim([10^(-2), 2*10^(-1)]);
%xlim([1000,10000]);
xlim([500,5000])
xtickformat('%1.0e');
set(gca,'TickLabelInterpreter','latex','FontSize',Fsize);
     
% label
xlabel('Number of Samples $T\times N$','Interpreter','latex','FontSize',Fsize);
ylabel('$\|\hat{G} - G\|/\|G\|$','Interpreter','latex','FontSize',Fsize);
% h = legend([h1,h2],'Multiple trajectories (all data points)','Multiple trajectories (final data point)',...
%    'Location','Northeast');
%set(h,'FontSize',Fsize,'Interpreter','latex','box','off')

% figure size
set(gcf,'Position',[250 150 300 350]);
print(gcf,'Fig2a','-painters','-dpng','-r600')


%% Figure 2b
load data_unstable
 % line
figure;
shadedErrorBar(Num*T, Err1, {@mean,@std}, 'lineprops', '-r'); hold on
shadedErrorBar(Num*T, Err2, {@mean,@std}, 'lineprops', '-b');
h1 = semilogy(Num*T,mean(Err1),'r','linewidth',lwidth,'marker','o','markersize',Msize,'markerfacecolor','r');
h2 = semilogy(Num*T,mean(Err2),'b','linewidth',lwidth,'marker','d','markersize',Msize,'markerfacecolor','b');

% limits     
%ylim([10^(-3), 1e-1]);
xlim([500,5000]);
xtickformat('%1.0e');
set(gca,'TickLabelInterpreter','latex','FontSize',Fsize);
     
% label
xlabel('Number of Samples $T\times N$','Interpreter','latex','FontSize',Fsize);
ylabel('$\|\hat{G} - G\|/\|G\|$','Interpreter','latex','FontSize',Fsize);

% figure size
set(gcf,'Position',[250 150 300 350]);
print(gcf,'Fig2b','-painters','-dpng','-r600')    