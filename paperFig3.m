
%% Draw figures

clc;close all;


load data_marginally_stable_varyingT


Msize  = 6;
lwidth = 1.2;
Fsize  = 10;

% line
for i = 1:5
    figure;
    shadedErrorBar(Tind, Err1{i}, {@mean,@std}, 'lineprops', '-r'); hold on
    shadedErrorBar(Tind, Err2{i}, {@mean,@std}, 'lineprops', '-b');
    shadedErrorBar(Tind, Err3{i}, {@mean,@std}, 'lineprops', '-g');
    h1 = semilogy(Tind,mean(Err1{i}),'r','linewidth',lwidth,'marker','o','markersize',Msize,'markerfacecolor','r');
    h2 = semilogy(Tind,mean(Err2{i}),'b','linewidth',lwidth,'marker','d','markersize',Msize,'markerfacecolor','b');
    h3 = semilogy(Tind,mean(Err3{i}),'g','linewidth',lwidth,'marker','s','markersize',Msize,'markerfacecolor','g');
    set(gca,'TickLabelInterpreter','latex');

    % label
    xlabel('Length of each experiment $T$','Interpreter','latex','FontSize',Fsize);
    ylabel('$\|\hat{G} - G\|/\|G\|$','Interpreter','latex','FontSize',Fsize);

  %  h = legend([h1,h2,h3],'Multiple trajectories (all data points)','Multiple trajectories (Sun et al.)',...
  %      'Single trajectory (filtering)','Location','Northeast');
  %  set(h,'FontSize',Fsize,'Interpreter','latex','box','off')

    % figure size
    set(gcf,'Position',[250 150 300 320]);
    figName = ['mstable_t',num2str(i)];
    print(gcf,figName,'-painters','-dpng','-r600')
    %imwrite(gcf, 'a.png', 'png', 'transparency') 

end
     