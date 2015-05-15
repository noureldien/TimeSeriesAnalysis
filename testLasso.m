load Data_TSReg4

[lassoBetas,lassoInfo] = lasso(X0,y0,'CV',10);

[hax,hfig] = lassoPlot(lassoBetas,lassoInfo,'PlotType','Lambda');
hax.XGrid = 'on';
hax.YGrid = 'on';
hax.GridLineStyle = '-';
hax.Title.String = '{\bf Lasso Trace}';
hax.XLabel.String = 'Regularisation Parameter';
hax.YLabel.String = 'Feature Weights';
hlplot = hax.Children;
hMSEs = hlplot(5:6);
htraces = hlplot(4:-1:1);
set(hlplot,'LineWidth',2);
set(hMSEs,'Color','m');
legend(htraces,predNames0,'Location','NW');
hfig.HandleVisibility = 'on';






