% Parse data from datafile
data = dlmread('out_file');
dim = sqrt(size(data,1));

% Do some magic on data
X=reshape(data(:,1),dim,dim);
Y=reshape(data(:,2),dim,dim);
Z=reshape(data(:,3),dim,dim);

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'PlotBoxAspectRatio',[1 1 1],'FontSize',12);
view(axes1,[-127.5 16]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create mesh
mesh(X,Y,Z,'Parent',axes1);

% Create title
title('Variational energy','FontSize',20);

% Create xlabel
xlabel('\alpha','HorizontalAlignment','right','FontSize',16);

% Create ylabel
ylabel('\beta','FontSize',16,'HorizontalAlignment','left');

% Create zlabel
zlabel('E','FontSize',16);

% Create colorbar
colorbar('peer',axes1);

% Find the minimum energy
[E_min, idx] = min(data(:,3));
sprintf('Minimum energy: %f +/- %f\nMinimun alpha: %f\nMinimum beta: %f',E_min,data(idx,4),data(idx,1),data(idx,2))

% Insert some information on the plot
annotation(figure1,'textbox',...
    [0.0588806386390572 0.788701054413174 0.196193265007321 0.10752688172043],...
    'String',{...
        strjoin({'Minimum energy:',num2str(E_min),'+/-',num2str(data(idx,4))}),...
        strjoin({'Minimum alpha: ',num2str(data(idx,1))}),...
        strjoin({'Minimum beta: ',num2str(data(idx,2))})...
    },...
    'FitBoxToText','on',...
    'FontSize',12);


% Do some cleaning
clear('dim','X','Y','Z','E_min','idx');