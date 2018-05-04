%Plotting

maindir=pwd;
addpath(genpath(strcat(maindir,'\CodeDWT')))
addpath(genpath(strcat(maindir,'\Matlab_functions')))

%%%%% Load summary file
load('SummaryModel2.mat')
%%%%% Set up variables from above needed for plots 
X_age_unique=(20:90)';
X_pressure_unique=[7,10,15,20,25,30,35,40,45]';
X_age_all=kron(ones(length(X_pressure_unique),1),X_age_unique);
X_pressure_all=kron(X_pressure_unique,ones(length(X_age_unique),1));
n_pressure=length(X_pressure_unique);
n_age=length(X_age_unique);
X_pressure=kron(X_pressure_unique,ones(34,1));

AUC_wgts=[5/2,4,5,5,5,5,5,5,5/2];
ngrid=639;
indx=1:ngrid; 
X_age=model.X(:,2);
T_long=length(Theta);
T_lat=length(Phi);
latitude=repmat(Phi,1,T_long);
longitude=repmat(Theta',T_lat,1);
xx=radius*sin(latitude).*cos(longitude);
yy=-radius*sin(latitude).*sin(longitude);  %%% Note 4/27 -- reverse y for picture since heatmap plot has negative y up not down.
x_pos=reshape(xx,1,T_long*T_lat);       
y_pos=reshape(yy,1,T_long*T_lat);

ages=[20,30,40,50,60,70,80,90];
pressure=unique(X_pressure(:,1));
interval=6;
n1=floor(length(Phi)/interval);

% Create directories to store plots
mkdir('Plots')
mkdir('Plots\age_by_pressure_with_df')
mkdir('Plots\Combo_plots')
mkdir('Plots\auc_circumferential')
mkdir('Plots\AUC_dAUC_by_age_by_region')
mkdir('Plots\Intra_corr')
mkdir('Plots\Intra_IOP')
mkdir('Plots\AUCvsAGE')



%%% Figure 2 - age by pressure with df
%%% Supplementary material plot: Generate all the file to produce the movie MPSvsAge-wave.mp4


row=indx((X_age_all==90).*(X_pressure_all==45)==1);
B0=fhat.mean(row,:); 

figure('units','normalized','outerposition',[0 0 1 1])
colormap('jet')
set(gcf, 'Renderer', 'zbuffer');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6,'fontWeight','bold')


figure_no=0;
for i1=1:(120/interval)
    for i2=1:(120/interval)
        
    t=(i2-1)*6*120+(i1-1)*6+1;
    
    for i=1:n_pressure
         pressure_i=((i-1)*n_age+1):(i*n_age);
         subplot(3,6,i+3*ceil(i/3))
         plot(X_age_unique,fhat.mean(pressure_i,t),'LineWidth',2)
         title(['Pressure=',num2str(pressure(i))],'FontSize',6)
         min1=min(min(fhat.Q025(:,t)));
         max1=max(max(fhat.Q975(:,t)));
         %min2=min(Y_simulated(:,t));
         %max2=max(Y_simulated(:,t));
         %ylim([min(min1,min2),max(max1,max2)]);
         ylim([min1,max1]);
         xlim([20,90])
         hold on
         plot(X_age_unique,fhat.Q025(pressure_i,t),'--r','LineWidth',1)
         plot(X_age_unique,fhat.Q975(pressure_i,t),'--r','LineWidth',1)
         plot(X_age_unique,fhat.sQ025(pressure_i,t),'r','LineWidth',2)
         plot(X_age_unique,fhat.sQ975(pressure_i,t),'r','LineWidth',2)
         %plot(X_age(X_pressure==pressure(i)),Y_simulated(X_pressure==pressure(i),t),'k.')
         xlabel('Age','FontSize',6);
         ylabel('MPS','FontSize',6);
         hold off
     end;
    
    subplot(2,2,1)
    plot_eye(B0,radius,Phi,Theta,[0.5,1.5],2000),colorbar
    set(gca,'xtick',[]),set(gca,'ytick',[])
    title('Fitted MPS, 90 year old, 45 mmHg'); %note: actually intercept.
    hold on
    plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    %plot(x_pos(i)/12+500,y_pos(i)/12+500,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    hold off
    subplot(2,2,3)
    plot_eye(df_age_mean,radius,Phi,Theta,[2,4],2000),colorbar
    colormap('jet')
    set(gca,'xtick',[]),set(gca,'ytick',[])
    title('Posterior mean, df(\phi,\theta) of Age Effect');
    hold on
    plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    hold off
 
    figure_no=figure_no+1;
    fname=strcat(maindir,'\Plots\age_by_pressure_with_df\age_by_pressure_with_df',num2str(figure_no),'.png');
    saveas(gcf,fname);
    end;
end;

%%% Figure 3 - Combo Plot 
%%% Supplementary Combo plot.mp4: generates file to produce the movie
n_eyes=34;
age_raw=kron(1/9*ones(1,9),eye(n_eyes))*X_age;
unique_idx = n_eyes:n_eyes:size(model.X,1);
X_IOP_unique= model.X(unique_idx,3:4);
IOP_levels=[7,10,15,20,25,30,35,40,45]';
IOP_grid=7:45;
X_IOP_grid=[ones(size(IOP_grid))',spline(IOP_levels,X_IOP_unique(:,1),IOP_grid)',spline(IOP_levels,X_IOP_unique(:,2),IOP_grid)'];
AUC_wgts_all=kron(AUC_wgts,eye(n_eyes));
%AUC_raw=AUC_wgts_all*Y_simulated;

[Phi_1,Kj1] = Get_DWT('db3',120,'per       ',0,5);
[Phi_2,Kj2] = Get_DWT('db3',120,'reflection',1,5);
Phi_all = kron(Phi_1,Phi_2); % DWT matrix 
IDWT = kron((Phi_1'*Phi_1)\Phi_1',(Phi_2'*Phi_2)\Phi_2');
H = size(model.m,2);
T = 14400;
keep_reorder=wavespecs_compressed.keep;                    %%% non-thresholded parameters in MCMC order
keep_reorder(wavespecs_compressed.reorder)=keep_reorder;
flag=keep_reorder==1;

close
figure('units','normalized','outerposition',[0 0 1 1])
colormap('jet')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6,'fontWeight','bold')

tic
figure_no=0;
for i1=1:(120/interval)
    for i2=1:(120/interval)
        t=(i2-1)*6*120+(i1-1)*6+1;
        %%% Compute serial correlation across IOP and eye-to-eye marginal
        %%% variance by IOP and scleral position
        covP=X_IOP_grid*diag(mvar2(t,1:3))*X_IOP_grid';
        sigP=diag(covP);
        corP=diag(sigP.^(-1/2))*covP*diag(sigP.^(-1/2));
        %%% Cpompute intrafunctional correlation 
        inter_cor = zeros(T,H);
        for h=1:H
            tmp_sig = sig_star(h,flag)';
            for t2=1:T              
                    inter_cov = (IDWT(t,flag).*IDWT(t2,flag))*tmp_sig;
                    inter_cor(t2,h) = inter_cov/sqrt(mvar2(t,h)*mvar2(t2,h));          
            end
        end   
        cor_lim = [-1,1]; 
        subplot(2,3,1)
            set(gca,'xtick',[])
            plot_eye(B0,radius,Phi,Theta,[0.5,1.5],2000)
            set(gca,'xtick',[]),set(gca,'ytick',[])
            title('(a) Fitted MPS, 90yr, 45mmHg');%,'FontSize',16); 
            hold on
            plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
            hold off
        subplot(2,3,2)
            plot(X_age_unique,fhat.AUC_mean(:,t),'LineWidth',2)
            title('(b) AUC vs. age')
            ylim([0,100])
            xlim([20,90])
            hold on
            %plot(age_raw,AUC_raw(:,t),'k.')
            plot(X_age_unique,fhat.AUC_Q025(:,t),'--r','LineWidth',1)
            plot(X_age_unique,fhat.AUC_Q975(:,t),'--r','LineWidth',1)
            plot(X_age_unique,fhat.sAUC_Q025(:,t),'r','LineWidth',2)
            plot(X_age_unique,fhat.sAUC_Q975(:,t),'r','LineWidth',2)
            xlabel('Age')%,'FontSize',14);
            ylabel('AUC')%,'FontSize',14);
            hold off
        subplot(2,3,3)
            plot_eye(df_age_mean,radius,Phi,Theta,[2,4],2000)
            set(gca,'xtick',[]),set(gca,'ytick',[])
            %title('(e) df(\theta, \phi) of f(age, \theta, \phi)')%,'FontSize',14);
            title('(c) df(\theta,\phi)')
            hold on
            plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
            hold off
        subplot(2,3,4)
            imagesc(IOP_grid,IOP_grid,corP,[0.4,1]) % make x and y axis 7-45
            title('(d) Corr\{Y_{ip_1}(\theta,\phi),Y_{ip_2}(\theta,\phi)\}')%,'FontSize',14)
            xlabel('IOP')%^,'FontSize',16)
            ylabel('IOP')%,'FontSize',16)
        subplot(2,3,5)
            plot_eye(inter_cor(:,1),radius,Phi,Theta,cor_lim,2000)
            set(gca,'xtick',[]),set(gca,'ytick',[])
            title({'(e) Intrafunctional Correlation,';'Eye Intercept Q_1'})%,'FontSize',14); 
            hold on
            plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
            hold off
        subplot(2,3,6)
            plot_eye(inter_cor(:,4),radius,Phi,Theta,cor_lim,2000)
            colormap('jet')
            set(gca,'xtick',[]),set(gca,'ytick',[])
            title({'(f) Intrafunctional Correlation,';'Residual S'})%,'FontSize',14); 
            hold on
            plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
            hold off     
        figure_no=figure_no+1;
        fname=strcat(maindir,'\Plots\Combo_plots\Combo_plot_',num2str(figure_no),'.png');
        saveas(figureHandle,fname);
    end;
    i1;toc
end;
toc

%%% Figure 4 - Aggregated AUC
close
figure('units','normalized','outerposition',[0 0 1 1])
set(gcf, 'Renderer', 'zbuffer');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',5,'fontWeight','bold')

imagesc(fhat_circum_mean(71:-1:1,:));colorbar;
colormap('jet');
xlabel('Distance from ONH (meridional degree)');
ylabel('Age');
title('Aggregated AUC over circumferential regions');
set(gca,'XTick',[0.5,24:24:120]);
set(gca,'XTickLabels',(9:3:24));
set(gca,'YTick',[0.5,10:10:70]);
set(gca,'YTickLabels',[90:-10:30,20]);
%xticks([0.5,24:24:120]); if you have newer matlab versions
%xticklabels(9:3:24); if you have newer matlab versions
%yticks([0.5,10:10:70]); if you have newer matlab versions
%yticklabels([90:-10:30,20]); if you have newer matlab versions

fname=strcat(maindir,'\Plots\auc_circumferential\auc_circum.png');
saveas(gcf,fname);

close
%%% Figure 5


figure('units','normalized','outerposition',[0 0 1 1])

for i=1:2
    subplot(2,2,i)
    plot(X_age_unique,AUC_region.mean(:,i),'LineWidth',2)
    title({'Area under MPS vs. IOP curve';['Region ',region_labels{i}]})
    min1=min(min(AUC_region.Q025));
    max1=max(max(AUC_region.Q975));
    %min2=min(min(AUC_raw_region));
    %max2=max(min(AUC_raw_region));
    %ylim([min(min1,min2),max(max1,max2)]);
    ylim([min1,max1]);
    xlim([20,90])
    xlabel('Age')
    ylabel('AUC')
    hold on
    plot(X_age_unique,AUC_region.Q025(:,i),'--r','LineWidth',1)
    plot(X_age_unique,AUC_region.Q975(:,i),'--r','LineWidth',1)
    plot(X_age_unique,AUC_region.sQ025(:,i),'r','LineWidth',2)
    plot(X_age_unique,AUC_region.sQ975(:,i),'r','LineWidth',2)
    %plot(age_raw,AUC_raw_region(:,i),'k.')
    plot(X_age_unique,AUC_linear_region.mean(:,i),'k')
    hold off
    subplot(2,2,i+2)
    plot(X_age_unique,dAUC_region.mean(:,i),'LineWidth',2)
    title({'Derivative of Area under MPS vs. IOP curve';['Region ',region_labels{i}]})
    min1=min(min(dAUC_region.Q025));
    max1=max(max(dAUC_region.Q975));
    ylim([min1,max1])
    xlim([20,90])
    xlabel('Age')
    ylabel('dAUC/dAge')
    hold on
    plot(X_age_unique,dAUC_region.Q025(:,i),'--r','LineWidth',1)
    plot(X_age_unique,dAUC_region.Q975(:,i),'--r','LineWidth',1)
    plot(X_age_unique,dAUC_region.sQ025(:,i),'r','LineWidth',2)
    plot(X_age_unique,dAUC_region.sQ975(:,i),'r','LineWidth',2)
    plot(X_age_unique,dAUC_linear_region.mean(:,i),'k')
    plot(X_age_unique,zeros(size(X_age_unique)),'y')
    hold off
end;
fname=strcat(maindir,'\Plots\AUC_dAUC_by_age_by_region\AUC_dAUC_by_age_by_region.png');
saveas(gcf,fname);
close
%%% Supplementary files - Intrafunctional Correlation

figure('units','normalized','outerposition',[0 0 1 1])
colormap('jet')
set(gcf, 'Renderer', 'zbuffer');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6,'fontWeight','bold')

figure_no=0;
flag=(keep_reorder==1);
tic
for i1=1:(120/interval)
    for i2=1:(120/interval)
        t=(i2-1)*interval*120+(i1-1)*interval+1;
        inter_cor = zeros(T,H);
        for h=1:H
            tmp_sig = sig_star(h,flag)';
            for t2=1:T
                    inter_cov = (IDWT(t,flag).*IDWT(t2,flag))*tmp_sig;
                    inter_cor(t2,h) = inter_cov/sqrt(mvar2(t,h)*mvar2(t2,h));
            end
        end    
        cor_lim = [-1,1]; 
        subplot(2,2,1)
        plot_eye(inter_cor(:,1),radius,Phi,Theta,cor_lim,2000),colorbar
        title('(a) Eye Random Effect'); 
        hold on
        plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
        hold off
        
        subplot(2,2,2)
        plot_eye(inter_cor(:,2),radius,Phi,Theta,cor_lim,2000),colorbar
        title('(b) Random Slope for IOP1'); 
        hold on
        plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
        hold off
        
        subplot(2,2,3)
        plot_eye(inter_cor(:,3),radius,Phi,Theta,cor_lim,2000),colorbar
        title('(c) Random Slope for IOP2'); 
        hold on
        plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
        hold off
        
        subplot(2,2,4)
        plot_eye(inter_cor(:,4),radius,Phi,Theta,cor_lim,2000),colorbar
        title('(d) Residual'); 
        hold on
        plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',5)
        hold off
        
        figure_no=figure_no+1;
        fname=strcat(maindir,'\Plots\Intra_corr\interfuntional_cor_',num2str(figure_no),'.png');
        saveas(gcf,fname);
        
        clear inter_cor
      
    end
end
close
%%% Supplementary files - Intra-IOP-corr

figure('units','normalized','outerposition',[0 0 1 1])
colormap('jet')
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6,'fontWeight','bold')
tic
figure_no=0;
for i1=1:(120/interval)
    for i2=1:(120/interval)
        t=(i2-1)*6*120+(i1-1)*6+1;
        covP=X_IOP_grid*diag(mvar2(t,1:3))*X_IOP_grid';
        sigP=diag(covP);
        corP=diag(sigP.^(-1/2))*covP*diag(sigP.^(-1/2));
        subplot(2,2,2)  
            plot(IOP_grid,sigP,'LineWidth',2);
            %title('Eye-to-Eye Variability (IOP, \theta,\phi)')%,'FontSize',16)
            title('Var(Y_{ip}|IOP, \theta,\phi)')%,'FontSize',14)
            xlabel('IOP')%,'FontSize',16)
            ylabel('Var(Y_{ip})')%,'FontSize',16)
            %xlabel('IOP','FontSize',16)
            %ylabel('Q_{eye}(IOP,\theta,\phi)')%,'FontSize',16)
            ylim([0,0.5])
            xlim([7,45])
        subplot(2,2,4)
            imagesc(IOP_grid,IOP_grid,corP,[0.4,1]) % make x and y axis 7-45
            title('Corr\{Y_{ip_1}(\theta,\phi),Y_{ip_2}(\theta,\phi)\}')%,'FontSize',14)
            %title('Intra-IOP Correlation (\theta,\phi)')%,'FontSize',16)
            xlabel('IOP')%,'FontSize',16)
            ylabel('IOP')%,'FontSize',16)
        subplot(1,2,1)
            plot_eye(B0,radius,Phi,Theta,[0.5,1.5],2000)
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            title('Fitted MPS, 90 year old, 45 mmHg')%,'FontSize',16); %note: actually intercept.
            hold on
            plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
            hold off
        figure_no=figure_no+1;
        fname=strcat(maindir,'\Plots\Intra_IOP\Intra_IOP_corr_',num2str(figure_no),'.png');
        saveas(figureHandle,fname);
    end;
    i1,toc
end;
toc
close
%%% Supplementary files: AUCvsAge

row=indx((X_age_all==90).*(X_pressure_all==45)==1);
B0=fhat.mean(row,:); 

figure('units','normalized','outerposition',[0 0 1 1])
colormap('jet')
set(gcf, 'Renderer', 'zbuffer');
figureHandle = gcf;
set(findall(figureHandle,'type','text'),'fontSize',6,'fontWeight','bold');

tic
figure_no=0;
for i1=1:(120/interval)
    for i2=1:(120/interval)
  
    t=(i2-1)*6*120+(i1-1)*6+1;
    
    i=6;
    pressure_i=((i-1)*n_age+1):(i*n_age);
    subplot(1,2,2)
    plot(X_age_unique,fhat.mean(pressure_i,t),'LineWidth',2)
    title('Nonparametric Age Effect on MPS','FontSize',6)
    min1=min(min(fhat.Q025(:,t)));
    max1=max(max(fhat.Q975(:,t)));
    %min2=min(Y_simulated(:,t));
    %max2=max(Y_simulated(:,t));
    %ylim([min(min1,min2),max(max1,max2)]);
    ylim([min1,max1]);
    xlim([20,90])
    hold on
    plot(X_age_unique,fhat.Q025(pressure_i,t),'--r','LineWidth',1)
    plot(X_age_unique,fhat.Q975(pressure_i,t),'--r','LineWidth',1)
    plot(X_age_unique,fhat.sQ025(pressure_i,t),'r','LineWidth',2)
    plot(X_age_unique,fhat.sQ975(pressure_i,t),'r','LineWidth',2)
    %plot(X_age(X_pressure==pressure(i)),Y_simulated(X_pressure==pressure(i),t),'k.')
    xlabel('Age','FontSize',6);
    ylabel('MPS','FontSize',6);
    hold off
    
    subplot(2,2,1)
    plot_eye(B0,radius,Phi,Theta,[0.5,1.5],2000),colorbar
    %axes('Color','none','XColor','none');
    title('Fitted MPS, 90 year old, 30 mmHg'); %note: actually intercept.
    
    hold on
    plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    %plot(x_pos(i)/12+500,y_pos(i)/12+500,'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    hold off
    subplot(2,2,3)
    plot_eye(df_age_mean,radius,Phi,Theta,[2,5],2000),colorbar
    title('Posterior mean, df(\phi,\theta) of Age Effect');
    hold on
    plot(x_pos(t),y_pos(t),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',10)
    hold off
    
    figure_no=figure_no+1;
    fname=strcat(maindir,'\Plots\AUCvsAGE\AUCvsAGE_age90p45_',num2str(figure_no),'.png');
    saveas(gcf,fname);
    end;
    i1,toc
end;