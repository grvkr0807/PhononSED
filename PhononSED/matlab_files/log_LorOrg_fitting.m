
clear
clc

%%%%%%%%%%%%%v%% Constants %%%%%%%%%%%
c_light= 2.99792458e10; %speed of light

%%%%%%%%%%%%% Inputs, Config %%%%%%%%%%%%%
fraction= 2e-2;    % fraction=0.1 means any value of SED less than
%  0.1*(peak_SED) will not be considered for fitting
fraction_str= '2e-2';
material='RDX_lammps_1ns_T300';      %file header based on Sim. conditions  
type= 'log_LorOrg_fit';    %fitting algorithm 
Natomsunitcell= 168;    % 2 for silicon
Nunitcells= 1; % no of unit cells
Nk= 1;    % no. of kpoints for which SED,eigV etc have been calculated
Natoms= Natomsunitcell*Nunitcells;
Neig= Natomsunitcell*3;   %no. of modes/eigenstate per kpoint
Nmodes= Neig*(Nk);    %Ignoring 3 acoustic modes at gamma point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Beginning to import SED and frequency data \n');
% m= 3*Natoms; % no. of branches

SED_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'_SED.dat');
data = importdata (SED_file);   %#1 or 2 represents index of k-point
[f,junk]= size(data);   %no. of discrete freqs and no. of modes


Tau_total= zeros(Nk,Neig);
w_mid_rad_total= zeros(Nk,Neig);
SED_fit_total= zeros(Nk,f,Neig);
index_total=zeros(Nk,2,Neig);   %index for lower limit and upper limit of SED/frew window
freqs_total= zeros(Nk,f,1);
SED_rms= zeros(Nk,Neig);
plot_name1_total= strings(Nk,1);


% Tau_total= zeros(Nmodes,1);
% w_mid_rad_total= zeros(Nmodes,1);
% SED_fit_total= zeros(f,Nmodes);
% index_total=zeros(Nmodes,2);   %index for lower limit and upper limit of SED/frew window
% freqs_total= zeros(Nmodes,f);
% plot_name1_total= strings(Nmodes,1);

kpoints_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'_kvectors.dat');
kpoints = importdata (kpoints_file);   %#1 or 2 represents index of k-point


for kk= 1:Nk
    
    k_no= num2str(kk); %ith kpoint convert to string for file naming
    SED_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_SED.dat');
    data = importdata (SED_file);   %#1 or 2 represents index of k-point
    frequency_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_frequencies.dat');
    w0 = importdata (frequency_file);   %#1 or 2 represents index of k-point
    w0= w0*c_light*2*pi;
    w0(1:3)=[0;0;0];   % zero frequency for acoustic modes
    fprintf('Done importing SED and frequency data \n');

    freqs= data(:,1);   % first row is frequency, domain of frequency
    SED= data(:,2:end); 
    SED(SED<1e-6)= 1e-6;
    
    log_SED= zeros(f,Neig);  
    SED_fit= zeros(f,Neig);
    gamma= zeros(1,Neig);  % to store gamma value, tau= 1/(2*gamma)
    max_SED= zeros(1,Neig);    % to store peak SED
    max_i= zeros(1,Neig);      % to store index of peak SED
    peak_freqs= zeros(1,Neig);  % to store peak frequency
    w_mid= zeros(1,Neig);  
    resnormal= zeros(1,Neig);  % normalized residual after fitting for each modes
    residual= zeros(f,Neig);   % for each SED value in every mode
    params= zeros(3,Neig);  
    norm_error= zeros(1,Neig);     % to store norm error of fitting for each mode
    index= zeros(2,Neig);
    ind_mid= zeros(2,Neig);
    

for i=1:Neig

        [max_SED(i), max_i(i)] = max(SED(:,i));    %max SED and index for ith mode
        peak_freqs(i)= freqs(max_i(i)); 
        [~, index(1,i)]= min(abs(SED(1:max_i(i),i)- (fraction)*max_SED(i)));   %this outputs the index of frequency at which SED is close to 10%of its value, lower limit
        [~, index(2,i)]= min(abs(SED(max_i(i):f,i)- (fraction)*max_SED(i)));   %this outputs the index of frequency at which SED is close to 10%of its value, upper limit
        [~, ind_mid(1,i)]= min(abs(SED(1:max_i(i),i)- (0.5)*max_SED(i)));   %this outputs the index of frequency at which SED is half of maxima
        [~, ind_mid(2,i)]= min(abs(SED(max_i(i):f,i)- (0.5)*max_SED(i)));   %this outputs the index of frequency at which SED is half of maxima
        index(2,i)= max_i(i) + index(2,i) - 1;
        ind_mid(2,i)= max_i(i) + ind_mid(2,i) - 1;
        log_SED(index(1,i):index(2,i),i)= log(SED(index(1,i):index(2,i),i));
        
        SED_rms(kk,i)= rms(log(SED(index(1,i):index(2,i),i)));
        P2= peak_freqs(i);               
        P3= 0.5*(freqs(ind_mid(2,i)) - freqs(ind_mid(1,i)));
        P1= (max_SED(i))*P3^0.5;
        P0= [P1,P2,P3];    %parameters for Lorentzian fit
        BOUNDS= [0,0,0;Inf,Inf,Inf];
        [SED_fit(index(1,i):index(2,i),i), params(:,i), resnormal(i), residual(index(1,i):index(2,i),i)] = SEDfit_func(freqs(index(1,i):index(2,i)),log_SED(index(1,i):index(2,i),i),P0,BOUNDS,'loglororg');
        norm_error(i)= rms(log_SED(index(1,i):index(2,i),i) - SED_fit(index(1,i):index(2,i),i))/ SED_rms(kk,i);
        [a,b]=max(SED_fit(:,i));
        w_mid(i)= freqs(b);     %frequency in cm-1
        

end

% % Calculation of gamma using interpolation of SED_fit 
% abs_SED_fit= exp(SED_fit);
% for i= 1:m
%     a= 0.5*max(abs_SED_fit(:,i));
%     ind= find(abs_SED_fit(:,i)>=a);
%     nl= min(ind);
%     nr= max(ind);
%     fl = (freqs(nl) - freqs(nl-1))*(a - abs_SED_fit(nl,i))/(abs_SED_fit(nl,i) - abs_SED_fit(nl-1,i)) + freqs(nl);
%     fr = (freqs(nr) - freqs(nr+1))*(a - abs_SED_fit(nr,i))/(abs_SED_fit(nr,i) - abs_SED_fit(nr+1,i)) + freqs(nr);
%     gamma(i)= c_light*abs(fr-fl)/2;
% end

plot_name1=strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/','SED plots_',material,'_',type,'_','N',num2str(Nunitcells),'_','k',k_no,'fraction',fraction_str); % to be used by plot_gen

gamma= c_light*(abs(params(3,:)).^0.5); % gamma is value of P3, unit converted
w_mid_rad_total(kk,:)= w_mid*c_light*2*pi;
w_mid_rad_total(1,1:3)= 1e-6*[1 1 1];
Tau_total(kk,:)= 1./(2*gamma);     %lifetime in s
Tau_total(1,1:3)= max(Tau_total(1,4:end))*[1 1 1];  % zero tau for acoustic modes at gamma point
SED_fit_total(kk,:,:)= SED_fit;
freqs_total(kk,:,:)= freqs;
index_total(kk,:,:)= index;
plot_name1_total(kk,:)= plot_name1;


% Plot of Tau vs Frequency
plot_name=strcat('Tau vs Freq',{' '},material,type,{' '},',',num2str(fraction),{' '},'N',num2str(Nunitcells),{' '},'Nk',num2str(Nk),{' '},'k',k_no);
figure;
set( gca,'FontName', 'Helvetica' );
% loglog(w_mid_rad_total(kk,:),Tau_total(kk,:),'o','MarkerSize',6,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.0 0.9 0.9]);
loglog(w0,Tau_total(kk,:),'o','MarkerSize',6,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.0 0.9 0.9]);
title(plot_name,'Interpreter','none')
axis([1e12,1e15,1e-13,1e-10]);
     xlabel('Frequency (rad/s)');
     ylabel('Lifetimes (s)');
     

error_total= rms(norm_error)

filename=strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_',type,'_','_','fraction',fraction_str,'_','N',num2str(Nunitcells),'_','k',k_no);

% Save the plot
saveas(gcf,strcat(filename,'_Tau_vs_Omega.fig'))

end

x= w_mid_rad_total';
y= Tau_total';

% filename=strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_',type,'_','ntest',num2str(ntest),'_','fraction',fraction_str,'_','N',num2str(Nunitcells),'_','k',num2str(Nk));
% 
% % Save Tau data in h5 format
% h5create(strcat(filename,'_Tau_total.h5'),'/Tau_total',[size(Tau_total,1) size(Tau_total,2)]);
% h5write(strcat(filename,'_Tau_total.h5'),'/Tau_total',Tau_total);
% 
% % Save resonance Frequency data in h5 format
% h5create(strcat(filename,'_w_mid_rad_total.h5'),'/w_mid_rad_total',[size(w_mid_rad_total,1) size(w_mid_rad_total,2)]);
% h5write(strcat(filename,'_w_mid_rad_total.h5'),'/w_mid_rad_total',w_mid_rad_total);
% 
% % Save SED fit data in h5 format
% h5create(strcat(filename,'_SED_fit_total.h5'),'/SED_fit_total',[size(SED_fit_total,1), size(SED_fit_total,2),size(SED_fit_total,3)]);
% h5write(strcat(filename,'_SED_fit_total.h5'),'/SED_fit_total',SED_fit_total);
% 
% % Save freqs data in h5 format
% h5create(strcat(filename,'_freqs_total.h5'),'/freqs_total',[size(freqs_total,1), size(freqs_total,2), size(freqs_total,3)]);
% h5write(strcat(filename,'_freqs_total.h5'),'/freqs_total',freqs_total);
% 
% % Save index data in h5 format
% h5create(strcat(filename,'_index_total.h5'),'/index_total',[size(index_total,1), size(index_total,2), size(index_total,3)]);
% h5write(strcat(filename,'_index_total.h5'),'/index_total',index_total);
% 
% plot_name1_total= convertStringsToChars(plot_name1_total);
% % Save plot name1 data 
% save(strcat(filename,'_plot_name1_total.mat'),'plot_name1_total');
% 
% % % Read tau and frequency data if calculation has already been performed and
% % % data have been saved
% % Tau_total= h5read(strcat(filename,'_Tau_total.h5'),'/Tau_total');
% % w_mid_rad_total= h5read(strcat(filename,'_w_mid_rad_total.h5'),'/w_mid_rad_total');
% % SED_fit_total= h5read(strcat(filename,'_SED_fit_total.h5'),'/SED_fit_total');
% % freqs_total= h5read(strcat(filename,'_freqs_total.h5'),'/freqs_total');
% % index_total= h5read(strcat(filename,'_index_total.h5'),'/index_total');
% % load(strcat(filename,'_plot_name1_total.mat'),'plot_name1_total');
% 
% 
% 
% 
% % plot_name2=strcat('Tau vs Freq',' ',material,' ',type,' ',num2str(ntest),',',num2str(fraction),' ','N',num2str(Nunitcells),' ','k',num2str(Nk));
% % close all
% % for i=1:Nk
% % scatter(log10(w_mid_rad_total(i,:)),log10(Tau_total(i,:)));
% % hold on
% % end
% % xlabel('Log(Frequency) (rad/s)')
% % ylabel('Log(Tau) (s)')
% % title('Tau vs Frequency for 2x2x2 RDX')
% % axis([12 15 -12.5 -9.5])
% % saveas(gcf,strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/Tau_vs_Frequency.fig'))
% 
% 
% 
% % %%%%%%%%% For plotting all the SED fits and making video
% % close all
% % opengl('software')% open basic version to avoid memory problem
% % clc
% % for kk=Nk:Nk
% %     k_no= num2str(kk); %ith k point calculation
% % SED_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_SED.dat');
% % data = importdata (SED_file);   %#1 or 2 represents index of k-point
% % frequency_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_frequencies.dat');
% % w0 = importdata (frequency_file);   %#1 or 2 represents index of k-point
% % fprintf('Done importing SED and frequency data \n');
% % 
% % freqs= data(:,1);   % first row is frequency, domain of frequency
% % SED= data(:,2:end); % SED data for all modes, usually 2:end, 5:end when ignoring first 3 acoustic modes
% % % therefore starting 5
% % SED(SED<1e-6)= 1e-6;
% % [f,m]= size(SED);   %no. of discrete freqs and no. of modes
% % 
% % clear F
% % clear G
% % 
% % loops = Neig;
% % F(loops) = struct('cdata',[],'colormap',[]);
% % G(loops) = struct('cdata',[],'colormap',[]);
% % % %  v1 = VideoWriter('silicon_SED_plots.avi');
% % % %  v2 = VideoWriter('silicon_6x6x6_Phonon_Lifetimes.avi');
% % 
% % v1 = VideoWriter(char(plot_name1_total(kk,:)));
% % % v2 = VideoWriter(plot_name2); %for lifetimes video
% % 
% %     v1.FrameRate = 2;
% % %     v2.FrameRate = 2;
% % counter=1;
% % for i=1:loops
% % figure(i);
% %     % when plotting the entire SED
% %     % plot(freqs,log(SED(:,i)),'Color',[0 (1-i/loops) 1.0],'LineWidth',1);
% %     plot(freqs(index_total(kk,1,i):index_total(kk,2,i)),(log(SED(index_total(kk,1,i):index_total(kk,2,i),i))),'Color',[0 (1-i/loops) 1.0],'Linewidth',1);
% %      set(gca,'FontName','Helvetica')
% %      xlabel('Frequency cm^-^1');
% %      ylabel(strcat('Log(SED)'));
% %      %axis([0,800,-15,12])   %for silicon
% % %     axis([0,4500,2,18])   %for RDX
% %      hold on
% %      plot(freqs(index_total(kk,1,i):index_total(kk,2,i)),(SED_fit_total(kk,index_total(kk,1,i):index_total(kk,2,i),i)),'r','Linewidth',2);
% %       legend('Actual SED', 'SED fit');
% %       plot_name=strcat('SED fit Mode-#',num2str(i));
% %       title(plot_name,'Interpreter','none');
% %       F(counter)=getframe(gcf);
% %       close 
% %       
% % %       figure(i);
% % % plot_name=strcat(plot_name2,'-#',num2str(i));
% % % set( gca,'FontName', 'Helvetica' );
% % % loglog(w_mid_rad(1:i),lifetimes(1:i),'o','MarkerSize',6,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.0 0.9 0.9]);
% % % title(plot_name,'Interpreter','none');
% % % axis([1e12,1e15,1e-14,1e-10]);
% % %      xlabel('Frequency (cm^-^1)');
% % %      ylabel('Lifetimes (s)');
% % %       G(counter)=getframe(gcf);
% % %       close 
% %       
% %       counter=counter+1;
% % end
% %       open(v1)
% %       writeVideo(v1,F)
% %       close(v1)
% %       
% % %       open(v2)
% % %       writeVideo(v2,G)
% % %       close(v2)
% % end
%       
% 