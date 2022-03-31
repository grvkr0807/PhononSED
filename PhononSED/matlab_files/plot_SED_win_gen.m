close all
opengl('software')% open basic version to avoid memory problem
clc
for kk=Nk:Nk
    k_no= num2str(kk); %ith k point calculation
SED_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_SED.dat');
data = importdata (SED_file);   %#1 or 2 represents index of k-point
frequency_file= strcat('/home/gkumar1/Phonon_Lifetimes/PhononSED/',material,'_N',num2str(Nunitcells),'_k',num2str(Nk),'/',material,'_N',num2str(Nunitcells),'_k',k_no,'_frequencies.dat');
w0 = importdata (frequency_file);   %#1 or 2 represents index of k-point
fprintf('Done importing SED and frequency data \n');

freqs= data(:,1);   % first row is frequency, domain of frequency
SED= data(:,2:end); % SED data for all modes, usually 2:end, 5:end when ignoring first 3 acoustic modes
% therefore starting 5
SED(SED<1e-6)= 1e-6;
[f,m]= size(SED);   %no. of discrete freqs and no. of modes

clear F
% clear G

loops = Neig;
F(loops) = struct('cdata',[],'colormap',[]);
% G(loops) = struct('cdata',[],'colormap',[]);
% %  v1 = VideoWriter('silicon_SED_plots.avi');
% %  v2 = VideoWriter('silicon_6x6x6_Phonon_Lifetimes.avi');

v1 = VideoWriter(char(plot_name1_total(kk,:)));
% v2 = VideoWriter(plot_name2); %for lifetimes video

    v1.FrameRate = 2;
%     v2.FrameRate = 2;
counter=1;
for i=1:loops
figure(i);
    % when plotting the entire SED
    % plot(freqs,log(SED(:,i)),'Color',[0 (1-i/loops) 1.0],'LineWidth',1);
    plot(freqs(index_total(kk,1,i):index_total(kk,2,i)),(log(SED(index_total(kk,1,i):index_total(kk,2,i),i))),'Linewidth',1);
     set(gca,'FontName','Helvetica')
     xlabel('Frequency cm^-^1');
     ylabel(strcat('Log(SED)'));
     %axis([0,800,-15,12])   %for silicon
%     axis([0,4500,2,18])   %for RDX
     hold on
     plot(freqs(index_total(kk,1,i):index_total(kk,2,i)),((SED_fit_total(kk,index_total(kk,1,i):index_total(kk,2,i),i))),'r','Linewidth',2);
      legend('Actual SED', 'SED fit');
      plot_name=strcat('SED fit Mode-#',num2str(i));
      title(plot_name,'Interpreter','none');
      F(counter)=getframe(gcf);
      close 
      
%       figure(i);
% plot_name=strcat(plot_name2,'-#',num2str(i));
% set( gca,'FontName', 'Helvetica' );
% loglog(w_mid_rad(1:i),lifetimes(1:i),'o','MarkerSize',6,'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.0 0.9 0.9]);
% title(plot_name,'Interpreter','none');
% axis([1e12,1e15,1e-14,1e-10]);
%      xlabel('Frequency (cm^-^1)');
%      ylabel('Lifetimes (s)');
%       G(counter)=getframe(gcf);
%       close 
      
      counter=counter+1;
end
      open(v1)
      writeVideo(v1,F)
      close(v1)
      
%       open(v2)
%       writeVideo(v2,G)
%       close(v2)
end
      