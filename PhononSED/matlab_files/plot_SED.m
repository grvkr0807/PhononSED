clc
i= 454;
x=c_light*2*pi*freqs;
y= SED(:,i);
tt=Tau_total(i);
figure
loglog(x, y)
hold on
xx= w_mid_rad_total(i);
line([xx xx],[1 1e8],'Color','red')
xlabel('Omega cm^-1')
ylabel('SED')
title(strcat('RDX 300 K, Mode index ',num2str(i)))

