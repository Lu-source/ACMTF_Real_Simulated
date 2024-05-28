function [rmin, rmax] = bar_wrange(A,B,l,x1,x2, lims)
% This function plots bar plots based on the averages of the columns of
% input matrices with lines showing the min and max value of each colum.
% Corresponding columns of input matrices are plotted together.

m    = [mean(A)' mean(B)'];
mmin = [min(A)' min(B)'];
mmax = [max(A)' max(B)'];
rmin = m-mmin;
rmax = mmax-m;
%f1= figure;
hb=bar(m, 'grouped'); hold on;
% x1 = 0.85:1:3.85;
% x2 = 1.15:1:4.15;
hb(2).FaceColor = [0 .7 .3];


errorbar(x1,m(:,1),rmin(:,1),rmax(:,1),'+r','Linewidth',4); hold on;
errorbar(x2,m(:,2),rmin(:,2),rmax(:,2),'+r','Linewidth',4); hold on;

axis(lims);
hold on;

if nargin > 2
    legend(l);
end
