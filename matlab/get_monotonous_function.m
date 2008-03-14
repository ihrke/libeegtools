function [ fct ] = get_monotonous_function( a, b, step )
% get_monotonous_function return a strictly monotonous function
% a    - a 1x2 matrix (lower point of the line)
% b    - a 1x2 matrix (upper point of the line)
% step - sampling step

x = a(1):step:b(1)-step;
nperiods = 10;
ps = (b(2)-b(1))/2;
noise_strength = ps;

% create a periodic, noisy, positive function and integrate
y = ps*sin(1/((b(1)-a(1))/(nperiods*pi))*x)+(ps+noise_strength);
y = y+noise_strength*rand(1,length(x));
y = cumsum(y);

fct = y/(y(end)/(b(2)-0.0001));
if ~all(diff(fct)>0), disp('Error: function not monotonic'); end;
return;
% plot(x,y);
% hold on;
% plot(x,fct,'r');
% hold off;
% 
% fct = [a(2)];
% x = a(1):step:b(1)-step;
% l = interp1([a(1) b(1)], [a(2) b(2)], a(1):step:b(1), 'linear');
% yrange = b(2)-a(2);
% xrange = b(1)-a(1);
% 
% num_crossings = round(size(x,2)/5);
% p = 0; % plotting? (debugging)
% 
% if p==1
%     close all;
%     figure;
%     plot([a(1) b(1)], [a(2) b(2)]);
%     axis([a(1)-100 b(1)+100 a(2)-100 b(2)+100]);
%     hold on;
% end;
% 
% 
% % crossings is a num_crossings x 2 matrix, containing the points on
% % the input line where it is crossed
% tmp = (rand(num_crossings,1).*xrange)+a(1);
% crossings = zeros(num_crossings,2);
% for c = 1:size(tmp,1)
%     [v i] = min(abs(x-tmp(c)));
%     crossings(c,:) = [x(i) l(i)];
% end;
% [tmp m] = unique(crossings(:,1));
% crossings = sort(crossings(m,:),1);
% crossings = [crossings; b(1) b(2)];
% 
% %fct((crossings(:,1)))=crossings(:,2);
% lborder = a(2); % moving lower border
% crange = diff(crossings(:,2));     % y-range between two crossing-points
% % find a random traversal through the crossing points
% for i=2:size(x,2)
%     if ~isempty(find(crossings(:,1)==x(i)))
%         lborder = crossings(find(crossings(:,1)==x(i)),2);
%         fct(i)=lborder;
%         continue;
%     end;
% 
%     j = find(x(i)<crossings(:,1),1);
%     cy = crossings(j,2);
%     avail = cy-lborder; % available range
%     assert(avail>0, 'no way, dude');
%     if j<=1, pcy = 0; else, pcy=crossings(j-1,2); end;
%     remainder = avail/(cy-pcy);
%     while 1   
%         if remainder*8>1, k = remainder*8; else, k=1; end;
%         r = (rand(1).^k)*avail; % map it to x^n to encourage low prob
%         %  r = (0.5.^(k))*avail;
%         if r<=0, continue; end;
%         if (cy-lborder-r)>0, break; end;
%     end;
%     lborder = lborder+r;
%     fct(i)=lborder;
% end;
%             
% if p==1
%     plot(x, fct, 'r.');
%     plot(fct, x, 'k');
%     plot(crossings(:,1), crossings(:,2), 'ro');
%     axis square;
%     hold off;
% end;
% 
% if(~all(diff(fct)>0)), disp('function not strictly monotonous'); end;
% 
% 
% return;