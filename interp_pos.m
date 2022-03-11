function pos = interp_pos(pos,range,c)

temp = pos(range);
x = find(~isnan(temp));
y = temp(x);
xq = [1:length(temp)]';

yq = round(interp1(x,y,xq,'linear','extrap'));

pos(range) = yq;
pos(~range) = c;

end
%
