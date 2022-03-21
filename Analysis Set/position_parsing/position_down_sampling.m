% position data down sampling function 
% convert 60Hz position data to 30Hz position data
% position_down_sampling(pos_file, epochST, epochED)
% ex) posfile = 'VT1.nvt', epochST and epochED in 1/1000000 second unit
% Code by Hyunwoo Lee, 2014-Dec-16

function [t, x, y, a] = position_down_sampling(t, x, y, a)

% [t x y a] = Nlx2MatVT(pos_file, [1 1 1 1 0 0], 0, 1, 0);

if abs(1 - (t(61) - t(1))/1000000) > 0.1
    disp([num2str((t(61)-t(1))/1000000) ', this is not 60Hz']);
    return;
end
    
    d = diff(t(1:10));
    if abs(d(1) - d(2)) <= 1
        even_line = 2;
    elseif abs(d(2) - d(3)) <= 1
        even_line = 1;
    end
    
    t(even_line : 2 : end) = [];
    x(even_line : 2 : end) = [];
    y(even_line : 2 : end) = [];
    a(even_line : 2 : end) = [];

% out_range = find(epochST > t | t > epochED);
% t(out_range) = [];
% x(out_range) = [];
% y(out_range) = [];
% a(out_range) = [];

disp('Position data down sampling completed.');

end