function v = cal_speed(time, Pos,bin_size)

dim = size(Pos,2);
bin_range = (1:1:bin_size) - ceil(bin_size/2);
        
        v = nan(length(time),dim);
        
        for d = 1:size(Pos,2)
            for t_iter = 1 : length(time)
                temp_range = t_iter + bin_range;
                temp_range(temp_range<1 | temp_range>length(time)) = [];
                if length(temp_range) < ceil(bin_size/2)+1, continue; end % speed not assigned at both ends of position data
                
                x_range = Pos(temp_range(1):temp_range(end),d);
                x_range(isnan(x_range))=[];
                if isempty(x_range), temp_distance=0;
                else
                    temp_distance =  x_range(end) - x_range(1);
                end
                temp_time = time(temp_range(end)) - time(temp_range(1));
                
                v(t_iter,d) = (temp_distance) / temp_time; % centimeter/sec
            end
        end
end
