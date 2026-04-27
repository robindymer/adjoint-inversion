    function receivers = placeReceivers(boundaries, spacing)
        receivers = struct;
        receivers.x = {};
        receivers.blockIds = [];
        
        x_o_l = boundaries.outer.x(1);
        x_o_r = boundaries.outer.x(2);
        y_o_l = boundaries.outer.y(1);
        y_o_r = boundaries.outer.y(2);
        
        x_i_l = boundaries.inner.x(1);
        x_i_r = boundaries.inner.x(2);
        y_i_l = boundaries.inner.y(1);
        y_i_r = boundaries.inner.y(2);

        
        dr = spacing;

        x = x_o_l:dr:x_o_r;
        x = x ((x  <= x_i_l) | (x  >= x_i_r));
        y = y_o_l:dr:y_o_r;
        for i = 1:length(x)
            for j = 1:length(y)
                if y(j) > 0
                    block_id = 2;
                else
                    block_id = 1;
                end
                receivers.x{end+1} = [x(i), y(j)];
                receivers.blockIds(1,end+1) = block_id;
            end
        end
        
        x = x_i_l+dr:dr:x_i_r-dr;
        %x = x_i_l:dr:x_i_r;
        y = y_o_l:dr:y_o_r;
        y = y((y <= y_i_l) | (y >= y_i_r));
        for i = 1:length(x)
            for j = 1:length(y)
                if y(j) > 0
                    block_id = 2;
                else
                    block_id = 1;
                end
                receivers.x{end+1} = [x(i), y(j)];
                receivers.blockIds(1,end+1) = block_id;
            end
        end
        fprintf('Number of receivers: %d\n', numel(receivers.x));
        fprintf('Receiver spacing: %d\n', dr);
    end