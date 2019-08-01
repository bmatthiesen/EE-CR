"aP:wnfunction c = plot_data_reduce(c, sw)
for ii=1:numel(c)
    min_x = floor(min(c{ii})); min_x = min_x(1);
    max_x = ceil(max(c{ii})); max_x = max_x(1);

    step = floor(size(c{ii},1)/length(min_x:sw:max_x));

    if step > 1
        c{ii} = c{ii}(1:step:end,:);
    end
end
end
