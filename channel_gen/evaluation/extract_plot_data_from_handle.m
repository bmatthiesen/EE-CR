"aP:wnfunction v = extract_plot_data_from_handle(h)
v = [get(h, 'XData')', get(h,'YData')'];
v(abs(v(:,1)) == inf, :) = [];
end
