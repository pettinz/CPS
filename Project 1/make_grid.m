function make_grid(x, y, xs, ys)
%MAKE_GRID Summary of this function goes here
%   Detailed explanation goes here
    plot(xs,ys,'.r'), grid on
    xticks(x)
    yticks(y)
    xlim([x(1) x(end)])
    ylim([y(1) y(end)])

end

