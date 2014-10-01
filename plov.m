%plov: plotting x,y function and if mov==1 plotting video of trajectory
plot(x,y)
if mov
    for n=1:10:numel(tt)
        plot(x,y,x(n),y(n),'*')
        pause(0.001)
    end
end