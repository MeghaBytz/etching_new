    function y = h_side(x, delx)
        ep = 1.5*delx(1);
        y = zeros(size(x));
        y(x<(-1*ep))= 0;
        y(x>ep) = 1;
        z = 0.5 + x./(2*ep) + (0.5/pi)*sin((pi*x)./(ep));
        y((x>-1*ep)&(x<ep)) = z((x>-1*ep)&(x<ep));%0.5 + x./(2*ep) + (0.5/pi)*sin((pi*x)./(ep));
    
        %{
        if (x<(-1*eps))
            y = 0;
        elseif (x>(eps))
            y = 1;
        else
            y = 0.5 + x/(2*eps) + (0.5/pi)*sin((pi*x)/(eps));
        end
        %}
    end