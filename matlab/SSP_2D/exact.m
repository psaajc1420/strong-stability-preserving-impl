function uExact = exact(xx,yy,t,type)

    switch type
        case 'gaussian'
            xc = 0;
            yc = 0;
            sigma = 1;
            A = 1/(2*sqrt(2));
            uExact  = A*exp(-(((xx-xc-t).^2 + (yy-yc-t).^2)./(2*sigma^2)));
        case 'sines'
            uExact = sin(2*pi*(xx-t)).*sin(2*pi*(yy-t));  
        case 'sinxy'
            uExact = sin(pi*(xx + yy - 2*t));            
    end

end