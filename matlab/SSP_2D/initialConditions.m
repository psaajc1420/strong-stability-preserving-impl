function u0 =  initialConditions(xx,yy,type)

    switch type
        case 'gaussian'
            xc = 0;
            yc = 0;
            sigma = 1;
            A = 1/(2*sqrt(2));
            u0  = A*exp(-(((xx-xc).^2 + (yy-yc).^2)./(2*sigma^2)));
        case 'sines'
            u0 = sin(2*pi*xx).*sin(2*pi*yy);  
        case 'sinxy'
            u0 = sin(pi*(xx + yy));
              
    end

end