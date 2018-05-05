function res = WENO5(w,flux,dflux,S,dx)
% Lax-Friedrichs Flux Splitting
% a=max(abs(dflux(w))); 
% v=0.5*(flux(w)+a*w); 
% u=circshift(0.5*(flux(w)-a*w),[0,-1]);
% 
% % Left Flux
% hp = lflux(u);
% 
% % Right Flux
% hn = rflux(v);
% 
% % Compute finite volume residual term, df/dx.
% res = -((hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]))/dx - S(w));

[pfx,nfx]=fluxsplitting(w,flux,dflux,'LF');
lfx = rflux(nfx,[0 -1]);
rfx = rflux(pfx,[0 +1]);

res = lfx - circshift(lfx,[0,-1]) + rfx - circshift(rfx,[0, 1]);

res = -res/dx + S(w);

end

function [vp,vn] = fluxsplitting(u,f,df,strategy)
% split the flux into right-going and left-going
% OUTPUT:
%   * vp: positive flux, v^{+}, which corresponds to f_{i+1/2}^{-}
%   * vn: negative flux  v^{-}, which corresponds to f_{i+1/2}^{+}

switch strategy
    case 'Upwind' % Godunov/scalar fluxsplit (non-conservative)
        vp = f((u + abs(u))./2); %flux^{+}
        vn = f((u - abs(u))./2); %flux^{-}
    case 'LLF' % Local Lax-Friedrichs
        v = f(u); alpha = abs(df(u));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    case 'LF' % (Global) Lax-Friedrichs
        v = f(u); alpha = max(max(abs(df(u))));
        vp = 0.5.*(v + alpha.*u); %flux^{+}
        vn = 0.5.*(v - alpha.*u); %flux^{-}
    otherwise
        error('only case not available')
end

end
