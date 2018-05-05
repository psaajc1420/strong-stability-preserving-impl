function dF = residual(q,f,df,g,dg,dx,dy,S,strategy)

% Along x
[pfx,nfx]=fluxsplitting(q,f,df,strategy);
rxfhat = rflux(pfx,[0  1]);
lxfhat = rflux(nfx,[0 -1]);

dqx =lxfhat - circshift(lxfhat,[0, -1]) + rxfhat - circshift(rxfhat,[0 1]);

% Along y
[pfy,nfy]=fluxsplitting(q,g,dg,strategy);
ryfhat = rflux(pfy,[ 1 0]);
lyfhat = rflux(nfy,[-1 0]);

dqy = lyfhat - circshift(lyfhat,[-1,0]) + ryfhat- circshift(ryfhat,[1,0]);

% The Residual
dF = -dqx/dx - dqy/dy + S(q);