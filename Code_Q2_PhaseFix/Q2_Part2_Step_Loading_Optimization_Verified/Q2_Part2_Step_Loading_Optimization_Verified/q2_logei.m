function logEI = q2_logei(bestTime,mu,sd)
%Q2_LOGEI Stable analytic log expected improvement for time minimization.
% Ament et al. (NeurIPS 2023), Eq. (8)-(9). This is the analytic LogEI,
% not log(max(naively-computed-EI,realmin)).
assert(isscalar(bestTime) && isfinite(bestTime));
assert(isequal(size(mu),size(sd)) && all(isfinite(mu(:))) && ...
    all(isfinite(sd(:))) && all(sd(:)>=0));
sd=max(double(sd),1e-12);
z=(bestTime-double(mu))./sd;
logH=zeros(size(z));
ordinary=z>-1;
if any(ordinary(:))
    zo=z(ordinary);
    phi=exp(-0.5*zo.^2)/sqrt(2*pi);
    Phi=0.5*erfc(-zo/sqrt(2));
    logH(ordinary)=log(phi+zo.*Phi);
end
moderate=~ordinary & z>-1/sqrt(eps);
if any(moderate(:))
    zm=z(moderate);
    q=log(erfcx(-zm/sqrt(2)).*abs(zm))+0.5*log(pi/2);
    q=min(q,-eps);
    logH(moderate)=-0.5*zm.^2-0.5*log(2*pi)+log1mexp(q);
end
extreme=~ordinary & ~moderate;
if any(extreme(:))
    ze=z(extreme);
    logH(extreme)=-0.5*ze.^2-0.5*log(2*pi)-2*log(abs(ze));
end
logEI=log(sd)+logH;
end

function out=log1mexp(x)
assert(all(x(:)<0));
out=zeros(size(x));
lo=x<-log(2);
out(lo)=log1p(-exp(x(lo)));
out(~lo)=log(-expm1(x(~lo)));
end
