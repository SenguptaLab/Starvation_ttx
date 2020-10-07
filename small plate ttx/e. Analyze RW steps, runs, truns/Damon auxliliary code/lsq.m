function [p,s]=lsq(FUN,p0,x,y,e,varargin);

if length(varargin)>1
    opts=varargin{2};
else
    opts=optimset;
end
[p,chi2]=fminsearch(@sqerr,p0,opts,FUN,x,y,e);
s.chi2=chi2;
s.df=length(x)-length(p);
P_fit=chisqp(chi2,s.df);
P_fit=2*min(P_fit,1-P_fit);
s.p_fit=P_fit;

%% add in stuff here for the error in p -- hessian, etc.

if length(varargin)
    if varargin{1}>0
        if varargin{1}==2
            figure;
        end
        ha=errorbar(x,y,e);
        for i=1:length(ha)
            set(ha,'color','black');
        end;
        hold on;
        yfit=feval(FUN,p,x);
        plot(x,yfit,'k-');
%         keyboard;
    end
end


function fiterr=sqerr(p,FUN,x,y,e);

yfit=feval(FUN,p,x);
fiterr=sum((y-yfit).^2./e.^2);