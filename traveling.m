%FUNCTION TRAVELING: Computes travel time between bent points.
function [ta,tra]=traveling(xtemp,ytemp,ztemp,gridD,v,pvel)
    
    distlim=5*(gridD(1,2)-gridD(1,1));
    n = length(xtemp);
    xd=xtemp(2)-xtemp(1);
    yd=ytemp(2)-ytemp(1);
    zd=ztemp(2)-ztemp(1);
    ds=sqrt(xd^2+yd^2+zd^2);
    if ds <= distlim
        tra=ds*((1/v(1)+1/v(2))/2);
        ta(2)=tra;
        for i=3:n
            i1=i-1;
            xd=xtemp(i)-xtemp(i1);
            yd=ytemp(i)-ytemp(i1);
            zd=ztemp(i)-ztemp(i1);
            ds=sqrt(xd^2+yd^2+zd^2);
            tra=tra+ds/((v(i)+v(i1))/2);
            ta(i)=tra;
        end
    else
        tra=0;
        ta=zeros(1,n);
        for i=2:n
            i1=i-1;
            xd=xtemp(i)-xtemp(i1);
            yd=ytemp(i)-ytemp(i1);
            zd=ztemp(i)-ztemp(i1);
            ds=sqrt(xd^2+yd^2+zd^2);
            xcos=xd/ds;
            ycos=yd/ds;
            zcos=zd/ds;
            is=fix(ds/distlim);
            ds=ds/is;
            for j=1:is
                xx=xtemp(i1)+(j-0.5)*ds*xcos;
                yy=ytemp(i1)+(j-0.5)*ds*ycos;
                zz=ztemp(i1)+(j-0.5)*ds*zcos;
                vp=velocity(xx,yy,zz,gridD,pvel);
                tra=tra+ds/vp;
            end
            ta(i)=tra;
        end
    end
end