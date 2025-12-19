function [V1,V2,a,Nrev,index]=Lambert_aro3090(R1,R2,T,mus,Nrev1,Hd)

%
%OUTPUTS
%V1 = velocity vector at t1 (columns)
%V2 = velocity vectot at t2 (columns)
%a = semi-major axis of orbit
%N = number of revolutions around central body, see N input
%index = list of input indeces corresponding to outputs
%INPUTS
%R1 = 3d position at t1 (columns)
%R2 = 3d position at t2 (columns)
%tof = time of flight from R1 to R2
%GM = GM of central body
%N = desired number of revolutions about central body (N > 0 for long period solution, N < 0 for short)
%    if N is imaginary, then all possible multiple revolution arcs up to abs(N) will be found, or
%    if specified N value is impossible for corresponding R1,R2, and tof, then 0 rev case is calculated
%    e.g. N=2i automatically looks for N=0,-1,1,-2,2 solutions & outputs are tracked via index
%H = general direction of anglular momentum, default = [0;0;1] ([0;0;-1] is retrograde)

%Multiple inputs/outputs expand the 2nd dim.
%so V1(:,n) corresponds to the nth input column R1(:,n), etc.

%Gooding, R.H., "A Procedure for the Solution of Lambert's Orbital Boundary-Value Problem," Celestial Mechanics and Dynamical Astronomy, v48, 145-165, 1990.
%Lancaster, E.R., and Blanchard, R.C., "A Unified Form of Lambert's Theorem," Technical Report NASA TN D-5368, NASA Technical Note, 1969. 

if numel(T)==1;T=repmat(T,[1 size(R1,2)]);end
if numel(mus)==1;mus=repmat(mus,size(T));end
if size(R1,2)==1;R1=repmat(R1,size(T));end
if size(R2,2)==1;R2=repmat(R2,size(T));end
if ~exist('Nrev1')|isempty(Nrev1);Nrev1=0;end
if any(imag(Nrev1)<0);warning('negative imaginary Nrev');ii=imag(Nrev1)<0;Nrev1(ii)=-Nrev1(ii);end
if numel(Nrev1)==1;Nrev1=Nrev1*ones(size(T));end
if ~exist('Hd')|isempty(Hd);Hd=[0;0;1];end
if size(Hd,1)==1;Hd(3,:)=Hd;Hd(1,:)=0;end
if numel(Hd)==3;Hd=Hd*ones(size(T));end

im=find(T<0);if any(im);Rm=R1(:,im);R1(:,im)=R2(:,im);R2(:,im)=Rm;T(im)=-T(im);end

%index
nn=numel(T);nm=2e6/min(max(abs(Nrev1(:))),1e1);%max number to run in main code at a time (memory)
if nn<nm
[V1,V2,a,Nrev,index]=LambertP2(R1,R2,T,mus,Nrev1,Hd);
else
in=round(nn/ceil(nn/nm));
V1=[];V2=[];a=[];Nrev=[];index=[];res=[];jjj=0;
while ~isempty(T)
    jj=min(numel(T),in);
    [V1a,V2a,aa,Nreva,indexa]=LambertP2(R1(:,1:jj),R2(:,1:jj),T(1:jj),mus(1:jj),Nrev1(:,1:jj),Hd(:,1:jj));
    R1(:,1:jj)=[];R2(:,1:jj)=[];Hd(:,1:jj)=[];T(1:jj)=[];mus(1:jj)=[];Nrev1(:,1:jj)=[];
    if nargout>2;a=[a aa];Nrev=[Nrev Nreva];jj=jjj+[1:jj];jjj=jjj+numel(jj);index=[index jj(indexa)];end;clear aa Nreva jj indexa
    V1=[V1 V1a];V2=[V2 V2a];clear V1a V2a
end
end
if any(im);im=ismember(index,im);Vm=V1(:,im);V1(:,im)=V2(:,im);V2(:,im)=Vm;end
return

function [V1,V2,a,Nrev,index]=LambertP2(R1,R2,T,mus,Nrev,Hd)
%see if nest helps mem
%Nrev=Nrev1;
nn=numel(T);
[R1 r1]=unit(R1);[R2 r2]=unit(R2);
q=sum(R1.*R2);ii=q>1;if any(ii);q(ii)=1+imag(q(ii))*1i;end;ii=q<-1;if any(ii);q(ii)=-1+imag(q(ii))*1i;end
%sq=sqrt(sum((R1-R2).^2).*(1+q)/2);%small angle
%maybe useful for small angles: tan(q/2)=mag(R1-R2)/mag(R1+R2),sq=mag(R1-R2)*mag(R1+R2)/2;
q=acos(q);
[H h]=unit(cross(R1,R2));
aa=find(h<1e-6);H(:,aa)=unit(Hd(:,aa)-repmat(sum(R1(:,aa).*Hd(:,aa)),3,1).*R1(:,aa));
aa=find(sum(H.*Hd)<0);q(aa)=2*pi-q(aa);H(:,aa)=-H(:,aa);
r1r2q=4*r1.*r2.*sin(q/2).^2;
c=sqrt((r1-r2).^2+r1r2q);
s=(r1+r2+c)/2;
mus=sqrt(mus.*s/2);
cs=c./s;
cq=cos(q/2);q=sqrt(r1.*r2).*cq./s;
rho=(r1-r2)./c;rho(c==0)=0;
iz=T==0;T=4*mus.*T./s.^2;
if nn>1e6;clear Hd c aa;end

[x,Nrev,index]=xlamb(T,q,Nrev,cs);
if nn~=numel(index);q=q(index);cs=cs(index);r1r2q=r1r2q(index);mus=mus(index);rho=rho(index);r1=r1(index);r2=r2(index);s=s(index);R1=R1(:,index);R2=R2(:,index);H=H(:,index);cq=cq(index);end
[spot qzm qz zq]=Tx(x,q,Nrev,cs,1);
%v2=mus.*zq.*sqrt(r1r2q)./c;
u1=mus.*(qzm-qz.*rho)./r1;
v2=mus.*zq.*sqrt(1-rho.^2);
v1=v2./r1;
u2=-mus.*(qzm+qz.*rho)./r2;
v2=v2./r2;
a=[-s/2./(x.^2-1);cq];
%save lamdat x

if any(iz);iz=iz(index);ii=iz&rho==0;u1(ii)=0;u2(ii)=0;v1(ii)=0;v2(ii)=0;ii=iz&rho~=0;u1(ii)=nan;u2(ii)=nan;v1(ii)=nan;v2(ii)=nan;end
if nn>1e6;clear T cs mus q qz qzm r1 r1r2q r2 rho zq s x;end
V1=[1;1;1]*u1.*R1+[1;1;1]*v1.*cross(H,R1);
if nn>1e6;clear u1 v1 R1;end
V2=[1;1;1]*u2.*R2+[1;1;1]*v2.*cross(H,R2);

return

function [x,Nrev,index,res]=xlamb(T,q,Nrev,cs)
%~N&Nrev<0 similar
%N=abs(Nrev);
%x=zeros(size(T));
index=1:numel(T);nn=numel(index);
c0=1.7;c1=.5;c2=.03;c3=.15;c41=1;c42=.24;
%thr2=atan22(cs,2*q)/pi;
%T0=Tx(0*x,q,abs(Nrev),cs);
%DT=T-T0;
%intitial guess
%multi-rev
iN=imag(Nrev);Nrev=real(Nrev);
xm=[];Tm=[];d2T=[];
aa=find(abs(Nrev)>0);if ~isempty(aa)
    [xm Tm d2T]=minT(q(aa),Nrev(aa),cs(aa),T(aa));
    if 1%delete infeasible soln's
        bb=find(T(aa)+1e-3<Tm);ab=aa(bb);
        index(ab)=[];T(ab)=[];q(ab)=[];Nrev(ab)=[];cs(ab)=[];iN(ab)=[];
        xm(bb)=[];Tm(bb)=[];d2T(bb)=[];
    else
        ii=T(aa)<Tm;Nrev(aa(ii))=0;Tm(ii)=[];xm(ii)=[];d2T(ii)=[];
    end
    %tack on complement to minN if range of Nrev
    aa=find(Nrev);bb=find(iN(aa));ab=aa(bb);
    Nrev=[Nrev -Nrev(ab)];
    index=[index index(ab)];T=[T T(ab)];q=[q q(ab)];cs=[cs cs(ab)];
    xm=[xm xm(bb)];Tm=[Tm Tm(bb)];d2T=[d2T d2T(bb)];
end
%Expand to max revs
aa=find(iN);Na=abs(Nrev);
while ~isempty(aa)
    aa=aa(iN(aa)>Na(aa));Na=Na+1;
    [xma Tma d2Ta]=minT(q(aa),Na(aa),cs(aa),T(aa));
    bb=find(T(aa)>=Tma);aa=aa(bb);if isempty(aa);break;end
    oaa=Na(aa).'*[-1 1];Nrev=[Nrev oaa(:).'];
    xma=xma(bb);xm=[xm xma xma];
    Tma=Tma(bb);Tm=[Tm Tma Tma];
    d2Ta=d2Ta(bb);d2T=[d2T d2Ta d2Ta];
    q=[q q(aa) q(aa)];
    cs=[cs cs(aa) cs(aa)];
    T=[T T(aa) T(aa)];
    index=[index index(aa) index(aa)];
end

if isempty(Nrev);x=[];res=[];return;end%nada

%feasible multi-rev
aa=find(abs(Nrev)>0);if ~isempty(aa)
    
    Na=abs(Nrev(aa));
    DTm=max(T(aa)-Tm,0);
    nn=numel(index);

    bb=find(Nrev(aa)>0);if ~isempty(bb)%x>xm
        DTmb=DTm(bb);xmb=xm(bb);Nb=Na(bb);d2Tb=d2T(bb);
        if any(abs(d2Tb)<eps);d2Tb(abs(d2Tb)<eps)=6*Nb(abs(d2Tb)<eps)*pi;end
        xb=sqrt(DTmb./(d2Tb/2+DTmb./(1-xmb).^2));
        W=xmb+xb;
        W=W*4./(4+DTmb)+(1-W).^2;
%        xb=xb.*(1-(1+Nb+c41.*(thr2(aa(bb))-.5))./(1+c3.*Nb).*xb.*(c1.*W+c2.*xb.*sqrt(W)))+xmb;
        xb=xb.*(1-(1+Nb+c41.*(atan22(cs(aa(bb)),2*q(aa(bb)))/pi-.5))./(1+c3.*Nb).*xb.*(c1.*W+c2.*xb.*sqrt(W)))+xmb;
        if any(xb>1);disp('X>1');end
        x(aa(bb))=xb;
    end
    bb=find(Nrev(aa)<0);if ~isempty(bb)%x<xm
        
        ab=aa(bb);xb=zeros(size(ab));T0b=Tx(xb,q(ab),Nrev(ab),cs(ab));DTb=T(ab)-T0b;
        cc=find(DTb>0);if ~isempty(cc)%if DT>0
            bc=bb(cc);DTc=DTb(cc);thr2c=atan22(cs(aa(bc)),2*q(aa(bc)))/pi;Nc=Na(bc);
            xc=-DTc./(DTc+4);
            W=xc+c0*sqrt(2*(1-thr2c));
            dd=find(W<0);if ~isempty(dd)
                xd=xc(dd);DTd=DTc(dd);
                xc(dd)=xd-sqrt(rt8(-W(dd))).*(xd+sqrt(DTd./(DTd+1.5*T0b(cc(dd)))));
            end
            W=4./(4+DTc);
            xb(cc)=xc.*(1+(1+Nc+c42.*(thr2c-.5))./(1+c3*Nc).*xc.*(c1*W-c2*xc.*sqrt(W)));
        end
        cc=find(DTb<=0);if ~isempty(cc);
            bc=bb(cc);DTmc=DTm(bc);xmc=xm(bc);d2T2c=d2T(bc)/2;
            xb(cc)=xmc-sqrt(DTmc./(d2T2c-DTmc.*(d2T2c./(T0b(cc)-Tm(bc))-1./xmc.^2)));
        end
        if any(xb<-1);disp('X<1');end
        x(ab)=xb;
    end
end
if nn>1e6;clear DTm DTb DTc DTd DTmb DTmc Na Nb Nc T0b Tm W aa bb ab bc cc dd d2T d2Tb d2T2c thr2c xb xc xd xm xmb xmc;end
%0 rev
aa=find(Nrev==0);if ~isempty(aa)
    xa=zeros(size(aa));T0a=Tx(xa,q(aa),Nrev(aa),cs(aa));DTa=T(aa)-T0a;
    bb=find(DTa>0);if ~isempty(bb)%if DT>0
        DTb=DTa(bb);
        xb=-DTb./(DTb+4);
        W=xb+c0*sqrt(2*(1-atan22(cs(aa(bb)),2*q(aa(bb)))/pi));
        cc=find(W<0);if ~isempty(cc)
            xc=xb(cc);DTc=DTb(cc);
            xb(cc)=xc-sqrt(rt8(-W(cc))).*(xc+sqrt(DTc/(DTc+1.5*T0a(bb(cc)))));
        end
        W=4./(4+DTb);
        xa(bb)=xb.*(1+xb.*(c1*W-c2*xb.*sqrt(W)));
    end
    bb=find(DTa<=0);if ~isempty(bb)
        xa(bb)=T0a(bb).*DTa(bb)./(-4*T(aa(bb)));
    end
    x(aa)=xa;
end
%if numel(x)~=nn;x(nn)=0;end
I=0;
if nn>1e6;clear DTa DTb DTc T0a W aa bb xa xb xc;end
while I<3
    I=I+1;
    [Tc dT d2T]=Tx(x,q,Nrev,cs);
    dT=(T-Tc).*dT./(dT.^2+(T-Tc).*d2T/2);
    x=x+dT;
end
res=abs([T-Tc;dT]);
return
        
function [T,dT,d2T,d3T]=Tx(x,q,N,cs,flag)
nn=numel(x);
if nargin<5;flag=0;end
sw=4e-1;
T=0;
if nargout>1;dT=T;d2T=T;end
if nargout>3;d3T=T;end
%u=(1-x).*(1+x);
aa=(find(N>0|x<0|abs((1-x).*(1+x))>sw));
if flag;aa=1:numel(x);end
if ~isempty(aa)
    xa=x(aa);
    csa=cs(aa);
    ua=(1-xa).*(1+xa);
	z=sqrt(csa+q(aa).^2.*xa.^2);
    qx=q(aa).*xa;
    qz=q(aa).*z;
    term=csa.*(q(aa).^2.*ua-xa.^2);
    if ~flag
        a=z-qx;b=qz-xa;
        bb=find(qx>0);if ~isempty(bb);
            a(bb)=csa(bb)./(z(bb)+qx(bb));
            b(bb)=term(bb)./(qz(bb)+xa(bb));
        end
    else
        b=qz-xa;a2=z+qx;b2=qz+xa;
%these lines don't do *
%        bb=find(qx<0);if ~isempty(bb);a2(bb)=csa(bb)./(z(bb)-qx(bb));b2(bb)=term(bb)./b(bb);end
%        bb=find(qx>0);if ~isempty(bb);b(bb)=term(bb)./b2(bb);end
        dT=b;d2T=b2;d3T=a2;return
    end

    g=zeros(size(aa));
    
    bb=find(qx.*ua>0);if ~isempty(bb);g(bb)=xa(bb).*z(bb)+q(aa(bb)).*ua(bb);end
%    bb=find(qx.*ua<=0);if ~isempty(bb);g(bb)=-term(bb)./csa(bb)./(xa(bb).*z(bb)-q(aa(bb)).*ua(bb));end
    bb=find(qx.*ua<=0);if ~isempty(bb);g(bb)=-(q(aa(bb)).^2.*ua(bb)-xa(bb).^2)./(xa(bb).*z(bb)-q(aa(bb)).*ua(bb));end
    f=a.*sqrt(sqrt(ua.^2));Ta=zeros(size(f));
    bb=find(xa<=1);if ~isempty(bb);Ta(bb)=(abs(N(aa(bb))))*pi+atan(f(bb)./g(bb))+pi*(g(bb)<0);end
    bb=find(xa>1&f>sw);if ~isempty(bb);Ta(bb)=log(f(bb)+g(bb));end
    bb=find(xa>1&f<=sw);if ~isempty(bb);        
        fg12=f(bb)./(g(bb)+1);
        term=2*fg12;
        fg12=fg12.^2;
        Tt=term;
        Td=1;
        while max(abs(term/Td))>3e-14
        Td=Td+2;
        term=term.*fg12;
        Tt=Tt+term/Td;
        end
        Ta(bb)=Tt;
    end
    Ta=2*(Ta./sqrt(sqrt(ua.^2))+b)./ua;T(aa)=Ta;
    if nargout>1
        z=max(z,eps);
        qz=q(aa)./z;qz2=qz.^2;qz=qz.*qz2;
        dTa=(3*xa.*Ta-4*(a+qx.*csa)./z)./ua;dT(aa)=dTa;
        if nn>1e6;clear z qx term a b f g;end
        if nargout>2;d2T(aa)=(3*Ta+5*xa.*dTa+4*qz.*csa)./ua;end
        if nargout>3;d3T(aa)=(8*dTa+7*xa.*d2T(aa)-12*qz.*qz2.*xa.*csa)./ua;end
    end
end
if nn>1e6;clear ua csa  xa dTa z qx term a b f g Ta;end
aa=find(~(N>0|x<0|abs((1-x).*(1+x))>sw));if ~isempty(aa)&~flag
    qa=q(aa);xa=x(aa);csa=cs(aa);q2a=qa.^2;x2a=xa.^2;ua=(1-xa).*(1+xa);%u(aa);
    Tqs=zeros(size(aa));
    if nargout>1;dTa=Tqs;d2Ta=Tqs;end
    if nargout>3;d3Ta=Tqs;end
    u0I=1;u1I=1;u2I=1;u3I=1;
    term=4;
    Tq=qa.*csa;
    I=0;
    bb=find(qa<.5);if ~isempty(bb);Tqs(bb)=1-qa(bb).*q2a(bb);end
    bb=find(qa>=.5);if ~isempty(bb);Tqs(bb)=(1./(1+qa(bb))+qa(bb)).*csa(bb);end
    ttmo=term/3;
    Ta=ttmo*Tqs;
    DT=1;
    while max(abs(DT))>3e-14|I<nargout
        I=I+1;P=I;
        u0I=u0I.*ua;
        if nargout>1&I>1;u1I=u1I.*ua;end
        if nargout>2&I>2;u2I=u2I.*ua;end
        if nargout>3&I>3;u3I=u3I.*ua;end
        term=term*(P-.5)/P;
        Tq=Tq.*q2a;
        Tqs=Tqs+Tq;
        tterm=term/(2*P+3);
        tqterm=tterm.*Tqs;
        DT=u0I.*((1.5*P+.25)*tqterm./(P.^2-.25)-ttmo*Tq);
        Ta=Ta-DT;
        ttmo=tterm;
        tqterm=tqterm*P;
        if nargout>1;dTa=dTa+tqterm.*u1I;end
        if nargout>2;d2Ta=d2Ta+tqterm.*u2I*(P-1);end
        if nargout>3;d3Ta=d3Ta+tqterm.*u3I*(P-1)*(P-2);end
    end
    if nargout>3;d3T(aa)=8.*xa.*(1.5.*d2Ta-x2a.*d3Ta);end
    if nargout>2;d2T(aa)=2*(2*x2a.*d2Ta-dTa);end
    if nargout>1;dT(aa)=-2*xa.*dTa;end
    T(aa)=Ta./x2a;
end
return

function [xm_ Tm_ d2T_]=minT(q,N,cs,T)
persistent Tm_data Tm_fit
aN=abs(N);Nx=max(aN);Tm_=2*pi*aN;xm_=Tm_;d2T_=Tm_;index=1:numel(q);
if nargin>3
ii=T>=Tm_;if ~any(ii);return;elseif sum(ii)/numel(ii)<.7;q=q(ii);N=N(ii);cs=cs(ii);T=T(ii);index=index(ii);aN=abs(N);end

if Nx>numel(Tm_fit)
for N_=1:Nx;if size(Tm_data,2)>=N_&&any(Tm_data(:,N_));continue;end
    q_=[-1 0];[xm Tm d2T]=minT(q_,[N_ N_],1-q_.^2);
    Tm_data(:,N_)=[Tm(2);fminbnd(@(q) (2*acos(q)+2*q.*sqrt(1-q.^2)+diff(Tm)-pi).^2,-1,.98,optimset('tolx',1e-14));0];
%     q_=q;ii=q<0;q_(ii)=-Tm_data(2,aN(ii)).*q(ii);T_=2*(acos(q_)+q_.*sqrt(1-q_.^2))+Tm_data(1,aN)-pi;Tm_(index)=T_;
%     ii=T>=T_;if ~any(ii);return;elseif any(~ii);q=q(ii);N=N(ii);cs=cs(ii);T=T(ii);index=index(ii);aN=abs(N);end

    %fit cubic xm/Tm/d2T + derivative wrt q
    qx=.99;q_=[];q2_=[-qx -.9 -.7 -.4 0 .4 .7 .9 qx];
    for mesh_count=1:9;if isempty(q2_);break;end
    q_=unique([q_ q2_]);
    s_=sin(acos(q_)).^2;[xm Tm d2T]=minT(q_,repmat(N_,size(q_)),s_);
    u=(1-xm).*(1+xm);y=sqrt(u);z=sqrt(s_+q_.^2.*xm.^2);qz=q_./z;f=y.*(z-q_.*xm);g=z.*xm+q_.*u;
    %xz=1./sqrt(s_./xm.^2+q_.^2);%x/z breaks down at q=1,xm=0
    zq=-qz.*u;fq=y.*(zq-xm);gq=zq.*xm+u;dq=(g.*fq-f.*gq)./(f.^2+g.^2);
    Tq=2.*(z+q_.*zq+dq./y)./u;dTq=(3*xm.*Tq+4*q_.^2.*xm./z.*(3-qz.*zq))./u;

    xq=-dTq./d2T;%dx/dq along dT/dx=0
    
    qz2=qz.^2;qz=qz.*qz2;zq=zq+q_.^2.*xm./z.*xq;dqz=3*qz2.*(1-q_./z.*zq)./z;
    d2Tq=(3*Tq+4*(dqz.*s_-2*q_.*qz)+2*d2T.*xm.*xq)./u;%d2T=(3*Tm+4*qz.*s_)./u;

    ii=1:numel(q_)-1;c0=xm(ii);c1=xq(ii);Dx=diff(xm);Ddx=diff(xq);Dq=diff(q_);Dq2=Dq.^2;
    a_=Dx-c1.*Dq;b_=Dq.*Ddx;c2=(3*a_-b_)./Dq2;c3=(-2*a_+b_)./Dq2./Dq;
    %n_=numel(q_)-1;c2(n_)=c2(n_-1)+3*c3(n_-1).*Dq2(n_-1);c3(n_)=(Dx(n_)-Dq(n_)*c1(n_)-Dq2(n_)*c2(n_))/Dq2(n_)/Dq(n_);
    pp.form='pp';pp.breaks=q_;pp.pieces=numel(q_)-1;pp.order=4;pp.dim=1;
    pp.coefs=[c3;c2;c1;c0].';%pp.coefs=[c3 ; c2-c3*3.*q1 ; c1+c3*3.*q12-c2*2.*q1 ; c0-c1.*q1+c2.*q12-c3.*q12.*q1].';%q1=q_(ii);q12=q1.^2;

    q2_=(q_(1:end-1)+q_(2:end))/2;
    [xm2 Tm2 d2T2]=minT(q2_,repmat(N_,size(q2_)),sin(acos(q2_)).^2);
    err=ppval(pp,q2_)-xm2;q2_(abs(err)<3e-7)=[];

    end%mesh_count
%     c0=[xm;Tm;d2T];c1=[xq;Tq;d2Tq];Dx=diff(c0')';Ddx=diff(c1')';ii=1:numel(q_)-1;c0=c0(:,ii);c1=c1(:,ii);Dq=repmat(diff(q_),3,1);Dq2=Dq.^2;
%     a_=Dx-c1.*Dq;b_=Dq.*Ddx;c2=(3*a_-b_)./Dq2;c3=(-2*a_+b_)./Dq2./Dq;
%     pp.form='pp';pp.breaks=q_;pp.pieces=numel(q_)-1;pp.order=4;pp.dim=3;pp.coefs=[c3(:) c2(:) c1(:) c0(:)];

    Tm_data(3,N_)=min(abs([qx q2_]));
    if isempty(Tm_fit);Tm_fit=pp;else;Tm_fit(N_)=pp;end

%     q_=linspace(0,pi,4e4+1);s_=sin(q_).^2;q_=cos(q_);[xm_ Tm_ d2T_]=minT(q_,repmat(N_,size(q_)),s_);
%     figure(1),plot(q_,xm_,q_,ppval(pp,q_))
%     figure(2),plot(q_,ppval(pp,q_)-xm_,pp.breaks,0,'*');set(gca,'xlim',qx*[-1 1])
%     xm=ppval(pp,q_);xm=xm(1,:);
%     u=(1-xm).*(1+xm);y=sqrt(u);z=sqrt(s_+q_.^2.*xm.^2);qz=q_./z;f=y.*(z-q_.*xm);g=z.*xm+q_.*u;
%     Tm=2*((N_*pi+acos(g./sqrt(f.^2+g.^2)))./y+q_.*z-xm)./u;d2T=(3*Tm+4*qz.*qz.*qz.*s_)./u;
%     figure(3),plot(q_,d2T-d2T_);
%     dT=(3*xm.*Tm+4*q_.*q_.^2.*xm./z-4)./u;d2T=d2T+5*xm.*dT./u;

end%for N_
end%if Nx

aq=abs(q);
for N_=min(aN):Nx
ii=aN==N_;if ~any(ii);continue;end
q_=q(ii);s_=cs(ii);xm=ppval(Tm_fit(N_),q_);
u=(1-xm).*(1+xm);y=sqrt(u);z=sqrt(s_+q_.^2.*xm.^2);qz=q_./z;f=y.*(z-q_.*xm);g=z.*xm+q_.*u;
Tm=2*((N_*pi+acos(g./sqrt(f.^2+g.^2)))./y+q_.*z-xm)./u;d2T=(3*Tm+4*qz.*qz.*qz.*s_)./u;
ii=index(ii);xm_(ii)=xm;Tm_(ii)=Tm;d2T_(ii)=d2T;
end
ii=aq>Tm_data(3,aN);if ~any(ii);return;end
q=q(ii);N=N(ii);cs=cs(ii);index=index(ii);xm=xm_(index);Tm_(index)=2*pi*aN(ii);
end%nargin>3

nn=numel(q);
if nargin<4
xm=1./(1.5*(N+.5)*pi);thr2=atan22(cs,2*q)/pi;
aa=find(thr2<.5);if ~isempty(aa);xm(aa)=rt8(2*thr2(aa)).*xm(aa);end
aa=find(thr2>.5);if ~isempty(aa);xm(aa)=(2-rt8(2-2*thr2(aa))).*xm(aa);end
end
Dxm=1;I=0;
if nn>1e6;clear aa thr2;end
while Dxm>3e-7^2&I<12
    I=I+1;
    if nn>1e6;clear Tm dT d2T d3T;end
    [Tm dT d2T d3T]=Tx(xm,q,N,cs);
    dT=dT.*d2T./(d2T.^2-dT.*d3T/2);
    %dT=dT./d2T;
    %dT_=dT;%if abs(dT_)>.5;dT_=sign(dT_)*.5;end
    xm=xm-dT;
    Dxm=max(abs(dT./max(xm,eps)));
end
ii=~(Tm>Tm_(index)-1e-14);%can be NaN
if any(ii);%warning('%d out of %d Tmin cases didn''t converge',sum(ii),numel(ii));
%    xm(ii)=xm_(index(ii));[Tm(ii),~,d2T(ii),~]=Tx(xm(ii),q(ii),N(ii),cs(ii));
    [xm(ii),Tm(ii),d2T(ii)]=minT(q(ii),N(ii),cs(ii));
    ii=~(Tm>Tm_(index)-1e-14);if any(ii);warning('%d out of %d Tmin cases didn''t converge',sum(ii),numel(ii));end
end
if nn<numel(xm_);xm_(index)=xm;Tm_(index)=Tm;d2T_(index)=d2T;
else;xm_=xm;Tm_=Tm;d2T_=d2T;end
return

% function [xm Tm d2T]=minT(q,N,cs)
% thr2=atan2(cs,2*q)/pi;
% xm=1./(1.5*(N+.5)*pi);
% aa=find(thr2<.5);if ~isempty(aa);xm(aa)=rt8(2*thr2(aa)).*xm(aa);end
% aa=find(thr2>.5);if ~isempty(aa);xm(aa)=(2-rt8(2-2*thr2(aa))).*xm(aa);end
% Dxm=1;I=0;
% while Dxm>3e-7&&I<12
%     I=I+1;
%     [Tm dT d2T d3T]=Tx(xm,q,N,cs);
%     dT=dT.*d2T./(d2T.^2-dT.*d3T/2);
%     xm=xm-dT;
%     Dxm=max(abs(dT./max(xm,eps)));
% end
% return

function a=atan22(y,x)
%a=2*atan((sqrt(x.^2+y.^2)-x)./y);
a=2*atan(y./(sqrt(x.^2+y.^2)+x));
ii=real(y)==0;if any(ii);a(ii&x>=0)=0;a(ii&x<0)=pi;end
return

function x=rt8(x)
x=sqrt(sqrt(sqrt(x)));
return
