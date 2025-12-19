function [V v dV]=unit(V)
%[U v dU]=unit(V)
%U  = unit vector of V
%v  = magnitude of V
%dU = derivative of U wrt V, dimension = [3 n 3]
d3=find(size(V)==3)==2;if d3;V=V.';end%work in column vectors
v=mag(V);V=V./repmat(v,size(V,1),1);
ii=v==0;if any(ii);V(1,ii)=1;V(2:3,ii)=0;end%no divide by 0
if nargout>2%derivative
dV=zeros([size(V) 3]);for ii=1:3;dV(ii,:,ii)=1./v;end
dV=dV-bsxfun(@times,V./repmat(v,size(V,1),1),permute(V,[3 2 1]));
end
if d3;V=V.';if nargout>2;dV=permute(dV,[3 2 1]);end;end%flip back to rows
return
