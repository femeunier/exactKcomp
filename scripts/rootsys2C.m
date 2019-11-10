function [S,prev]=rootsys2C(file,varargin)

if nargin <2
    [~,seg_info]=read_rootsys(file);
else
    seg_info=varargin{1};
end
Nr=size(seg_info,1);
seg_info(:,end+1)=1:Nr;
[~,prev]=ismember(seg_info(:,5),seg_info(:,1));
Im1 = spalloc(Nr+1,Nr,2*Nr);
Im1=[[-1,sparse(zeros(1,Nr-1))];speye(Nr)];
Im1=spsubsasgn (Im1,[prev(2:end,1)+1],[2:Nr]',-1);
Im2=speye(Nr);
S=[Im1,[sparse(zeros(1,Nr));Im2]];


end