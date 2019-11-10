%draw_root
% plot the root system in white (usually on a prevsiouly dran FEM 3D plot)
% Javaux, M., 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_root(rootsys,varargin)
% input arg (rootsys,color = yes(1) or no(0) periodic domain = yes (dimensions), no (0), position of the seed)

[info_tip,seg_info,age,info_grw,n_br,n_tip,n_root,xyroot]=read_rootsys(rootsys);
Xr=seg_info(:,2);Yr=seg_info(:,3);Zr=seg_info(:,4);
Xg=info_tip(:,2);Yg=info_tip(:,3);Zg=info_tip(:,4);
diam=seg_info(:,9)./seg_info(:,8);

diam_m=mean(diam);

switch nargin
    case 1
        varargin={0,0,[0,0],'def'} ; %color continu position type
    case 2
        varargin={varargin{1},0,[0,0],'def'} ; 
    case 3
        varargin={varargin{1},varargin{2},[0,0],'def'} ; 
    case 4
        varargin={varargin{1},varargin{2},varargin{3},'def'};
    case 5
    otherwise
        error('wrong number of input parameters (1-5)')     
end

lcontinu=varargin{2};
if lcontinu(1)~=0
   dim_x=lcontinu(1);
   dim_y=lcontinu(2);
   lcontinu=1;
end

pos_x=varargin{3}(1);
pos_y=varargin{3}(2);

                     
hold on;

type=varargin{4};


color=varargin{1};
if color==0
    map=zeros(length(unique(info_tip(:,6))),3);   
else
    [info_tip,seg_info]=change_order(info_tip,seg_info,type);
%     info_tip(isnan(info_tip(:,6)),:)=[];
%     map=flipud(jet(length(unique(info_tip(:,6))))); 
%     N=rand(1,3);
%     map=repmat(N./norm(N),3,1);
    map=hsv(3);
end

% map=repmat([0.4 0.6 1],size(map,1),1);
% map=[0.15 0.15 0.7;0.4 0.6 1];
orders=(unique(info_tip(:,6)));
%relates tip and brenches
% n_br=size(info_tip,1);
switch lcontinu
    case 0
for in=1:n_br
    %disp(in)
    clear brench;clear num;
    OK=find(seg_info(:,7)==in);
    brench(1,:)=[Xr(OK)',Xg(in)]+pos_x;
    brench(2,:)=[Yr(OK)',Yg(in)]+pos_y;
    brench(3,:)=[Zr(OK)',Zg(in)];
    num=[seg_info(OK,1)',info_tip(in,1)];
    matseg=seg_info(OK,:);
    %find embranchment
    if min(seg_info(OK,5))~=0 %not connected to the seed
         OKbrench=find(seg_info(matseg(:,5),7)~=in);
         ibrench=matseg(OKbrench,5);
         brench(1:3,2:end+1)=brench;
         brench(1:3,1)=[Xr(ibrench)+pos_x;Yr(ibrench)+pos_y;Zr(ibrench)];
         num(2:end+1)=num;
         num(1)=ibrench;
    elseif in==1 & sum(brench(:,1)-brench(:,2))==0 %Seed
        %plot3(brench(1,1),brench(2,1),brench(3,1),'og');hold on;
        brench(:,1)=[];
        num(1)=[];
    end
    for j=1:n_root
        xi=brench(1,:)+xyroot(j,2);
        yi=brench(2,:)+xyroot(j,3);
        zi=brench(3,:);
            l_line(in)=line(xi,yi,zi,...
                'color',map(find(orders==info_tip(in,6)),:),'linewidth',1.2);
            hold on;%
    end
%     
    clear brench;clear num;
end
   
    case 1
        
for in=1:n_br
    OK=find(seg_info(:,7)==in);
    brench(1,:)=[Xr(OK)',Xg(in)];
    brench(2,:)=[Yr(OK)',Yg(in)];
    brench(3,:)=[Zr(OK)',Zg(in)];
    num=[seg_info(OK,1)',info_tip(in,1)];
    matseg=seg_info(OK,:);
    %find embranchment
    if min(seg_info(OK,5))~=0 %not connected to the seed
         OKbrench=find(seg_info(matseg(:,5),7)~=in);
         ibrench=matseg(OKbrench,5);
         brench(1:3,2:end+1)=brench;
         brench(1:3,1)=[Xr(ibrench);Yr(ibrench);Zr(ibrench)];
         num(2:end+1)=num;
         num(1)=ibrench;
    elseif in==1 & sum(brench(:,1)-brench(:,2))==0 %Seed
        plot3(brench(1,1),brench(2,1),brench(3,1),'og');hold on;
        brench(:,1)=[];
        num(1)=[];
    end
    for j=1:n_root
        xi=brench(1,:)+xyroot(j,2);
        tr_x=ceil((xi-dim_x./2)./(dim_x));
        xi=-tr_x*dim_x/2+(xi-(tr_x*dim_x/2))+pos_x;
        yi=brench(2,:)+xyroot(j,3);
        tr_y=ceil((yi-dim_y./2)./(dim_y));
        yi=-tr_y*dim_y/2+(yi-(tr_y*dim_y/2))+pos_y;
        zi=brench(3,:);
        diffe=abs(diff(tr_x))+abs(diff(tr_y));
        pos=[0,find(diffe~=0),length(tr_x)];
        for iii=1:length(pos)-1
            l_line(in)=line(xi(pos(iii)+1:pos(iii+1)),yi(pos(iii)+1:pos(iii+1)),zi(pos(iii)+1:pos(iii+1)),...
                   'color',map(info_tip(in,6),:),'linewidth',1.2);
 %              'color',map(find(orders==info_tip(in,6)),:),'linewidth',1.2); 
                 hold on;%
        end
    end
%     
    clear brench;clear num;
end        
        
end
hold off;grid on;
axis equal;%tight
xlabel('X');ylabel('Y');zlabel('Z');
view([45 25])

if color ~= 0
for i=1:length(orders)

    temp=find(info_tip(:,6)==orders(i));
    position(i)=temp(1);
    legende{i}=['Order #',num2str(orders(i))];
end

legend(l_line(position),legende)

end

end
