function [info_tip,seg_info]=change_order(info_tip,seg_info,type)

if nargin <3
   type='def'; 
end

type_max=max(seg_info(:,6));
if strcmp(type,'def')
    switch type_max
        case 20
            seg_info(seg_info(:,6)<19,6)=1;
            seg_info(seg_info(:,6)>=19,6)=2;
            info_tip(info_tip(:,6)<19,6)=1;
            info_tip(info_tip(:,6)>=19,6)=2;
        case 13
            seg_info(seg_info(:,6)<12,6)=1;
            seg_info(seg_info(:,6)>=12,6)=2;
            info_tip(info_tip(:,6)<12,6)=1;
            info_tip(info_tip(:,6)>=12,6)=2;
        case 12
            seg_info(seg_info(:,6)<12,6)=1;
            seg_info(seg_info(:,6)>=12,6)=2;
            info_tip(info_tip(:,6)<12,6)=1;
            info_tip(info_tip(:,6)>=12,6)=2;
        case 10
            seg_info(seg_info(:,6)<9,6)=1;
            seg_info(seg_info(:,6)>=9,6)=2; 
            info_tip(info_tip(:,6)<9,6)=1;
            info_tip(info_tip(:,6)>=9,6)=2; 
        case 5
            seg_info(seg_info(:,6)<=2,6)=1;
            seg_info(seg_info(:,6)>=3 & seg_info(:,6)<=5,6)=2 ; 
            info_tip(info_tip(:,6)<=2,6)=1;
            info_tip(info_tip(:,6)>=3 & info_tip(:,6)<=5,6)=2 ; 
        case 4
            seg_info(seg_info(:,6)<2,6)=1;
            seg_info(seg_info(:,6)>=2 & seg_info(:,6)<3,6)=2 ; 
            seg_info(seg_info(:,6)>=3,6)=3 ;  
            
            info_tip(info_tip(:,6)<2,6)=1;
            info_tip(info_tip(:,6)>=2 & info_tip(:,6)<3,6)=2 ; 
            info_tip(info_tip(:,6)>=3,6)=3 ;  
        case 7
            seg_info(seg_info(:,6)<2,6)=1;
            seg_info(seg_info(:,6)>=2 & seg_info(:,6)<4,6)=2 ;
            seg_info(seg_info(:,6)>=4 & seg_info(:,6)<8,6)=3 ;  
            info_tip(info_tip(:,6)<2,6)=1;
            info_tip(info_tip(:,6)>=2 & info_tip(:,6)<4,6)=2 ;
            info_tip(info_tip(:,6)>=4 & info_tip(:,6)<8,6)=3 ;  
        case 6
            seg_info(seg_info(:,6)<2,6)=1;
            seg_info(seg_info(:,6)>=2 & seg_info(:,6)<4,6)=2 ;
            seg_info(seg_info(:,6)>=4 & seg_info(:,6)<8,6)=3 ; 
            info_tip(info_tip(:,6)<2,6)=1;
            info_tip(info_tip(:,6)>=2 & info_tip(:,6)<4,6)=2 ;
            info_tip(info_tip(:,6)>=4 & info_tip(:,6)<8,6)=3 ; 
        case 1
        case 2
        case 3
        otherwise
            error('Problem with number type max')
    end
elseif strcmp(type,'rootswms')==1
    global lstem
    seg_info(:,6)=seg_info(:,6)-lstem;
    info_tip(:,6)=info_tip(:,6)-lstem;
    seg_info(seg_info(:,6)>=3,6)=3;
    info_tip(info_tip(:,6)>=312,6)=3;
    seg_info(seg_info(:,6)==0,6)=1;
    info_tip(info_tip(:,6)==0,6)=1;
elseif strcmp(type,'Mil_sixtine')==1
    seg_info(:,6)=seg_info(:,6);
    info_tip(:,6)=info_tip(:,6);
else
    switch type
        case 'sathyan'
            seg_info(seg_info(:,6)==4,6)=2;
            seg_info(seg_info(:,6)==5,6)=3;
            info_tip(info_tip(:,6)==4,6)=2;
            info_tip(info_tip(:,6)==5,6)=3;
            
        case 'MRT'
            seg_info(seg_info(:,6)<12,6)=1;
            seg_info(seg_info(:,6)>=12,6)=2;
            info_tip(info_tip(:,6)<12,6)=1;
            info_tip(info_tip(:,6)>=12,6)=2;
        case 'maize_RB'
            seg_info(seg_info(:,6)>3,6)=1;
            info_tip(info_tip(:,6)>3,6)=1;
        case 'wheat'
            seg_info(seg_info(:,6)<19,6)=1;
            seg_info(seg_info(:,6)>=19,6)=2;
            info_tip(info_tip(:,6)<19,6)=1;
            info_tip(info_tip(:,6)>=19,6)=2;
        case 'maize'
            seg_info(seg_info(:,6)<3,6)=1;
            seg_info(seg_info(:,6)>=3,6)=2;
            info_tip(info_tip(:,6)<3,6)=1;
            info_tip(info_tip(:,6)>=3,6)=2;
            
        case 'fibrous'
            seg_info(:,end)=distance2tip(seg_info,info_tip); 
        case 'taproot'
            seg_info(:,end)=distance2tip(seg_info,info_tip); 
        case 'lolium'
            seg_info(seg_info(:,6)<=2,6)=1;
            seg_info(seg_info(:,6)>=3 & seg_info(:,6)<=5,6)=2 ; 
            info_tip(info_tip(:,6)<=2,6)=1;
            info_tip(info_tip(:,6)>=3 & info_tip(:,6)<=5,6)=2 ; 
        case 'simple_lat' 
            seg_info(seg_info(:,6)>=3,6)=2 ; 
            info_tip(info_tip(:,6)>=3,6)=2 ;    
        case 'arab' 
            
            seg_info(seg_info(:,6)<3,6)=1 ; 
            seg_info(seg_info(:,6)>=3 & seg_info(:,6)<=4,6)=2 ;
            seg_info(seg_info(:,6)>4,6)=3 ; 

            info_tip(info_tip(:,6)<3,6)=1 ; 
            info_tip(info_tip(:,6)>=3 & info_tip(:,6)<=4,6)=2 ; 
            info_tip(info_tip(:,6)>4,6)=3 ;    
            
        case 'any'
            seg_info(seg_info(:,6)>2 & seg_info(:,6)<=type_max,6)=3 ; 
            info_tip(info_tip(:,6)>2 & info_tip(:,6)<=type_max,6)=3 ; 
            
        case 'xray'
            seg_info(seg_info(:,6)<3,6)=1 ; 
            seg_info(seg_info(:,6)>=3 & seg_info(:,6)<=3,6)=2 ;
            seg_info(seg_info(:,6)>3,6)=3 ; 

            info_tip(info_tip(:,6)<3,6)=1 ; 
            info_tip(info_tip(:,6)>=3 & info_tip(:,6)<=3,6)=2 ; 
            info_tip(info_tip(:,6)>3,6)=3 ;    
                       
        case 'xray2'
            seg_info(seg_info(:,6)<3,6)=1 ; 
            seg_info(seg_info(:,6)>=3 & seg_info(:,6)<=3,6)=2 ;
            seg_info(seg_info(:,6)>3,6)=3 ; 

            info_tip(info_tip(:,6)<3,6)=1 ; 
            info_tip(info_tip(:,6)>=3 & info_tip(:,6)<=3,6)=2 ; 
            info_tip(info_tip(:,6)>3,6)=3 ;                     
        case 'maize_rootbox'
            seg_info(seg_info(:,6)>3,6)=1 ; 
            info_tip(info_tip(:,6)>3,6)=1 ; 
            
        case 'maize_mutez'
            seg_info(seg_info(:,6)==4,6)=1 ; 
            info_tip(info_tip(:,6)==4,6)=1 ; 
            
            seg_info(seg_info(:,6)==6,6)=2 ; 
            info_tip(info_tip(:,6)==6,6)=2 ; 
            
            seg_info(seg_info(:,6)==5,6)=3 ; 
            info_tip(info_tip(:,6)==5,6)=3 ; 
            
        case 'maize_mutez_mohsen'
            
            seg_info_mod=seg_info;
            tip_info_mod=info_tip;
            
            seg_info_mod(seg_info(:,6)==1,6)=1 ; % primary
            tip_info_mod(info_tip(:,6)==1,6)=1 ; 
            
            seg_info_mod(seg_info(:,6)==4,6)=1 ; % Seminal
            tip_info_mod(info_tip(:,6)==4,6)=1 ; 
            
            seg_info_mod(seg_info(:,6)==2,6)=4 ; % Short laterals
            tip_info_mod(info_tip(:,6)==2,6)=4 ; 
            
            seg_info_mod(seg_info(:,6)==6,6)=4 ; % Long laterals
            tip_info_mod(info_tip(:,6)==6,6)=4 ; 
            
            for i=1:size(info_tip,1)
               if info_tip(i,6)==5;
                   pos=find(seg_info(:,7)==i);
                   if (seg_info(pos(1),4))<0
                       tip_info_mod(i,6)=2;
                       seg_info_mod(pos,6)=2;
                   else
                       tip_info_mod(i,6)=3;
                       seg_info_mod(pos,6)=3;
                   end
               end
            end
                
            seg_info=seg_info_mod;
            info_tip=tip_info_mod;
            
        case 'no_change'    
            
    end
    
end

end