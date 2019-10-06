% This is a function used in Mathematical Problems in Engineering
% Volume 2009, Article ID 762652, 22 pages
% This function is to produce random bit sequnce using logistic map


function image = Boru2009MPE(p,para)
% pixel_count=6666;

% the control parameters and initial values of the logistic map and tent
% map serve as the key
% it is proposed to be produced by secret bits, but the production process is note presented.
% However, it is not important,thus we direct use their valid values

miu=3.9993;
x0=0.66;
tent_x0=0.77;
tent_p=0.44;

[row_count,colum_count]=size(p);

add_mod=@(a,b)uint8(mod(double(a)+double(b),256));% modulo addition
inv_add_mod=@(c_sum,b)uint8(mod(double(c_sum)-double(b),256)); % inverse of modulo addition

%% This the permutation part 
% preparation for the permutation
order=row_count:-1:2;
gi=zeros(1,row_count-1);
ji=zeros(1,row_count-1);
ki=zeros(1,row_count-1);

for i=1:row_count-1
    ji(1,i)=floor(log2(i))+1;
end
for i=3:row_count-1
    for s=2:i-1
        ki(1,i)=ki(1,i)+floor(log2(s))+1;
    end
end
required_inter=ki(1,row_count-1)+ji(1,row_count-1);

% this is the iteration of the chaotic map
rad_bits_logi=zeros(1,required_inter);

for i=1:300
    x0=miu*x0*(1-x0);
end
for i=1:required_inter   
    x0=miu*x0*(1-x0);
    rad_bits_logi(1,i)=x0>0.6;
end
rad_bits_logi=double(rad_bits_logi);

% this is the production of gi
gi(1,1)=1;
gi(1,2)=floor((2*rad_bits_logi(1)+rad_bits_logi(2))/3)+1;
for i=3:row_count-1
    tj=ji(i);
    tk=ki(i);
    tmp_fenzi=0;
    for kk=1:tj
        tmp_fenzi=tmp_fenzi+2^(tj-kk)*rad_bits_logi(tk+kk-1);
    end
    tmp_fenzi=tmp_fenzi*(i-1);
    tmp_fenzi=floor(tmp_fenzi/(2^tj-1));
    gi(1,i)=tmp_fenzi+1;
end

% this is the production of the permutation matrix
per_vec=1:1:row_count;
%Tompkins-Paige Algorithm
for i=1:row_count-1
    per_vec_front=per_vec(1:i-1);
    per_vec_rear=per_vec(i:row_count);
    per_vec_rear=circshift(per_vec_rear,[0,gi(row_count-i)]);
    per_vec=[per_vec_front per_vec_rear];
end

%% This is the substitution part
for i=1:300
    if tent_x0<tent_p
        tent_x0=tent_x0/tent_p;
    else
        tent_x0=(1-tent_x0)/(1-tent_p);
    end
end

mask=zeros(row_count,colum_count);
for i=1:row_count
    for j=1:colum_count 
        if tent_x0<tent_p
        tent_x0=tent_x0/tent_p;
        else
        tent_x0=(1-tent_x0)/(1-tent_p);
        end
        mask(i,j)=floor(tent_x0*256);
    end
end

%% start encryption or decryption
image=zeros(row_count,colum_count);
img_tmp=zeros(row_count,colum_count);
switch para
    case 'en'  % this is encryption
        % rwo shuffling
        for i=1:row_count
            img_tmp(i,:)=p(per_vec(i),:);
        end
        % column shuffling
        for j=1:colum_count
            image(:,j)=img_tmp(:,per_vec(j));
        end
        
        % substitution
        for i=1:row_count
            for j=1:colum_count
                image(i,j)=add_mod(image(i,j),mask(i,j));
            end
        end        
        
    case 'de'
        % inverse of substitution
        for i=1:row_count
            for j=1:colum_count
                p(i,j)=inv_add_mod(p(i,j),mask(i,j));
            end
        end   
        
        % inverse of column shuffling
        for j=1:colum_count
            img_tmp(:,per_vec(j))=p(:,j);
        end
        
        for i=1:row_count
            image(per_vec(i),:)=img_tmp(i,:);
        end
end

image=uint8(image);
