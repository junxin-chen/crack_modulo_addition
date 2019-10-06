function varargout = Hua_2018SP_LSC_dec( iMat,rKey )
% 解密过程
% iMat      已加密的图像矩阵
% m,n       图像矩阵大小
% rKey      解密密钥

[m,n]=size(iMat);
p = m*n;
% if(max(iMat(:))>1)
%     F = 256;
% else
%     F = 2;
% end
F = 256;


iMat = double(iMat);


[x, y, r1, r2, r3, r4] = GenKey( rKey );
seq1 = SequenceGenerator( p,r1,x,y );
seq2 = SequenceGenerator( p,r2,seq1(1,p),seq1(2,p));
seq3 = SequenceGenerator( p,r3,seq2(1,p),seq1(2,p));
seq4 = SequenceGenerator( p,r4,seq3(1,p),seq1(2,p));

chaoMP1 = reshape(seq1(1,:),m,n);
chaoMD1 = reshape(seq1(2,:),m,n);
chaoMP2 = reshape(seq2(1,:),m,n);
chaoMD2 = reshape(seq2(2,:),m,n);
chaoMP3 = reshape(seq3(1,:),m,n);
chaoMD3 = reshape(seq3(2,:),m,n);
chaoMP4 = reshape(seq4(1,:),m,n);
chaoMD4 = reshape(seq4(2,:),m,n);

u1 = Diffusion( chaoMD4, iMat, 'decryption',F );
u2 = Permutation( chaoMP4, u1, 'decryption' );
u3 = Diffusion( chaoMD3, u2, 'decryption',F );
u4 = Permutation( chaoMP3, u3, 'decryption' );
u5 = Diffusion( chaoMD2, u4, 'decryption',F );
u6 = Permutation( chaoMP2, u5, 'decryption' );
u7 = Diffusion( chaoMD1, u6, 'decryption',F );
plain = Permutation( chaoMP1, u7, 'decryption' );





switch F
    case 2
        P = logical(plain);
    case 256
        P = uint8(plain);
end

varargout{1} = P;



end




function [x, y, r1, r2, r3, r4] = GenKey( key )
% Generate the parameters of the chaotic system


calculate = @(key,st,ed) sum(key(st:ed).*2.^(-(1:(ed-st+1))));
integer = @(key,st,ed) sum(key(ed:-1:st).*2.^(0:(ed-st)));

if length(key)==256
    x = calculate(key,1,52);
    y = calculate(key,53,104);
    r = calculate(key,105,156);
    A1 = integer(key,157,181);
    A2 = integer(key,182,206);
    A3 = integer(key,207,231);
    A4 = integer(key,232,256);
    
    r1 = mod(r*A1,1);
    r2 = mod(r*A2,1);
    r3 = mod(r*A3,1);
    r4 = mod(r*A4,1);
    

end


end

function chaoSequence = SequenceGenerator( n,r,x,y )
% 混沌序列生成
% n       序列总长度
% r       混沌系统参数
% x,y     混沌系统初值
format long eng

chaoSequence = zeros(2,n);

chaoSequence(1,1) = x;
chaoSequence(2,1) = y;
for i=2:n
    chaoSequence(1,i) = sin(pi*(4*r*x*(1-x)+(1-r)*sin(pi*y)));
    chaoSequence(2,i) = sin(pi*(4*r*y*(1-y)+(1-r)*sin(pi*x)));
    x = chaoSequence(1,i);
    y = chaoSequence(2,i);
end
end



function t = Permutation( chaoM, iMat, para )
% 置乱操作
% chaoM      混沌辅助矩阵
% iMat       图像像素矩阵
% para       参数


[m,n] = size(iMat);
t = zeros(m,n);
ctemp = zeros(1,n);
[~,col] = sort(chaoM,1);
switch para
    case 'encryption'

        for i=1:m
            for j=1:n
                ctemp(1,j) = chaoM(col(i,j),j);
            end
            [~,indc] = sort(ctemp,2);
            for k=1:n
                %t(col(i,indc(1,k)),indc(1,k)) = iMat(col(i,k),k);
                t(col(i,k),k) = iMat(col(i,indc(1,k)),indc(1,k));
            end

        end
    case 'decryption'

        for i=1:m
            for j=1:n
                ctemp(1,j) = chaoM(col(i,j),j);
            end
            [~,indc] = sort(ctemp,2);
            for k=1:n
                %t(col(i,k),k) = iMat(col(i,indc(1,k)),indc(1,k));
                t(col(i,indc(1,k)),indc(1,k)) = iMat(col(i,k),k);
            end

        end
end

% rtemp = zeros(m,1);
% 
% for i=1:n
%     for j=1:m
%         rtemp(j,1) = chaoM(j,row(j,i));
%     end
%     [~,indr] = sort(rtemp,1);
%     for k=1:m
%         cipher(indr(k,1),row(indr(k,1),i)) = t(k,row(k,i));
%     end
%     
% end

        
end

function cipher = Diffusion( chaoM, iMat, para, F )
%  扩散操作  
%  chaoM         混沌辅助矩阵
%  iMat          图像矩阵
%  para          参数
%  F             图像类别

[m,n] = size(iMat);
cipher = zeros(m,n);
cSupMat = mod(floor(chaoM.*2^32),F);
switch para
    case 'encryption'
        t(1,:) = mod(iMat(1,:)+iMat(m,:)+iMat(m-1,:)+cSupMat(1,:),F);
        t(2,:) = mod(iMat(2,:)+t(1,:)+iMat(m,:)+cSupMat(2,:),F);
        for i=3:m
            t(i,:) = mod(iMat(i,:)+t(i-1,:)+t(i-2,:)+cSupMat(i,:),F);
        end
        cipher(:,1) = mod(t(:,1)+t(:,n)+t(:,n-1)+cSupMat(:,1),F);
        cipher(:,2) = mod(t(:,2)+cipher(:,1)+t(:,n)+cSupMat(:,2),F);
        for i=3:n
            cipher(:,i) = mod(t(:,i)+cipher(:,i-1)+cipher(:,i-2)+cSupMat(:,i),F);
        end

    case 'decryption'
        for i=n:-1:3
            t(:,i) = mod(iMat(:,i)-iMat(:,i-1)-iMat(:,i-2)-cSupMat(:,i),F);
        end
        t(:,2) = mod(iMat(:,2)-t(:,n)-iMat(:,1)-cSupMat(:,2),F);
        t(:,1) = mod(iMat(:,1)-t(:,n)-t(:,n-1)-cSupMat(:,1),F);
        
        for i=m:-1:3
            cipher(i,:) = mod(t(i,:)-t(i-1,:)-t(i-2,:)-cSupMat(i,:),F);
        end
        cipher(2,:) = mod(t(2,:)-cipher(m,:)-t(1,:)-cSupMat(2,:),F);
        cipher(1,:) = mod(t(1,:)-cipher(m,:)-cipher(m-1,:)-cSupMat(1,:),F);
        
     
        
end


end


