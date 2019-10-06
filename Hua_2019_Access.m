function varargout = Hua_2019_Access( P,para,K )

%% 1. Initialization
% 1.1. Genereate Random Key if K is not given
if ~exist('K','var') && strcmp(para,'encryption') 
    K = round(rand(1,256));
    varOutN = 2;
elseif ~exist('K','var') && strcmp(para,'decryption') 
    error('Cannot Complete Decryption without Encryption Key')
else
    varOutN = 1;
end
KPA=zeros(2,120);
S=sum(K(241:256));
K2=K(121:240);
K2(105:120)=xor(K(225:240),K(241:256));
KPA(1,:)=xor(K(1:120),circshift(K2',S)');
K1=K(1:120);
K1(105:120)=xor(K(105:120),K(241:256));
KPA(2,:)=xor(K(121:240),circshift(K1',S)');

format long eng
C = double(P);
switch para
    case 'encryption'
        for iter = 1:2
             C= permutation(C,'encryption',KPA(iter,:));
            C = Diffusion1(C,'encryption',KPA(iter,:));          
        end
        
      case 'decryption'
        for iter = 2:-1:1
             C = Diffusion1(C,'decryption',KPA(iter,:));
              C = permutation(C,'decryption',KPA(iter,:)); 
        end
end
%% 4. Output

       C = uint8(C);


switch varOutN
    case 1
        varargout{1} = C;
    case 2
        varargout{1} = C;
        varargout{2} = K;
end


function output = permutation( P,para ,K )
C = P;
[M,N]=size(P);
A=1:M;
B=1:N;
RowIndex=1:M;
MP=mod(bi2de(K(1:8),'left-msb'),M)+1;
NP=mod(bi2de(K(9:16),'left-msb'),N)+1;
MStep=bi2de(K(17:20),'left-msb')+1;
NStep=bi2de(K(21:24),'left-msb')+1;

for k=1:M
RowIndex(k)=A(MP);
 A(MP)=[];
 MP=mod(MP-2+MStep,M-k)+1;
 MStep=MStep+1;
end
switch para
    case 'encryption'
        for i=1:M
            B=1:N;
            for j=1:N
                column=B(NP);
                row=mod(RowIndex(mod(column-1,M)+1)+i-2,M)+1;
                B(NP)=[];
                NP=mod(NP-2+NStep,N-j)+1;
                NStep=NStep+1;
                C(row,column)=P(i,j);
            end
                NP=column;
        end      
    case 'decryption'
        for i=1:M
            B=1:N;
            for j=1:N
%                 row=mod(RowIndex(mod(j-1,M)+1)+i-2,M)+1;
                column=B(NP);
                row=mod(RowIndex(mod(column-1,M)+1)+i-2,M)+1;
                B(NP)=[];
                NP=mod(NP-2+NStep,N-j)+1;
                NStep=NStep+1;
                C(i,j)=P(row,column);
            end
            NP=column;
        end
end

        output = C;

        function C = Diffusion1( P, para, K )
 A=zeros(1,8);
 %transFrac = @(K,st,ed) bi2de( xor(K(st:fix((st+ed)/2)),K(fix((st+ed)/2)+1:ed)),'left-msb') ;
 transFrac = @(K,st,ed) bi2de(K(st:ed),'left-msb') ;
        for i=1:8
             A(i) = transFrac(K,(24+12*(i-1)+1),24+12*i);
         end

 [B,IX]=sort(A);
 W1=A(1)+find(IX==1);
 W2=A(2)+find(IX==2);
 W3=A(3)+find(IX==3);
 W4=A(4)+find(IX==4);
 W5=A(5)+find(IX==5);
 W6=A(6)+find(IX==6);
 W7=A(7)+find(IX==7);
 W8=A(8)+find(IX==8);
%     mask =[2^32, 2^32+1,2^32+3;2^32+4,2^32+5,2^32+6;2^32+7,2^32+8,1]
%     invertmask=[2^32,2^32+1,2^32+3;2^32+4,2^32+5,2^32+6;2^32+7,2^32+8,0]
  mask=[W1,W2,W3;W4,W5,W6;W7,W8,1];
  invertmask=[W1,W2,W3;W4,W5,W6;W7,W8,0];
format long eng
switch para
    case 'encryption'
        CI= double(P);
        [M,N]=size(CI);
        for i= 1:M
            for j= 1:N
                  block=[CI(mod(i-3,M)+1,mod(j-3,N)+1),CI(mod(i-3,M)+1,mod(j-2,N)+1),CI(mod(i-3,M)+1,mod(j-1,N)+1);
                      CI(mod(i-2,M)+1,mod(j-3,N)+1),CI(mod(i-2,M)+1,mod(j-2,N)+1),CI(mod(i-2,M)+1,mod(j-1,N)+1);
                      CI(mod(i-1,M)+1,mod(j-3,N)+1),CI(mod(i-1,M)+1,mod(j-2,N)+1),CI(mod(i-1,M)+1,mod(j-1,N)+1)];
                  CI(i,j) = mod(sum(sum(double(block).*double(mask)))+j,256);
            end
        end
        C=uint8(CI);
    case 'decryption'
        D=double(P);
        [M,N]=size(D);
        for i= M:-1:1
            for j= N:-1:1
                  block=[D(mod(i-3,M)+1,mod(j-3,N)+1),D(mod(i-3,M)+1,mod(j-2,N)+1),D(mod(i-3,M)+1,mod(j-1,N)+1);
                     D(mod(i-2,M)+1,mod(j-3,N)+1),D(mod(i-2,M)+1,mod(j-2,N)+1),D(mod(i-2,M)+1,mod(j-1,N)+1);
                     D(mod(i-1,M)+1,mod(j-3,N)+1),D(mod(i-1,M)+1,mod(j-2,N)+1),D(mod(i-1,M)+1,mod(j-1,N)+1)];
                  D1=D(i,j);
                  D(i,j) = mod(D1-sum(sum(double(block).*double(invertmask)))-j,256);
            end
        end
        C=uint8(D);
end



