clear all
A=imread('1.jpg');%�Ҷ�ͼƬ
B=imread('2.jpg');%��ɫͼƬ
A=A(:,:,1);
A=histeq(A);
[l1,m1,s1]=match(A,B);
C=lms2rgb(l1,m1,s1);
subplot(1,3,1),imshow(A)
title('�Ҷ�Ŀ��ͼ��')
subplot(1,3,2),imshow(C);
title('Welsh�㷨����')
subplot(1,3,3),imshow(B)
title('��ɫԴͼ��')
%%
function [l,m,s]=rgb2lms(A) %RGBתl����
R=A(:,:,1);
G=A(:,:,2);
B=A(:,:,3);
L=double(0.3811*R+0.5783*G+0.0402*B);
M=double(0.1967*R+0.7244*G+0.0782*B);
S=double(0.0241*R+0.1288*G+0.8444*B);
l=0.5774*L+0.5774*M+0.5774*S;
m=0.4082*L+0.4082*M-0.8165*S;
s=0.7071*L-0.7071*M;
end
%%
function P=LD(l)  %��(��L+��D)/2
[M,N]=size(l);
l=double(l);
L=zeros(M,N);D=L;
for i=3:M-2
    for j=3:N-2
        sq1=l(i-2:i+2,j-2:j+2);
        L(i,j)=mean(sq1(:));
        D(i,j)=std(sq1(:),1);
    end
    sq2=l(i-2:i+2,1:3);
    L(i,1)=mean(sq2(:));
    D(i,1)=std(sq2(:),1);
    sq3=l(i-2:i+2,1:4);
    L(i,2)=mean(sq3(:));
    D(i,2)=std(sq3(:),1);
    sq4=l(i-2:i+2,end-3:end);
    L(i,end-1)=mean(sq4(:));
    D(i,end-1)=std(sq4(:),1);
    sq5=l(i-2:i+2,end-2:end);
    L(i,1)=mean(sq5(:));
    D(i,1)=std(sq5(:),1);
end
for j=3:N-2
    sq6=l(1:3,j-2:j+2);
    L(1,j)=mean(sq6(:));
    D(1,j)=std(sq6(:),1);
    sq7=l(1:4,j-2:j+2);
    L(2,j)=mean(sq7(:));
    D(2,j)=std(sq7(:),1);
    sq8=l(end-3:end,j-2:j+2);
    L(end-1,j)=mean(sq8(:));
    D(end-11,j)=std(sq8(:),1);
    sq9=l(end-2:end,j-2:j+2);
    L(end,j)=mean(sq9(:));
    D(end,j)=std(sq9(:),1);
end
%����
a=l(1:3,1:3);L(1,1)=mean(mean(a));D(1,1)=std(a(:),1);
a=l(1:4,1:3);L(2,1)=mean(mean(a));D(4,1)=std(a(:),1);
a=l(1:3,1:4);L(1,2)=mean(mean(a));D(1,2)=std(a(:),1);
a=l(1:4,1:4);L(2,2)=mean(mean(a));D(2,2)=std(a(:),1);
%����
a=l(M-2:M,1:3);L(M,1)=mean(mean(a));D(M,1)=std(a(:),1);
a=l(M-3:M,1:3);L(M-1,1)=mean(mean(a));D(M-1,1)=std(a(:),1);
a=l(M-2:M,1:4);L(M,2)=mean(mean(a));D(M,1)=std(a(:),1);
a=l(M-3:M,1:4);L(M-1,2)=mean(mean(a));D(M-1,2)=std(a(:),1);
%����
a=l(1:3,N-2:N);L(1,N)=mean(mean(a));D(1,N)=std(a(:),1);
a=l(1:4,N-2:N);L(2,N)=mean(mean(a));D(2,N)=std(a(:),1);
a=l(1:3,N-3:N);L(1,N-1)=mean(mean(a));D(1,N-1)=std(a(:),1);
a=l(1:4,N-3:N);L(2,N-1)=mean(mean(a));D(2,N-1)=std(a(:),1);
%����
a=l(M-2:M,N-2:N);L(M,N)=mean(mean(a));D(M,N)=std(a(:),1);
a=l(M-3:M,N-2:N);L(M-1,N)=mean(mean(a));D(M-1,N)=std(a(:),1);
a=l(M-2:M,N-3:N);L(M,N-1)=mean(mean(a));D(M,N-1)=std(a(:),1);
a=l(M-3:M,N-3:N);L(M-1,N-1)=mean(mean(a));D(M,N-1)=std(a(:),1); 

P=(L+D)/2;
end
%%
function [l1,m1,s1]=match(A,B) %ƥ�����ص㣬AΪ�Ҷ�ͼ��BΪ��ɫͼ
[M,N]=size(A);
[l,m,s]=rgb2lms(B);
A=double(A);
meanA=mean(mean(A));meanB=mean(mean(l));
stdA=std(A(:),1);stdB=std(l(:),1);
A1=stdB/stdA*(A-meanA)+meanB;
l1=A1;m1=zeros(M,N);s1=m1;
PA=LD(A1);PB=LD(l);
for i=1:M
    for j=1:N
        PAB=abs(PB-PA(i,j));
        [x,y]=find(PAB==min(min(PAB)));
        m1(i,j)=m(x(round(1)),y(round(1)));
        s1(i,j)=s(x(round(1)),y(round(1)));
    end
end
end
%%
function C=lms2rgb(l,m,s)  %l����תRGB
L=0.5774*l+0.4082*m+0.7071*s;
M=0.5774*l+0.4082*m-0.7071*s;
S=0.5774*l-0.8165*m;
R=4.4679*L-3.5873*M+0.1193*S;
G=-1.2186*L+2.3809*M-0.1624*S;
B=0.0497*L-0.2439*M+1.2045*S;
[m,n]=size(l);
C=zeros(m,n,3,'uint8');
C(:,:,1)=uint8(R);
C(:,:,2)=uint8(G);
C(:,:,3)=uint8(B);
end