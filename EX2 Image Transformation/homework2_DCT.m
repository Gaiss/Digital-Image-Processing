clear all
A=imread('lena.png');
[C,P,Q] = DCT(A);
C_A1 = DCTshift(C);
C_A2 = fftshift(dct2(A));
CP_A = atan2(imag(C_A1),real(C_A1)); % 或FP_A = angle(F_A1)
figure(1)
subplot(2,3,1),imshow(A),title('原图')
subplot(2,3,2),imshow(log(1+abs(C_A1)),[]),title('DCT自编函数')
subplot(2,3,3),imshow(log(1+abs(C_A2)),[]),title('DCT')
subplot(2,3,4),imshow(log(1+abs(real(C_A1))),[]),title('实部')
subplot(2,3,5),imshow(log(1+abs(imag(C_A1))),[]),title('虚部')
subplot(2,3,6),imshow(CP_A),title('相位图')
%% 低频重构
w1 = [8 16 32 64];
[M N] = size(C);
T1 = zeros(size(C));
figure(2)
for i=1:length(w1)
    t = M/2-w1(i)/2+1:M/2+w1(i)/2;
    T1(t,t) = C_A1(t,t);
    T_A = inv(P)*T1*inv(Q);
    subplot(2,2,i),imshow(abs(T_A),[]),title(['窗口大小',num2str(w1(i)),'×',num2str(w1(i))])
end
suptitle('低频重构')
%% 高频重构
w2 = [256 128 64 32];
[M N] = size(C);
T2 = C_A1;
figure(3)
for i=1:length(w2)
    t = M/2-w2(i)/2+1:M/2+w2(i)/2;
    T2(t,t) = 0;
    T_A = inv(P)*T2*inv(Q);
    subplot(2,2,i),imshow(abs(T_A),[]),title(['挖去窗口大小',num2str(w2(i)),'×',num2str(w2(i))])
end
suptitle('高频重构')
%%
function [C,P,Q] = DCT(f)
f=double(f);
[M N] = size(f);
if M~=N
    error('图像必须为方阵！');
end
x = 0:M-1; v = x;
y = (0:N-1)'; u = y;
C = zeros(size(f));
a = sqrt(2/N)*eye(M); a(1,1)=sqrt(1/N);
P = a*cos(u*(2*x+1)*pi/2/M);
Q = a*cos((2*y+1)*v*pi/2/N);
C = P*f*Q;
end
%%
function A_shift = DCTshift(A)
A=double(A);
M = size(A,1);
N = size(A,2);
if log2(M)~=fix(log2(M)) || log2(N)~=fix(log2(N))
    error('图像尺寸必须为2的幂！')
end
A_shift = zeros(size(A));
A_shift(1:M/2,1:N/2) = A(M/2+1:M,N/2+1:N);
A_shift(1:M/2,N/2+1:N) = A(M/2+1:M,1:N/2);
A_shift(M/2+1:M,1:N/2) = A(1:M/2,N/2+1:N);
A_shift(M/2+1:M,N/2+1:N) = A(1:M/2,1:N/2);
end