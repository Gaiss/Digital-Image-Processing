%% ��һ warp����
clear all
A = imread('Colorful Rose.jpg');
xx = 0:0.01:10;
yy = -10:0.01:10;
[x1 y1] = meshgrid(xx,yy);
z1 = x1.^2 + y1.^2;
z1(z1<10) = NaN;
z1(z1>100) = NaN;
% 
% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2),warp(x,y,z,A);
% view(3)
% axis square
% grid on
% title('warp����')
% 
% figure()
% subplot(2,2,1),warp(x,y,z,A)
% title('����ͼ')
% view(0,0)
% axis square
% subplot(2,2,2),warp(x,y,z,A)
% title('����ͼ')
% view(90,0)
% axis square
% subplot(2,2,3),warp(x,y,z,A)
% title('����ͼ')
% view(2)
% axis square
%% ���� CData���� ǰ��ӳ��
% clear all
% A = imread('Colorful Rose.jpg');
% [m n ~] = size(A);
% xx = linspace(0,10,n);
% yy = linspace(10,-10,m);
% [x2 y2] = meshgrid(xx,yy);
% z2 = x2.^2 + y2.^2;
% B = z2<=100 & z2>=10;
% z2(B == 0) = NaN;
% 
% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A); % ������ͼ
% view(3)
% axis square
% title('CData���� ǰ��ӳ��')
% 
% figure()
% subplot(2,2,1)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(0,0)
% axis square
% subplot(2,2,2)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(90,0)
% axis square
% subplot(2,2,3)
% h=surf(x1,y1,z1);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(2)
% axis square
%% ���� CData���� ����ӳ��
% clear all
A = imread('Colorful Rose.jpg');
[m n ~] = size(A);
theta = linspace(-pi/2,pi/2,110) ;
r = linspace(sqrt(10),10,73);
[theta r] = meshgrid(theta,r);
x = r.*cos(theta) ;
y = r.*sin(theta) ;
z = x.^2 + y.^2;

% figure()
% subplot(1,2,1),imshow(A)
% title('Original')
% subplot(1,2,2)
% mesh(x,y,z,'EdgeColor','k')
% axis off
% shading interp
% set(h,'CData',A); % ������ͼ
% view(3)
% axis square
% title('CData���� ����ӳ��')
% 
% figure()
% subplot(2,2,1)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(0,0)
% axis square
% subplot(2,2,2)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(90,0)
% axis square
% subplot(2,2,3)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A);
% title('����ͼ')
% view(2)
% axis square
%%
% figure()
% subplot(1,3,1),warp(x1,y1,z1,A);
% view(2)
% axis square
% grid on
% title('warp����')
% subplot(1,3,2)
% h=surf(x2,y2,z2);
% shading interp
% set(h,'CData',A); % ������ͼ
% view(2)
% axis square
% title('CData���� ǰ��ӳ��')
% subplot(1,3,3)
% h=surf(x,y,z);
% shading interp
% set(h,'CData',A); % ������ͼ
% view(2)
% axis square
% title('CData���� ����ӳ��')

% w=5; %���Ʒֿ鴰�ڵĴ�С��16��16
% M=m/w;N=n/w;  %���м�������
% % �����ָ�
% t1=(0:M-1)*w+1;t2=(1:M)*w;
% t3=(0:N-1)*w+1;t4=(1:N)*w;
% figure() 
% imshow(A); 
% hold on
% for i=1:M
%     for j=1:N
%         x=t1(i):t2(i);
%         y=t3(j):t4(j);
%         rectangle ('Position', [t3(j) t1(i) length(x) length(y)], ...,
%             'EdgeColor','r','LineWidth',0.05); %���ƾ��ηָ���ӣ���ɫ
%     end
% end