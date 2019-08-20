clear all
A = imread('plate_number.png');
P = im2bw(A,0.7);  % ��ֵ��
figure(),imshow(P)
title('��ֵ�����','FontSize',12)
P = bwareaopen(P,2000);  % ��С���ȥ��
I1 = bwareaopen(P,5000); 
P = P-I1;  % �Ѵ����ȥ��
figure(),imshow(P)
title('ȥ����','FontSize',12)
[x y] = find(P ==1);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
rects = [ymin-1,xmin-1,ymax-ymin+2,xmax-xmin+2];  % ����������ı߽��������浽rects
figure, imshow(A)
rectangle('Position',rects,'EdgeColor','r','LineWidth',1);  % ��λ�������򣬲��ú�ɫ�Ŀ���
title('���ƶ�λ','FontSize',12)
P = imcrop(A,rects);  % ���պ��߿��и������
figure()
subplot(1,2,1),imshow(P)
title('У��ǰ')
% ˮƽ�������
T=affine2d([0 1 0;1 0 0;0 0 1]);
I2=imwarp(P,T);  % ͼ��ת�ã�˳ʱ����ת90�����ˮƽ����
theta = -20 : 20;  % ������б�Ƕȵķ�Χ
r1 = radon(I2, theta);  % radon�任ȷ����б��
result1 = sum(abs(diff(r1)), 1);  % ����е�������ֵ���ۼӺͣ����Ķ�Ӧ��б��
rot1 = find(result1==max(result1))-21;
P = imrotate(P, rot1);
% ��ֱ�������
r2 = radon(P, theta);
result2 = sum(abs(diff(r2)), 1);
rot2 = (find(result2==max(result2))-21)/57.3;  % ����ֵתΪ�Ƕ�
if rot2>0
    T1 = affine2d([1 0 0 ; -tan(rot2) 1 0 ; size(P, 1) * tan(rot2) 0 1]);
else
    T1 = affine2d([1 0 0 ; tan(-rot2) 1 0 ; size(P, 1) * tan(-rot2) 0 1]);
end
P = imwarp(P, T1);
subplot(1,2,2), imshow(P)
title('У����')
suptitle('��бУ��')
P = im2bw(P,0.7);  % ��ֵ��
se = [1;1;1];
P = imerode(P,se);
[x y] = find(P ==1);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
rects = [ymin,xmin,ymax-ymin,xmax-xmin];
P = imcrop(P,rects);  % ��ֵ��
P = P(3:end-2,6:end-5);  % �ñ�
figure, imshow(P)
title('��ȷ��λ','FontSize',12)
[m n] = size(P);
m = ceil(m/2); n = ceil(n/4);
B1 = P(m+1:end,1:n);
B2 = P(m+1:end,n+1:2*n);
B3 = P(m+1:end,2*n+1:3*n);
B4 = P(m+1:end,3*n+1:end);
figure()
subplot(1,4,1),imshow(B1)
subplot(1,4,2),imshow(B2)
subplot(1,4,3),imshow(B3)
subplot(1,4,4),imshow(B4)
suptitle('���ַָ�')
B = [B1(:) B2(:) B3(:) B4(:)];
temp = zeros(28,28);
C = zeros(784,4);
figure()
for i = 1:4
    img = B(:,i);
    img = reshape(img,size(B1));
    [x y] = find(img == 0);
    img = img(min(x):max(x),min(y):max(y));
    img = imresize(img,[20,20]);
    temp(5:24,5:24) = 1 - img;
    subplot(1,4,i),imshow(temp)
    C(:,i) = temp(:);
end
suptitle('��һ��')
[train,trainlabels,test,testlabels] = read_minst();
train = imbinarize(train,0.3); 
[~,I] = pdist2(double(train'),C','euclidean','Smallest',25);
knn_label = trainlabels(I);
predict = mode(knn_label);
figure(),imshow(A)
title(['ʶ������',int2str(predict)],'FontSize',12)