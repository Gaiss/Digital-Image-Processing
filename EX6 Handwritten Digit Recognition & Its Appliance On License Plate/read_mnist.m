function [train_img,train_label,test_img,test_label] = read_mnist()
train_img = []; test_img = [];
train_label = zeros(60000,1); test_label = zeros(10000,1);
fid = fopen('train-images.idx3-ubyte','r'); % 训练集
fread(fid,16,'uint8');  % 前十六位是说明信息
for i = 1:60000
    tmp = fread(fid,28*28,'uint8=>double');
    tmp = reshape(tmp,28,28);
    tmp = tmp'; 
    train_img = [train_img,tmp(:)];
end
fclose(fid);

fid = fopen('train-labels.idx1-ubyte','r'); % 训练集标签
fread(fid,8,'uint8');
for i = 1:60000
   train_label(i) = fread(fid,1,'uint8=>double');
end
fclose(fid);

fid = fopen('t10k-images.idx3-ubyte','r'); %测试集
fread(fid,16,'uint8'); 
for i = 1:10000
    tmp = fread(fid,28*28,'uint8=>double');
    tmp = reshape(tmp,28,28);
    tmp = tmp'; 
    test_img = [test_img,tmp(:)];
end
fclose(fid);

fid = fopen('t10k-labels.idx1-ubyte','r'); % 测试集标签
fread(fid,8,'uint8');
for i = 1:10000
    test_label(i) = fread(fid,1,'uint8=>double');
end
fclose(fid);