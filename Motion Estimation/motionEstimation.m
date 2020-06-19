function result = motionEstimation(frame178, frame179, QP)



[M,N] = size(frame178);
quant = zeros(M,N); % array that holds quant operators of hole image
result = zeros(M,N);


%split frame 178 into 4x4 blocks
numFullBlocksX = floor(N/4); 
numFullBlocksY = floor(M/4);
xBlocks = [repmat(4,numFullBlocksX,1); mod(N,4)*ones(mod(N,4)>0)];
yBlocks = [repmat(4,numFullBlocksY,1); mod(M,4)*ones(mod(M,4)>0)];

splited_frame178 = mat2cell(frame178,yBlocks,xBlocks);%its the array that contains 4x4 blocks of frame178
int_DCT = mat2cell(zeros(M,N),yBlocks,xBlocks); %array that hold integer dct of each block 
block_quants = mat2cell(zeros(M,N),yBlocks,xBlocks); %array that hold quantization of each block
inverse_quant_array = mat2cell(zeros(M,N),yBlocks,xBlocks); 
inverse_integer_array = mat2cell(zeros(M,N),yBlocks,xBlocks);

[X,Y] = size(splited_frame178);




c1 = 1;
c2 = 1;
c3 = 4;
c4 = 4;

%calculate here the integer DCT of each block
for i = 1:X
	for j = 1:Y
		int_DCT{i,j} = integer_transform(splited_frame178{i,j});
	end
	
end

%calculate here the quantization of each block after integer DCT
for i = 1:X
	for j = 1:Y
		block_quants{i,j} = quantization(int_DCT{i,j}, QP);
	end
	
end

%calculate here quant operators for hole image
for i = 1:X
	for j = 1:Y
		quant(c1:c3, c2:c4) = block_quants{i,j};
		if j < Y
			c2 = c2 + 4;
			c4 = c4 + 4;
		end
	end
	c1 = c1 + 4;
	c3 = c3 + 4;
	c2 = 1;
	c4 = 4;
end


%calculate entropy here
entropy(uint8(abs(quant)))

%do inverse quantization here
for i = 1:X
	for j = 1:Y
		inverse_quant_array{i,j} = inv_quantization(block_quants{i,j}, QP);
	end
	
end

%do inverse integer transform here
for i = 1:X
	for j = 1:Y
		inverse_integer_array{i,j} = inv_integer_transform(inverse_quant_array{i,j});
	end
	
end



c1 = 1;
c2 = 1;
c3 = 4;
c4 = 4;

%union the image here after the transforamtions and inverse transformation.
for i = 1:X
	for j = 1:Y
		result(c1:c3, c2:c4) = inverse_integer_array{i,j};
		if j < Y
			c2 = c2 + 4;
			c4 = c4 + 4;
		end
	end
	c1 = c1 + 4;
	c3 = c3 + 4;
	c2 = 1;
	c4 = 4;
end


result = round(result/64);
result = uint8(result);
frame178 = uint8(frame178);

figure, imshow(frame178)
title({"Original frame178"})
figure, imshow(result)
title({"recreated frame178 after inverse transformations"})

%calculate PSNR of recreated frame178 here

MSE = sum((frame178(:) - result(:)).^2)/prod(size(frame178));
PSNR = 10*log10(255*255/MSE);
disp('to PSNR tou anadhmiourghmenou frame178 einai:');
disp(PSNR);

%ERWTHMA 3
[M,N] = size(frame179);
frame179_prediction = zeros(M,N);
numFullBlocksX = floor(N/16); 
numFullBlocksY = floor(M/16);
xBlocks = [repmat(16,numFullBlocksX,1); mod(N,16)*ones(mod(N,16)>0)];
yBlocks = [repmat(16,numFullBlocksY,1); mod(M,16)*ones(mod(M,16)>0)];

splited_frame179 = mat2cell(frame179,yBlocks,xBlocks);
[X,Y] = size(splited_frame179);
mvx = zeros(X,Y);
mvy = zeros(X,Y);

%calculate here motion vectors using full search technique
for row = 1:16:M-16+1 %for each macroblock 16x16 in frame 179
	for column = 1:16:N-16+1
		min_SAD_value = 256*16*16; i = 0; j = 0;
		iblk = 0;jblk = 0;
		for k = -8:8 
			if row+k < 1 || row+k+16-1 > M %checking here for out of bounds
				continue
			end
			for l = -8:8
				if column+l < 1 || column+l+16-1 > N %checking here for out of bounds
					continue;
				
				end	
				SAD = 0;
				for a = 1:16 %for every pixel of mcroblock and block 
					for b = 1:16
						SAD = SAD + abs((frame179(row:row+16-1,column:column+16-1)(a,b) - double(result(row+k:row+k+16-1,column+l:column+l+16-1)(a,b))));
					end
				end
				if SAD < min_SAD_value
					min_SAD_value = SAD; i = l; j = k;
				end	
				
			end
		end
		%copy the best maching block into frame179 prediction now
		frame179_prediction(row:row+16-1,column:column+16-1) = result(row+j:row+j+16-1,column+i:column+i+16-1); 
		iblk = (floor)(row-1)/16 + 1; jblk = (floor)(column-1)/16 + 1;
		mvx(iblk,jblk)=i; mvy(iblk,jblk)=j;%storing motion vectors in 2 different arrays, 
		%that menas that each x of a motion vector will be saved in mvx and corresponding y in mvy
	end
end			

%disp(frame179_prediction);	

frame179_prediction = uint8(frame179_prediction);
figure, imshow(frame179_prediction);
title({"prediction of frame179 using frame178 and motion vectors"})

%ERWTHMA 4

prediction_fault = zeros(M,N);
quant = zeros(M,N);
result = zeros(M,N);
prediction_fault = frame179-frame179_prediction;
prediction_fault = double(prediction_fault);
%disp(prediction_fault);
%figure, imshow(prediction_fault)

numFullBlocksX = floor(N/4); 
numFullBlocksY = floor(M/4);
xBlocks = [repmat(4,numFullBlocksX,1); mod(N,4)*ones(mod(N,4)>0)];
yBlocks = [repmat(4,numFullBlocksY,1); mod(M,4)*ones(mod(M,4)>0)];

splited_prediction_fault = mat2cell(prediction_fault,yBlocks,xBlocks);
int_DCT_prediction_fault = mat2cell(zeros(M,N),yBlocks,xBlocks);
quantization_prediction_fault = mat2cell(zeros(M,N),yBlocks,xBlocks);
inverse_int_DCT_prediction_fault = mat2cell(zeros(M,N),yBlocks,xBlocks);
inverse_quantization_prediction_fault = mat2cell(zeros(M,N),yBlocks,xBlocks);

[X,Y] = size(splited_prediction_fault);

for i = 1:X
	for j = 1:Y
		int_DCT_prediction_fault{i,j} = integer_transform(splited_prediction_fault{i,j});
	end
	
end

for i = 1:X
	for j = 1:Y
		quantization_prediction_fault{i,j} = quantization(int_DCT_prediction_fault{i,j}, QP);
	end
	
end

c1 = 1;
c2 = 1;
c3 = 4;
c4 = 4;

for i = 1:X
	for j = 1:Y
		quant(c1:c3, c2:c4) = quantization_prediction_fault{i,j};
		if j < Y
			c2 = c2 + 4;
			c4 = c4 + 4;
		end
	end
	c1 = c1 + 4;
	c3 = c3 + 4;
	c2 = 1;
	c4 = 4;
end

entropy(uint8(abs(quant)))

for i = 1:X
	for j = 1:Y
		inverse_quantization_prediction_fault{i,j} = inv_quantization(quantization_prediction_fault{i,j}, QP);
	end
	
end

for i = 1:X
	for j = 1:Y
		inverse_int_DCT_prediction_fault{i,j} = inv_integer_transform(inverse_quantization_prediction_fault{i,j});
	end
	
end


c1 = 1;
c2 = 1;
c3 = 4;
c4 = 4;


for i = 1:X
	for j = 1:Y
		result(c1:c3, c2:c4) = inverse_int_DCT_prediction_fault{i,j};
		if j < Y
			c2 = c2 + 4;
			c4 = c4 + 4;
		end
	end
	c1 = c1 + 4;
	c3 = c3 + 4;
	c2 = 1;
	c4 = 4;
end

result = round(result/64);
result = uint8(result);
figure, imshow(result);
title({"recreated prediction fault after inverse transformations"})


%recreated_frame179 = zeros(M,N);

recreated_frame179 = result + frame179_prediction;

recreated_frame179 = uint8(recreated_frame179);
figure, imshow(recreated_frame179);
title({"recreated frame179"})
frame179  = uint8(frame179);
MSE = sum((frame179(:) - recreated_frame179(:)).^2)/prod(size(frame179));
PSNR = 10*log10(255*255/MSE);
disp('to PSNR tou anadhmiourghmenou frame179 einai:');
disp(PSNR);

end
