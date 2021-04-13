function varargout = immontage(imgCellArray, gridSize, border, boxLen, fcnHandle, padVal)
%IMMONTAGE Display multiple image frames as rectangular montage.
% Syntax:
%   out = immontage(imgCellArray, gridSize, border, boxLen, fcnHandle, padVal)
%
%   IMG_ARRAY can be either cell array or 4D array

if ~iscell(imgCellArray)
    imgArray = imgCellArray;
    numImg = size(imgArray, 4);
    imgCellArray = cell(1, numImg);
    for iImg = 1:numImg
        imgCellArray{iImg} = imgArray(:,:,:,iImg);
    end
end
if ~exist('gridSize', 'var') || isempty(gridSize)
    [numRow, numCol] = size(imgCellArray);
else
    numRow = gridSize(1);
    numCol = gridSize(2);
end
if ~exist('border', 'var') || isempty(border)
    border = 0;
end
if ~exist('boxLen', 'var') || isempty(boxLen)
    boxLen = 128;
end
if ~exist('fcnHandle', 'var') || isempty(fcnHandle)
    fcnHandle = @(x) x;
end
if ~exist('padVal', 'var') || isempty(padVal)
    padVal = 255;
end
imgCellArray = cellfun(fcnHandle, imgCellArray, 'UniformOutput', false);
imgClass = class(imgCellArray{1});
imgOut = padVal*ones(numRow*boxLen, numCol*boxLen, 3, imgClass);
numImg = numel(imgCellArray);
boxLen_ = boxLen-2*border;
for iRow = 1:numRow
    for iCol = 1:numCol
        ii = sub2ind([numRow numCol], iRow, iCol);
        if ii > numImg
            break
        end
        img = imgCellArray{ii};
        [height, width, nChnl] = size(img);
        if height < width
            img = imresize(img, boxLen_/width);
            [height, width, ~] = size(img);
            offset_h = fix((boxLen_-height)/2);
            offset_w = 0;
        elseif height > width
            img = imresize(img, boxLen_/height);
            [height, width, ~] = size(img);
            offset_w = fix((boxLen_-width)/2);
            offset_h = 0;
        else
            img = imresize(img, boxLen_/width);
            height = boxLen_;
            width = boxLen_;
            offset_w = 0;
            offset_h = 0;
        end
        if nChnl == 1
            img = repmat(img, [1 1 3]);
        end
        imgOut((iRow-1)*boxLen+1+offset_h+border:(iRow-1)*boxLen+height+offset_h+border, ...
            (iCol-1)*boxLen+1+offset_w+border:(iCol-1)*boxLen+width+offset_w+border,:) = img;
    end
end
if nargout > 0
    varargout{1} = imgOut;
else
    imshow(imgOut);
end
