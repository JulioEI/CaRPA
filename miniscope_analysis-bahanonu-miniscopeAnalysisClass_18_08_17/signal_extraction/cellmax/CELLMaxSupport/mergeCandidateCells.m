function [cellImage, scaledProbability] = mergeCandidateCells(cellImage1, cellImage2, scaledProbability1, scaledProbability2)

cellImage = cellImage1 + cellImage2;
cellImage = cellImage/sum(cellImage(:));

scaledProbability = scaledProbability1 + scaledProbability2;