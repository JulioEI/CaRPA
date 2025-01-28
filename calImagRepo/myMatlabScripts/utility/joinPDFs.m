function joinPDFs(outputPDF,inputPDF)

%Interface for using the append_pdfs toolbox. Calling the function with no
%arguments or only first argument brings up a gui. Inputing a folder string as
%the inputPDF will join all PDFs in the folder. Otherwise inputPDF, if
%entered, should be a cell of strings (with paths included). OutputPF, if
%entered, should be the save string, with path included.

if nargin < 2
    [file,path] = uigetfile('*.pdf','MultiSelect','on');
    if ~iscell(file)
        file = {file};
    end
    inputPDF = cellfun(@(x) [path,filesep,x],file,'UniformOutput',0);
    
    if nargin < 1
        [fileOut,pathOut] = uiputfile('*.pdf');
        outputPDF = [pathOut,filesep,fileOut];
    end
end

if isstr(inputPDF) %the input is a folder, join all files
    pdfFiles = dir([inputPDF,filesep,'*.pdf']);
    inputPDF = cellfun(@(x) [pdfFiles(1).folder,'\',x],{pdfFiles.name},'UniformOutput',0);
end

append_pdfs(outputPDF,inputPDF{:},'-f')
    
end

