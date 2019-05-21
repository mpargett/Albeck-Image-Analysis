%XLSREAD Backwards Compatibility Overload
%   Usage mimics the original XLSREAD, but uses READCELL if MATLAB version
%   is 2019a or later.
%   
%   This compatibility wrapper neglects the use of a custom function on the
%   data in the spreadsheet.
%
%Built-in XLSREAD Help:
% XLSREAD Read Microsoft Excel spreadsheet file.
%   [NUM,TXT,RAW]=XLSREAD(FILE) reads data from the first worksheet in
%   the Microsoft Excel spreadsheet file named FILE and returns the numeric
%   data in array NUM. Optionally, returns the text fields in cell array
%   TXT, and the unprocessed data (numbers and text) in cell array RAW.
%
%   [NUM,TXT,RAW]=XLSREAD(FILE,SHEET) reads the specified worksheet.
%
%   [NUM,TXT,RAW]=XLSREAD(FILE,SHEET,RANGE) reads from the specified SHEET
%   and RANGE. Specify RANGE using the syntax 'C1:C2', where C1 and C2 are
%   opposing corners of the region. Not supported for XLS files in BASIC
%   mode.
%
%   [NUM,TXT,RAW]=XLSREAD(FILE,SHEET,RANGE,'basic') reads from the
%   spreadsheet in BASIC mode, the default on systems without Excel
%   for Windows. RANGE is supported for XLSX files only.
%
%   [NUM,TXT,RAW]=XLSREAD(FILE,RANGE) reads data from the specified RANGE
%   of the first worksheet in the file. Not supported for XLS files in
%   BASIC mode.
%
%   The following syntaxes are supported only on Windows systems with Excel
%   software:
%
%   [NUM,TXT,RAW]=XLSREAD(FILE,-1) opens an Excel window to select data
%   interactively.
%
%   Input Arguments:
%
%   FILE    String that specifies the name of the file to read.
%   SHEET   Worksheet to read. One of the following:
%           * String that contains the worksheet name.
%           * Positive, integer-valued scalar indicating the worksheet
%             index.
%   RANGE   String that specifies a rectangular portion of the worksheet to
%           read. Not case sensitive. Use Excel A1 reference style.
%           If you do not specify a SHEET, RANGE must include both corners
%           and a colon character (:), even for a single cell (such as
%           'D2:D2').
%   'basic' Flag to request reading in BASIC mode, which is the default for
%           systems without Excel for Windows.  In BASIC mode, XLSREAD:
%           * Reads XLS, XLSX, XLSM, XLTX, and XLTM files only.
%           * Does not support an xlRange input when reading XLS files.
%             In this case, use '' in place of xlRange.
%           * For XLS files, requires a string to specify the SHEET,
%             and the name is case sensitive.
%           * Does not support function handle inputs.
%           * Imports all dates as Excel serial date numbers. Excel
%             serial date numbers use a different reference date than
%             MATLAB date numbers.
%   -1      Flag to open an interactive Excel window for selecting data.
%           Select the worksheet, drag and drop the mouse over the range
%           you want, and click OK. Supported only on Windows systems with
%           Excel software.
%
%   Notes:
%
%   * On Windows systems with Excel software, XLSREAD reads any file
%     format recognized by your version of Excel, including XLS, XLSX,
%     XLSB, XLSM, and HTML-based formats.
%
%   * If your system does not have Excel for Windows, XLSREAD operates in
%     BASIC mode (see Input Arguments).
%
%   * XLSREAD imports formatted dates as strings (such as '10/31/96'),
%     except in BASIC mode. In BASIC mode, XLSREAD imports all dates as
%     serial date numbers. Serial date numbers in Excel use different
%     reference dates than date numbers in MATLAB. For information on
%     converting dates, see the documentation on importing spreadsheets.
%
%   Examples:
%
%   % Create data for use in the examples that follow:
%   values = {1, 2, 3 ; 4, 5, 'x' ; 7, 8, 9};
%   headers = {'First', 'Second', 'Third'};
%   xlswrite('myExample.xls', [headers; values]);
%   moreValues = rand(5);
%   xlswrite('myExample.xls', moreValues, 'MySheet');
%
%   % Read data from the first worksheet into a numeric array:
%   A = xlsread('myExample.xls')
%
%   % Read a specific range of data:
%   subsetA = xlsread('myExample.xls', 1, 'B2:C3')
%
%   % Read from a named worksheet:
%   B = xlsread('myExample.xls', 'MySheet')
%
%   % Request the numeric data, text, and a copy of the unprocessed (raw)
%   % data from the first worksheet:
%   [ndata, text, alldata] = xlsread('myExample.xls')
%
%   See also XLSWRITE, XLSFINFO, DLMREAD, IMPORTDATA, TEXTSCAN.

%   Copyright 1984-2015 The MathWorks, Inc.
%=============================================================================


function [ynum, ytxt, yraw] = xlsread(fp, varargin)
%Check MATLAB Version
vv = ver('MATLAB');                                 %Get version
rr = regexpi(vv.Release, 'R(\d*)\D?', 'tokens');    %Get release year
isnew = str2num(rr{1}{1}) >= 2019;                  %Check if 2019 or later

%Choose procedure based on Version
if isnew %IF Version 2019a or later
    %Parse input
    nin = numel(varargin);
    %Check if inputs looks like a RANGE
    isr = cellfun(@isempty, ...
        regexpi(varargin(1:min(2,nin)), '\D+\d+:\D+\d+') );
    pn = {'SHEET', 'RANGE'}; %Prepare parameter names
    
    %Use readcell, with SHEET, RANGE parameters matched up
    switch nin
        case 0;     yraw = readcell(fp);
        case 1;     yraw = readcell(fp, pn{1+isr}, varargin{1});
        otherwise;  yraw = readcell(fp, pn{1+isr(1)}, varargin{1}, ...
                                        pn{1+isr(2)}, varargin{2});
    end
    %Pre-define function to get ranges from logical matrix
    get_bnds = @(inp,dm)find(any(inp,dm), 1, 'first') : ...
                        find(any(inp,dm), 1, 'last');
    
    %   Get numeric only array
    chk = cellfun(@(x)isnumeric(x)&&~isnan(x), yraw);   %Check for numbers
    b1 = get_bnds(chk,2); b2 = get_bnds(chk,1);	%Get numeric ranges  
    ynum = nan(numel(b1),numel(b2));           	%Initialize numeric outputs
    ynum(chk(b1,b2)) = [yraw{chk}];             %Assign numeric output
    
    %   Get text only cell array
    chk = cellfun(@ischar, yraw);               %Check for text
    b1 = get_bnds(chk,2); b2 = get_bnds(chk,1);	%Get text ranges  
    ytxt = cell(numel(b1),numel(b2));         	%Initialize text outputs
    ytxt(chk(b1,b2)) = yraw(chk);               %Assign text output    
    
    
else %IF Version prior to 2019a
    %Pass everything to original XLSREAD
    %   Get path to original XLSREAD function (since this overloads it)
    xlspath = which('xlsread','-all');
    xlspath = xlspath{~cellfun(@isempty, regexpi(xlspath, '\\matlab\\R'))};
    xlspath = regexprep(xlspath, '\\xlsread.m', '\\');
    wd = cd; cd(xlspath);           %Store working directory, and change
    xlsread_handle = @xlsread;      %Create handle to built-in XLSREAD
    cd(wd);                         %Return to the working directory
    
    %Run XLSREAD via the direct handle
    [ynum, ytxt, yraw] = xlsread_handle(fp, varargin{:});
end


end