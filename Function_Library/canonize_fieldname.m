function fno = canonize_fieldname(fni, varargin)
%This function attempts to create properly formatted field names. It's not
%yet complete, if I run into newer surprising characters or character
%combinations I will add them here.
%
%Usage:
%  OutputFieldName = canonize_fieldname(InputFieldName, ...);
%
%Input:
%  InputFieldName: string, field name to clean up
%
%Optional input argument
%  'ReplaceMap': Nx2 cellstr, describes replacement of invalid characters.
%    N is the number of characters to take care of, first column is what to
%    replace, second column is what to use instead. Default:
%      ReplaceMap = {' ' '_'; '-' '_'; '(' ''; ')' ''};
%  'PreChar': char, if InputFieldName does not begin with a letter, this
%    character is added to prepend the name. Default: 'a'
%Output:
%  OutputFieldName: string, hopefully valid field name
%
%Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parameters
repmap = {' ' '_'; '-' '_'; '(' ''; ')' ''; '.' '_'};

p = inputParser;
addRequired(p, 'fni', @isstr)
addParamValue(p, 'ReplaceMap', repmap, @iscell); %#ok<*NVREPL>
addParamValue(p, 'PreChar', 'a', @isstr);
parse(p, fni, varargin{:});

repmap = p.Results.ReplaceMap;

%% Main
fno = fni;
for ridx = 1:size(repmap, 1)
    fno = strrep(fno, repmap{ridx,1}, repmap{ridx,2});
end

if isstrprop(fno(1), 'digit')
    fno = [p.Results.PreChar fno];
end

end

