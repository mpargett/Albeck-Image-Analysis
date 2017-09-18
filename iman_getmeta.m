%IMAN_GETMETA
%   Collect MetaData associated with a Data Access Object


function GMD = iman_getmeta(dao, bkmd, op)
%Version check provision
if strcmpi(dao,'version'); GMD = 'v1.0'; return; end

%Search for MetaData, depending on data type
switch lower(dao.type)
    case 'nd2'
        GMD = nd2_meta(dao.read, dao.rinfo, bkmd, op);
    case 'image'
        GMD = struct();  %Initialize, but leave be fill with backup
        %FIXME: Assert that bkmd is complete, prompt with specifics needed?
end

%Enact override of file metadata, if applicable
%   Enforce backup values as requested
if ~isempty(op.mdover);
    for s = 1:size(op.mdover,1)
        GMD.(op.mdover{s,1}).(op.mdover{s,2}) = ...
            bkmd.(op.mdover{s,1}).(op.mdover{s,2});
    end
end

%Ensure all backup field are included, if not found in MetaData search
GMD = subf_rfill(bkmd, GMD);

end


%Recursive filling
function [dst] = subf_rfill(src, dst)
%Get fieldnames for source and destination structures
fns = fieldnames(src)';  
%   IF destination is not a structure, use empty fieldnames
if isstruct(dst); fnd = fieldnames(dst)'; else fnd = {}; end
%Get indices of fields missing from destination
fnmiss = ~ismember(fns, fnd);

%Fill missing entries from source
for s = fns(fnmiss);  dst.(s{1}) = src.(s{1});  end

%Recursively search for missing entries within matching sub-structures
for s = fns(~fnmiss)
    if isstruct(src.(s{1}))
        dst.(s{1}) = subf_rfill(src.(s{1}), dst.(s{1})); 
    end
end


end


