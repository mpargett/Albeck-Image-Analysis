%IMAN_GETFRAME
%
%
%

%FIXME: Write header



function im = iman_getframe(dao, frm)
%Version check provision
if strcmpi(dao,'version'); im = 'v1.0'; return; end

%Index order
iid = struct('t',1, 'c', 2, 'xy', 3, 'z', 4);

% --- ----- --- Data Access Procedures  --- ----- --- 
switch dao.type
    case 'nd2'    
        %% For an ND2 file
        %IF XY position is stored as Series, set proper series and rinfo
        if dao.info.xySeries
            dao.read.setSeries(frm(iid.xy) - 1);
            rii = dao.rinfo(frm(iid.xy));
            %IF XY stored, in Z, replace any Z value with XY
        else iid.z = iid.xy; dao.read.setSeries(0); rii = dao.rinfo(1); 
        end
        
        %Assert order of indices via index structure (or hardcode?)
        im = double( bfGetPlane( dao.read, ...
            rii.lblr(frm(iid.z),frm(iid.t),frm(iid.c)) ) ); %z,t,c
        
    case 'image'  
        %% For an image file
        %Assert order of indices
        oi = cellfun(@(x)iid.(x), fieldnames(dao.info.imax));
        if isempty(oi); oi = []; end  %Fixes size of empties
        
        if dao.info.multipage
            %Determine file index to desired frame
            mpdn = fieldnames(dao.info.mpdim);
            %Assemble MultiPage Index from frame indices
            n = 0;
            for s = 1:numel(mpdn)
                n = (frm(iid.(mpdn{s}))-1) + n*dao.info.mpdim.(mpdn{s});
            end;    n = n + 1;
            %Read image data from MultiPage image (name assembled)
            im = double( imread(dao.read(frm(oi)), dao.info.ftype, ...
                'Index', n, 'Info', dao.info.imfi{frm(oi)}) );
        else
            if dao.info.multichan
                %Do we want to ever use multi-channel (RGB)?
            else
                %Load image data (assembling image name)
                im = double( imread(dao.read(frm(oi)), dao.info.ftype) );
            end
        end
        
end


end



