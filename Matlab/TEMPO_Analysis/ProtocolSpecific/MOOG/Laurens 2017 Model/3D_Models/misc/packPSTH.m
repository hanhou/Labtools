function data = packPSTH(psth)
    dim = size(psth);
    np = psth(1,1,:);
    sp = psth(1,dim(2),:);
    ev = psth(1:dim(1)-1,2:dim(2)-1,:);
    data = [np(:); sp(:); ev(:)];
end