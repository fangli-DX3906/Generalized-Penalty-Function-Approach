function IRFStruct = mkIRFStruct(pemodel, IRFmat, uci, lci, alpha)

irHoriz = size(IRFmat, 2);
nvar = size(IRFmat, 1);
model = pemodel.basicModel;

IRFStruct.irf = zeros(irHoriz, nvar, nvar);

for i = 1:nvar
    for j = 1:irHoriz
        IRFStruct.irf(j, :, i) = IRFmat(i, j, :);
    end
end

info.h = irHoriz;
info.whos = model.SeriesNames;
info.model = model.Description;
info.alpha = alpha;
IRFStruct.uci = uci;
IRFStruct.lci = lci;
IRFStruct.info = info;

end

