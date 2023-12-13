[~,ia,ic] = unique(UnitPair_field.CA1.UID1);

UF1 = UnitPair_field.CA1(ia,:);

[~,ia,ic] = unique(UnitPair_field.SUB.UID1);

UF0 = UnitPair_field.SUB(ia,:);


    UF0.Nsp_M = nanmax(abs([UF0.Nsp_L,UF0.Nsp_R,UF0.Nsp_C]),[],2);
    UF1.Nsp_M = nanmax(abs([UF1.Nsp_L,UF1.Nsp_R,UF1.Nsp_C]),[],2);
 
    Reg='CA1';
for u=1:size(UnitsTable_field.(Reg),1)
    idx = find(strcmp(UnitsTable_field.(Reg).ID(u),UnitPair_field.(Reg).UID1),1);
    if ~isempty(idx)
    UnitsTable_field.(Reg).RPR(u) = UnitPair_field.(Reg).p(idx) / UnitPair_field.(Reg).p0(idx);
    else
        UnitsTable_field.(Reg).RPR(u)=nan;
    end
end
%%
UF1 = UnitsTable_field.CA1;
UF0 = UnitsTable_field.SUB;

figure;

subplot(1,2,1)
x=abs(UF0.RDI_RScene); y=UF0.RPR;
x(~(y<=1))=[]; y(~(y<=1))=[];
y(isnan(x))=[]; x(isnan(x))=[];
scatter(x,y,40,CList(1,:),'filled')
ylim([0 1]); xlim([0 1.2])
X = [ones(size(x)) x];
B = X\y;
Rsq = 1 - sum((y - X*B).^2)/sum((y - mean(y)).^2);
title(num2str(Rsq))

subplot(1,2,2)
x=abs(UF1.RDI_RScene); y=UF1.RPR;
x(~(y<=1))=[]; y(~(y<=1))=[];
y(isnan(x))=[]; x(isnan(x))=[];
scatter(x,y,40,CList(2,:),'filled')
ylim([0 1]); xlim([0 1.2])


corrcoef(x,y)
X = [ones(size(x)) x];
B = X\y;
Rsq = 1 - sum((y - X*B).^2)/sum((y - mean(y)).^2);
title(num2str(Rsq))