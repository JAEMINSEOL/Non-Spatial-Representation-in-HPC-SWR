ROOT.SaveRip = [ROOT.Processed '\ripples_mat\R0'];
ROOT.Rip0 = [ROOT.SaveRip];

Rip_CA1 = readtable([ROOT.Rip0 '\RipplesList_CA1' '.xlsx']);
Rip_SUB = readtable([ROOT.Rip0 '\RipplesList_SUB' '.xlsx']);

Rip_SUB.RippleDuration = Rip_SUB.EDtime - Rip_SUB.STtime;
Rip_CA1.RippleDuration = Rip_CA1.EDtime - Rip_CA1.STtime;

Rip_SUB(Rip_SUB.RippleDuration<0.04,:)=[];
Rip_CA1(Rip_CA1.RippleDuration<0.04,:)=[];

R1 = Rip_CA1;
R0 = Rip_SUB;
%%
for rc = 1:size(R1,1)
    thisR = R0(R0.rat==R1.rat(rc) & R0.session==R1.session(rc),:);
    st = R1.STtime(rc);
    ed = R1.EDtime(rc);

    if ~isempty(thisR)
        R1.Overlap(rc) = 0;
    else
        R1.Overlap(rc) = nan;
    end

    for rs=1:size(thisR,1)
        s = thisR.STtime(rs);
        e = thisR.EDtime(rs);

        if (st<=e && ed>=e) || (st<=s && ed>=s)
            R1.Overlap(rc) = R1.Overlap(rc)+1;
        end
    end
end

sum(R1.Overlap>0)


for rs = 1:size(R0,1)
    thisR = R1(R1.rat==R0.rat(rs) & R1.session==R0.session(rs),:);
    st = R0.STtime(rs);
    ed = R0.EDtime(rs);

    if ~isempty(thisR)
        R0.Overlap(rs) = 0;
    else
        R0.Overlap(rs) = nan;
    end

    for rc=1:size(thisR,1)
        s = thisR.STtime(rc);
        e = thisR.EDtime(rc);

        if (st<=e && ed>=e) || (st<=s && ed>=s)
            R0.Overlap(rs) = R0.Overlap(rs)+1;
        end
    end
end

sum(R0.Overlap>0)

%%

writetable(R0, [ROOT.SaveRip '\RipplesList_SUB' '.xlsx'],'writemode','replacefile')
writetable(R1, [ROOT.SaveRip '\RipplesList_CA1' '.xlsx'],'writemode','replacefile')

R0 = R0(~(R0.Overlap==0),:);
R1 = R1(~(R1.Overlap==0),:);

writetable(R0, [ROOT.SaveRip '\RipplesList_SUB_overlap' '.xlsx'],'writemode','replacefile')
writetable(R1, [ROOT.SaveRip '\RipplesList_CA1_overlap' '.xlsx'],'writemode','replacefile')

R0 = R0(R0.ensemble>2,:);
R1 = R1(R1.ensemble>2,:);