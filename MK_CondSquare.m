function wh = MK_CondSquare(var1,var2)
wh(1,1) = sum(var1 & var2);
wh(1,2) = sum(var1 & ~var2);
wh(2,1) = sum(~var1 & var2);
wh(2,2) = sum(~var1 & ~var2);