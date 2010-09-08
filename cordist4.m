function D = cordist4(z1, z2)
		    
D = cordist2(z1, z2);
cog1 = cors_cog(z1);
cog2 = cors_cog(z2);

dx = abs(cog1(1)-cog2(1)); 
dp = abs(angle(exp(i*(cog1(2)-cog2(2)))));
dy = abs(cog1(3)-cog2(3));
dc = abs(cog1(4)-cog2(4));
w = min(cog1(4),cog2(4))/max(cog1(4),cog2(4));
D = [D; [dx; dy; dp]*w; dc*(1-w)];
