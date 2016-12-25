function d=sine_noah(out,x,y,frqf)

d=abs(y-(out(1)*sin(frqf*2*pi*x+out(2)/180*pi)+out(3)));
return
