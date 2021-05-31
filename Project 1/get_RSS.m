function [RSS] = get_RSS(measured,sens,Pt,dev_stand)

dist = (vecnorm((measured-sens)'))';
d = (dist<=8);
RSS = d.*(Pt-40.2-20.*log10(dist)+dev_stand*randn()) + (1-d).*(Pt-58.5-33.*log10(dist)+dev_stand*randn());

end