%3-2-1 Transformaiton

syms th ph ps

r_ph=[cos(ph) sin(ph) 0; -sin(ph) cos(ph) 0; 0 0 1;];
r_th=[cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th);];
r_ps=[1 0 0; 0 cos(ps) sin(ps); 0 -sin(ps) cos(ps);];

%multiply in order
R=r_ps*r_th*r_ph;

%substititutions
r_ph_sub = subs(r_ph, ph, 0.4);
r_th_sub = subs(r_th, th, 0.6);
r_ps_sub = subs(r_ps, ps, pi/2);

n=4;

r_ph_num = vpa(r_ph_sub,n);
r_th_num = vpa(r_th_sub,n);
r_ps_num = vpa(r_ps_sub,n);

R_num=r_ps_num*r_th_num*r_ph_num