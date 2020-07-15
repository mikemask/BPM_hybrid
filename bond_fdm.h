void bond_fdm(Bond *bond, Part parti, Part partj)
{
    double *pos_i, *pos_j;
    double *pos0_i, *pos0_j;
    double *q_i, *q_j;
    double kn, ks, kb, kt;

//Getting bonds values

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();
//    kb = bond -> getStiffBend();
//    kt = bond -> getStiffTor();

//Getting particles values

    pos_i = parti.getPos();
    pos_j = partj.getPos();

    pos0_i = parti.getPosIn();
    pos0_j = partj.getPosIn();

    q_i = parti.getQuat();
    q_j = partj.getQuat();

    double *q_inv_i, *q_inv_j;

    q_inv_i = inverse(q_i);
    q_inv_j = inverse(q_j);

    /*Relative translation*/

    double r_f[3], r_0[3];

    r_f[0] = pos_i[0] - pos_j[0];
    r_f[1] = pos_i[1] - pos_j[1];
    r_f[2] = pos_i[2] - pos_j[2];

    r_0[0] = pos0_i[0] - pos0_j[0];
    r_0[1] = pos0_i[1] - pos0_j[1];
    r_0[2] = pos0_i[2] - pos0_j[2];

    double q_rf[4], *q_rc, r_c[3];

    q_rf[0] = 0.;
    q_rf[1] = r_f[0];
    q_rf[2] = r_f[1];
    q_rf[3] = r_f[2];

    q_rc = rotation(q_inv_j, q_rf, q_j, 3);

    r_c[0] = q_rc[1];
    r_c[1] = q_rc[2];
    r_c[2] = q_rc[3];

    double abs_rc = abs(r_c);
    double abs_r0 = abs(r_0);


    /*Normal and shear force*/

    double f_n[3], f_st[3], gamma;

    f_n[0] = kn*(abs_rc - abs_r0)/abs_rc*r_c[0];
    f_n[1] = kn*(abs_rc - abs_r0)/abs_rc*r_c[1];
    f_n[2] = kn*(abs_rc - abs_r0)/abs_rc*r_c[2];

    bond -> setForce_n(f_n);

    gamma = acos((r_0[0]*r_c[0]+r_0[1]*r_c[1]+r_0[2]*r_c[2])/(abs_rc*abs_r0));

    double S[3]={};

    S[0] = r_c[1]*(r_c[0]*r_0[1]-r_0[0]*r_c[1]) - r_c[2]*(r_0[0]*r_c[2]-r_c[0]*r_0[2]);
    S[1] = r_c[2]*(r_c[1]*r_0[2]-r_0[1]*r_c[2]) - r_c[0]*(r_0[1]*r_c[0]-r_c[1]*r_0[0]);
    S[2] = r_c[0]*(r_c[2]*r_0[0]-r_0[2]*r_c[0]) - r_c[1]*(r_0[2]*r_c[1]-r_c[2]*r_0[1]);

    double abs_S = abs(S);

    if (abs_S < 1.e-30)
    {
        f_st[0] = 0.;
        f_st[1] = 0.;
        f_st[2] = 0.;
    }

    else
    {
        double s[3];

        s[0] = S[0]/abs_S;
        s[1] = S[1]/abs_S;
        s[2] = S[2]/abs_S;

        f_st[0] = ks*abs_r0*gamma*s[0];
        f_st[1] = ks*abs_r0*gamma*s[1];
        f_st[2] = ks*abs_r0*gamma*s[2];
    }

    bond -> setForce_st(f_st);


    /*Torque due to shear force from relative translation*/

    double T[3];

    T[0] = r_0[1]*r_c[2] - r_0[2]*r_c[1];
    T[1] = r_0[2]*r_c[0] - r_0[0]*r_c[2];
    T[2] = r_0[0]*r_c[1] - r_0[1]*r_c[0];

    double abs_T = abs(T);
    double abs_fst = abs(f_st);

    double t[3], t_st[3];

    if (abs_T < 1.e-30)
    {
        t_st[0] = 0.;
        t_st[1] = 0.;
        t_st[2] = 0.;
    }

    else
    {
        t[0] = T[0]/abs_T;
        t[1] = T[1]/abs_T;
        t[2] = T[2]/abs_T;

        t_st[0] = 0.5*abs_rc*abs_fst*t[0];
        t_st[1] = 0.5*abs_rc*abs_fst*t[1];
        t_st[2] = 0.5*abs_rc*abs_fst*t[2];
    }


    bond -> setTorque_st(t_st);


//    /*Relative rotation*/
//
//    double *g_0, h[4];
//    bool flag = false;
//
//    g_0 = rotation(q_j, q_inv_j, q_i, 2);
//
//    if (r_0[0] < 1.e-30 && r_0[1] < 1.e-30)
//    {
//        h[0] = (sqrt(2)/2.)*sqrt((abs_r0 + r_0[2])/abs_r0);
//        h[1] = 0.;
//        h[2] = 0.;
//        h[3] = 0.;
//
//        flag = true;
//    }
//
//    else
//    {
//        h[0] = (sqrt(2)/2.)*sqrt((abs_r0 + r_0[2])/abs_r0);
//        h[1] = -(sqrt(2)/2.)*sqrt((abs_r0 - r_0[2])/abs_r0)*(r_0[1]/sqrt(r_0[0]*r_0[0]+r_0[1]*r_0[1]));
//        h[2] = (sqrt(2)/2.)*sqrt((abs_r0 - r_0[2])/abs_r0)*(r_0[0]/sqrt(r_0[0]*r_0[0]+r_0[1]*r_0[1]));
//        h[3] = 0.;
//    }
//
//
//    double *g, *h_inv, sinphi, cosphi, theta, psi;
//
//    h_inv = inverse(h);
//    g = rotation(h_inv, g_0, h, 3);
//
//    if(flag)
//    {
//        sinphi = 0.;
//        cosphi = 1.;
//    }
//
//    else
//    {
//        sinphi = (g[2]*g[3] - g[0]*g[1])/(sqrt((g[0]*g[0]+g[3]*g[3])*(g[1]*g[1]+g[2]*g[2])));
//        cosphi = (g[1]*g[3] + g[0]*g[2])/(sqrt((g[0]*g[0]+g[3]*g[3])*(g[1]*g[1]+g[2]*g[2])));
//    }
//
//
//    theta = acos(g[0]*g[0] - g[1]*g[1] - g[2]*g[2] + g[3]*g[3]);
//
//    psi = 2.*acos(g[0]/sqrt(g[0]*g[0] + g[3]*g[3]));
//
//    double g_prime_l[4], g_prime_r[4], *g_prime;
//
//    g_prime_l[0] = cos(theta/4.);
//    g_prime_l[1] = -sin(theta/4.)*sinphi;
//    g_prime_l[2] = sin(theta/4.)*cosphi;
//    g_prime_l[3] = 0.;
//
//    g_prime_r[0] = cos(psi/4.);
//    g_prime_r[1] = 0.;
//    g_prime_r[2] = 0.;
//    g_prime_r[3] = sin(psi/4.);
//
//    g_prime = rotation(q_i, g_prime_l, g_prime_r, 2);
//
//
//    /*Torques and forces relative to auxiliar coordinate frame*/
//
//    double tau_b_dp[4], tau_t_dp[4];
//
//    tau_b_dp[0] = 0.;
//    tau_b_dp[1] = -kb*theta*sinphi;
//    tau_b_dp[2] = kb*theta*cosphi;
//    tau_b_dp[3] = 0.;
//
//    tau_t_dp[0] = 0.;
//    tau_t_dp[1] = 0.;
//    tau_t_dp[2] = 0.;
//    tau_t_dp[3] = kt*psi;
//
//    double *tau_b_p, *tau_t_p, *g_prime_inv;
//
//    g_prime_inv = inverse(g_prime);
//
//    tau_b_p = rotation(g_prime, tau_b_dp, g_prime_inv, 3);
//    tau_t_p = rotation(g_prime, tau_t_dp, g_prime_inv, 3);
//
//    double f_s_p[4], t_s_p[4];
//
//    f_s_p[0] = 0.;
//    f_s_p[1] = -ks*abs_rc*theta*cosphi*0.5;
//    f_s_p[2] = -ks*abs_rc*theta*sinphi*0.5;
//    f_s_p[3] = 0.;
//
//    t_s_p[0] = 0.;
//    t_s_p[1] = ks*abs_rc*abs_rc*theta*sinphi*0.25;
//    t_s_p[2] = -ks*abs_rc*abs_rc*theta*cosphi*0.25;
//    t_s_p[3] = 0.;
//
//
//    /*Forces and torques relative to global coordinate system*/
//
//    double *t_br, *t_tr, *f_sr_step, *f_sr, *t_sr_step, *t_sr;
//    double g_dp[4], *g_dp_inv;
//
//    g_dp[0] = cos(gamma/2.);
//    g_dp[1] = t[0]*sin(gamma/2.);
//    g_dp[2] = t[1]*sin(gamma/2.);
//    g_dp[3] = t[2]*sin(gamma/2.);
//
//    g_dp_inv = inverse(g_dp);
//
//    t_br = rotation(h, tau_b_p, h_inv, 3);
//    t_tr = rotation(h, tau_t_p, h_inv, 3);
//
//    f_sr_step = rotation(h, f_s_p, h_inv, 3);
//    f_sr = rotation(g_dp, f_sr_step, g_dp_inv, 3);
//
//    t_sr_step = rotation(h, t_s_p, h_inv, 3);
//    t_sr = rotation(g_dp, t_sr_step, g_dp_inv, 3);

//    bond -> setTorque_br(t_br);
//    bond -> setTorque_tr(t_tr);
//    bond -> setForce_sr(f_sr);
//    bond -> setTorque_sr(t_sr);
}

