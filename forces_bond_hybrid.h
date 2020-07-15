void bond_potyondy(Bond *bond, Part parti, Part partj, double dt)
{
    double *omi, *omj;
    double *posi, *posj;
    double *posi_in, *posj_in;
    double kn, ks, radius, *inertia, area;
//    double *tn, *tsr, *ttr;

    /*Getting bond values*/

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

//    tn = bond -> getTorque_n();
//    tsr = bond -> getTorque_sr();
//    ttr = bond -> getTorque_tr();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    omi = parti.getRot();
    omj = partj.getRot();


    /*Computing relative motion and contact plane versors*/

    double r_f[3], r_0[3], r_c[3];

    r_f[0] = posi[0] - posj[0];
    r_f[1] = posi[1] - posj[1];
    r_f[2] = posi[2] - posj[2];

    r_0[0] = posi_in[0] - posj_in[0];
    r_0[1] = posi_in[1] - posj_in[1];
    r_0[2] = posi_in[2] - posj_in[2];

    double x_cont[3];

    x_cont[0] = 0.5*(posi[0] + posj[0]);
    x_cont[1] = 0.5*(posi[1] + posj[1]);
    x_cont[2] = 0.5*(posi[2] + posj[2]);

    int e;

    for (int i=0; i<3; i++)
    {
        double perm = 0.;

        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                if (i==j || i==k || j==k || i==j==k)
                    e = 0;
                else if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
                    e = 1;
                else
                    e = -1;

                perm += e*(omi[j]*(x_cont[k] - posi[k]) - omj[j]*(x_cont[k] - posj[k]));
            }
        }

        r_c[i] = r_f[i] + perm*dt;
    }

    double abs_rc, abs_r0, theta, abs_S, abs_T;
    double S[3], s[3], T[3], t[3];

    abs_r0 = abs(r_0);
    abs_rc = abs(r_c);

    S[0] = r_c[1]*(r_c[0]*r_0[1]-r_0[0]*r_c[1]) - r_c[2]*(r_0[0]*r_c[2]-r_c[0]*r_0[2]);
    S[1] = r_c[2]*(r_c[1]*r_0[2]-r_0[1]*r_c[2]) - r_c[0]*(r_0[1]*r_c[0]-r_c[1]*r_0[0]);
    S[2] = r_c[0]*(r_c[2]*r_0[0]-r_0[2]*r_c[0]) - r_c[1]*(r_0[2]*r_c[1]-r_c[2]*r_0[1]);

    abs_S = abs(S);

    if (abs_S < 1.e-50)
    {
        s[0] = 0.;
        s[1] = 0.;
        s[2] = 0.;
    }

    else
    {
        s[0] = S[0]/abs_S;
        s[1] = S[1]/abs_S;
        s[2] = S[2]/abs_S;
    }

    theta = acos((r_0[0]*r_c[0]+r_0[1]*r_c[1]+r_0[2]*r_c[2])/(abs_rc*abs_r0));

    T[0] = r_0[1]*r_c[2] - r_0[2]*r_c[1];
    T[1] = r_0[2]*r_c[0] - r_0[0]*r_c[2];
    T[2] = r_0[0]*r_c[1] - r_0[1]*r_c[0];

    abs_T = abs(T);

    if (abs_T < 1.e-50)
    {
        t[0] = 0.;
        t[1] = 0.;
        t[2] = 0.;
    }

    else
    {
        t[0] = T[0]/abs_T;
        t[1] = T[1]/abs_T;
        t[2] = T[2]/abs_T;
    }

    double delta_theta[3];

    delta_theta[0] = (omi[0] - omj[0])*dt;
    delta_theta[1] = (omi[1] - omj[1])*dt;
    delta_theta[2] = (omi[2] - omj[2])*dt;

    double delta_theta_n, delta_theta_t, delta_theta_s;

    delta_theta_n = (delta_theta[0]*r_c[0] + delta_theta[1]*r_c[1] + delta_theta[2]*r_c[2])/abs_rc;
    delta_theta_t = delta_theta[0]*t[0] + delta_theta[1]*t[1] + delta_theta[2]*t[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];

    /*Computing forces and torques*/

    double fn[3], fs[3], tst[3];
    double tn[3], tsr[3], ttr[3];

    fn[0] = kn*area*(abs_rc - abs_r0)*r_c[0]/abs_rc;
    fn[1] = kn*area*(abs_rc - abs_r0)*r_c[1]/abs_rc;
    fn[2] = kn*area*(abs_rc - abs_r0)*r_c[2]/abs_rc;

    fs[0] = ks*area*abs_r0*sin(theta)*s[0];
    fs[1] = ks*area*abs_r0*sin(theta)*s[1];
    fs[2] = ks*area*abs_r0*sin(theta)*s[2];

    double abs_fs = abs(fs);

    tn[0] = ks*inertia[1]*delta_theta_n*r_c[0]/abs_rc;
    tn[1] = ks*inertia[1]*delta_theta_n*r_c[1]/abs_rc;
    tn[2] = ks*inertia[1]*delta_theta_n*r_c[2]/abs_rc;

    tst[0] = 0.5*abs_rc*abs_fs*t[0];
    tst[1] = 0.5*abs_rc*abs_fs*t[1];
    tst[2] = 0.5*abs_rc*abs_fs*t[2];

    ttr[0] = kn*inertia[0]*delta_theta_t*t[0];
    ttr[1] = kn*inertia[0]*delta_theta_t*t[1];
    ttr[2] = kn*inertia[0]*delta_theta_t*t[2];

    tsr[0] = kn*inertia[0]*delta_theta_s*s[0];
    tsr[1] = kn*inertia[0]*delta_theta_s*s[1];
    tsr[2] = kn*inertia[0]*delta_theta_s*s[2];

    bond -> setForce_n(fn);
    bond -> setForce_s(fs);
    bond -> setTorque_st(tst);
    bond -> setTorque_n(tn);
    bond -> setTorque_sr(tsr);
    bond -> setTorque_tr(ttr);
}
