void bond_potyondy(Bond *bond, Part parti, Part partj, double gamma, double dt)
{
    double *omi, *omj;
    double *posi, *posj;
    double *posi_in, *posj_in;
    double kn, ks, radius, area;

    /*Getting bond values*/

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

//    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

    /*Getting particles positions and velocities*/

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

                perm += e*(omj[j]*(x_cont[k] - posj[k]) - omi[j]*(x_cont[k] - posi[k]));
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

    if (abs_S < 1.e-30)
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

    if (abs_T < 1.e-30)
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

    /*Computing forces and torques*/

    double fn[3], fs[3], ts[3];

    fn[0] = kn*area*(abs_rc - abs_r0)*r_c[0]/abs_rc;
    fn[1] = kn*area*(abs_rc - abs_r0)*r_c[1]/abs_rc;
    fn[2] = kn*area*(abs_rc - abs_r0)*r_c[2]/abs_rc;


    fs[0] = ks*area*abs_r0*theta*s[0];
    fs[1] = ks*area*abs_r0*theta*s[1];
    fs[2] = ks*area*abs_r0*theta*s[2];

    double abs_fs = abs(fs);

    ts[0] = 0.5*abs_rc*abs_fs*t[0];
    ts[1] = 0.5*abs_rc*abs_fs*t[1];
    ts[2] = 0.5*abs_rc*abs_fs*t[2];

    bond -> setForce_n(fn);
    bond -> setForce_s(fs);
    bond -> setTorque_s(ts);
}
