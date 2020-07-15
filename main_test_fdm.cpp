#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include "quat_inv.h"
#include "abs.h"
#include "bond.h"
#include "part.h"
#include "distance.h"
#include "rotation.h"
#include "bond_fdm.h"

int main()
{
    double const rho = 7000.;
    double const d = 0.001;
    double const E = 2.7e6;
    double const nu = 0.4;
    double const dt = 0.0000001;
    double const eps = d/1000.;
    double const lambda = 1.1;
    double const gamma = 1.;
    int bonds = 0;
    int count = 0;
    int const n = 100000000;
    int const n_part = 2;

    std::ofstream out_file;
    out_file.open("data_two_rot.csv");
    out_file << "timestep,ID,x,y,z,u,v,w,fx,fy,fz,omx,omy,omz,tx,ty,tz\n";
    int const n_out = 100000;
    bool flag = false;

    Part vec_part[n_part];

    Bond vec_bond;

////Initial particle configuration

    double const m = rho*(4./3.)*3.14159265*pow((d/2.),3.);

    double pos1[3], pos2[3];

    pos1[0] = 0.;
    pos1[1] = 0.;
    pos1[2] = d+eps;

    pos2[0] = 0.;
    pos2[1] = 0.;
    pos2[2] = 0.;

    vec_part[0].setPos(pos1);
    vec_part[0].setPosIn(pos1);
    vec_part[1].setPos(pos2);
    vec_part[1].setPosIn(pos2);

    for (int i=0; i<n_part; i++)
    {
        vec_part[i].setMass(m);
        vec_part[i].setDiameter(d);
    }


//    for (int i=0; i<n_part; i++)
//    {
//        std::cout << "particle:" << i << "  " << "X:" << vec_part[i].getPos()[0] << ","
//               << vec_part[i].getPos()[1] << "," << vec_part[i].getPos()[2] <<"\n" ;
//    }

////Bond initialization

    double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;
    int IDS[2];

    posi = vec_part[0].getPos();
    radi = vec_part[0].getDiameter()/2.;
    radi_cr = radi*lambda;

    posj = vec_part[1].getPos();
    radj = vec_part[1].getDiameter()/2.;
    radj_cr = radj*lambda;

    dist = distance(posi, posj);

    if (dist < (radi_cr+radj_cr))
    {
        double radius;
        //double const kn = E/(radi+radj);
        //double const ks = kn/2.5;
        double const kn = sqrt(2.)*E*radi/(2.*(1.-2.*nu));
        double const ks = kn/2.5;
//        double const kb = (1-2.*nu)*radi*radi/(24.*(1.-nu))*kn;
//        double const kt = (1-2.*nu)*radi*radi/(24.*(1.-nu))*ks;

        radius = lambda*std::min(radi,radj);

        IDS[0] = 0;
        IDS[1] = 1;

        vec_bond.setRad(radius);
        vec_bond.setStiffNorm(kn);
        vec_bond.setStiffShear(ks);
//        vec_bond.setStiffBend(kb);
//        vec_bond.setStiffTor(kt);
        vec_bond.setIds(IDS);
    }


    /*Starts the dynamics*/

    for (int i=0; i<n; i++)
    {
        /*first timestep, initialization*/

        if(i==0)
        {
            double *x;

            x = vec_part[0].getPos();

            x[0] += 1.e-6;

            vec_part[0].setPos(x);
        }

        else
        {
            /*New positions with integration of EOM due to bond forces*/

            for (int j=0; j<n_part; j++)
            {
                double *x, *v, *f, *q, *om, *t;
                double *v_old, *f_old, *om_old, *t_old, *q_old;
                double x_new[3], v_new[3], q_new[4], om_new[3];
                double v_act[3], q_act[3], om_act[3];
                double m, d, inertia;

                x = vec_part[j].getPos();
                v = vec_part[j].getVel();
                f = vec_part[j].getForce();

                q = vec_part[j].getQuat();
                om = vec_part[j].getRot();
                t = vec_part[j].getTorque();

                v_old = vec_part[j].getVelOld();
                f_old = vec_part[j].getForceOld();
                om_old = vec_part[j].getRotOld();
                q_old = vec_part[j].getQuatOld();
                t_old = vec_part[j].getTorqueOld(),

                m = vec_part[j].getMass();
                d = vec_part[j].getDiameter();
                inertia = (2./5.)*m*pow((d/2.),2.);

                if(i==1)
                {
                    q_new[0] = q[0] - 0.5*dt*(om[0]*q[1] + om[1]*q[2] + om[2]*q[3]);
                    q_new[1] = q[1] + 0.5*dt*(om[0]*q[0] - om[1]*q[3] + om[2]*q[2]);
                    q_new[2] = q[2] + 0.5*dt*(om[0]*q[3] + om[1]*q[0] - om[2]*q[1]);
                    q_new[3] = q[3] + 0.5*dt*(-om[0]*q[2] + om[1]*q[1] + om[2]*q[0]);

                    q_act[0] = q[0];

                    for (int k=0; k<3; k++)
                    {
                        q_act[k+1] = q[k+1];
                        om_act[k] = om[k];
                        v_act[k] = v[k];

                        x_new[k] = x[k] + v[k]*dt;
                        v_new[k] = v[k] + (dt/m)*f[k];
                        om_new[k] = om[k] + (dt/inertia)*t[k];
                    }
                }

                else
                {
                    q_new[0] = q[0] - 0.25*dt*(3.*(om[0]*q[1] + om[1]*q[2] + om[2]*q[3]) - (om_old[0]*q_old[1] + om_old[1]*q_old[2] + om_old[2]*q_old[3]));
                    q_new[1] = q[1] + 0.25*dt*(3.*(om[0]*q[0] - om[1]*q[3] + om[2]*q[2]) - (om_old[0]*q_old[0] - om_old[1]*q_old[3] + om_old[2]*q_old[2]));
                    q_new[2] = q[2] + 0.25*dt*(3.*(om[0]*q[3] + om[1]*q[0] - om[2]*q[1]) - (om_old[0]*q_old[3] + om_old[1]*q_old[0] - om_old[2]*q_old[1]));
                    q_new[3] = q[3] + 0.25*dt*(3.*(-om[0]*q[2] + om[1]*q[1] + om[2]*q[0]) - (-om_old[0]*q_old[2] + om_old[1]*q_old[1] + om_old[2]*q_old[0]));

                    q_act[0] = q[0];

                    for (int k=0; k<3; k++)
                    {
                        q_act[k+1] = q[k+1];
                        om_act[k] = om[k];
                        v_act[k] = v[k];

                        x_new[k] = x[k] + 0.5*dt*(3.*v[k] - v_old[k]);
                        v_new[k] = v[k] + 0.5*(dt/m)*(3.*f[k] - f_old[k]);
                        om_new[k] = om[k] + 0.5*(dt/inertia)*(3.*t[k] - t_old[k]);
                    }
                }

                vec_part[j].setPos(x_new);
                vec_part[j].setVel(v_new);
                vec_part[j].setQuat(q_new);
                vec_part[j].setRot(om_new);

                vec_part[j].setQuatOld(q_act);
                vec_part[j].setRotOld(om_act);
                vec_part[j].setTorqueOld(t);
                vec_part[j].setVelOld(v_act);
                vec_part[j].setForceOld(f);
            }
        }

        //Forces on bonds

        Bond *bondino;

        bondino = &vec_bond;
        bond_fdm(bondino, vec_part[0], vec_part[1]);
        //f_bonds(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt, gamma);

        if (i==1)
            std::cout << "normal force: " << vec_bond.getForce_n()[0] << "," << vec_bond.getForce_n()[1] << "," << vec_bond.getForce_n()[2] << "\n"
                      << "shear-trans force: " << vec_bond.getForce_st()[0] << "," << vec_bond.getForce_st()[1] << "," << vec_bond.getForce_st()[2] << "\n"
//                      << "shear-rot force: " << vec_bond.getForce_sr()[0] << "," << vec_bond.getForce_sr()[1] << "," << vec_bond.getForce_sr()[2] << "\n"
                      << "shear-trans torque: " << vec_bond.getTorque_st()[0] << "," << vec_bond.getTorque_st()[1] << "," << vec_bond.getTorque_st()[2] << "\n";
//                      << "shear-rot torque: " << vec_bond.getTorque_sr()[0] << "," << vec_bond.getTorque_sr()[1] << "," << vec_bond.getTorque_sr()[2] << "\n"
//                      << "bend-rot torque: " << vec_bond.getTorque_br()[0] << "," << vec_bond.getTorque_br()[1] << "," << vec_bond.getTorque_br()[2] << "\n"
//                      << "tors-rot torque: " << vec_bond.getTorque_tr()[0] << "," << vec_bond.getTorque_tr()[1] << "," << vec_bond.getTorque_tr()[2] << "\n";

        //Forces on particles

        for (int j=0; j<n_part; j++)
        {
            double F_bond_rel[4], T_bond[3], F_bond[3], *F_bond_quat;

            //Bonds contribution

            int *IDS;

            IDS = vec_bond.getIds();

            double  *q, *q_inv;

            q = vec_part[IDS[1]].getQuat();
            q_inv = inverse(q);

            if (j == IDS[0])
            {
                F_bond_rel[0] = 0.;
                F_bond_rel[1] = -(vec_bond.getForce_n()[0] + vec_bond.getForce_st()[0]);
                F_bond_rel[2] = -(vec_bond.getForce_n()[1] + vec_bond.getForce_st()[1]);
                F_bond_rel[3] = -(vec_bond.getForce_n()[2] + vec_bond.getForce_st()[2]);

                F_bond_quat = rotation(q, F_bond_rel, q_inv, 3);

                F_bond[0] = F_bond_quat[1];
                F_bond[1] = F_bond_quat[2];
                F_bond[2] = F_bond_quat[3];

                T_bond[0] = vec_bond.getTorque_st()[0];
                T_bond[1] = vec_bond.getTorque_st()[1];
                T_bond[2] = vec_bond.getTorque_st()[2];
            }

            else if (j == IDS[1])
            {
                F_bond_rel[0] = 0.;
                F_bond_rel[1] = (vec_bond.getForce_n()[0] + vec_bond.getForce_st()[0]);
                F_bond_rel[2] = (vec_bond.getForce_n()[1] + vec_bond.getForce_st()[1]);
                F_bond_rel[3] = (vec_bond.getForce_n()[2] + vec_bond.getForce_st()[2]);

                F_bond_quat = rotation(q, F_bond_rel, q_inv, 3);

                F_bond[0] = F_bond_quat[1];
                F_bond[1] = F_bond_quat[2];
                F_bond[2] = F_bond_quat[3];

                T_bond[0] = vec_bond.getTorque_st()[0];
                T_bond[1] = vec_bond.getTorque_st()[1];
                T_bond[2] = vec_bond.getTorque_st()[2];
            }


            vec_part[j].setForce(F_bond);
            vec_part[j].setTorque(T_bond);
        }


        for(int j=0; j<n_part; j++)
        {
            if (i == 1 || i == 2)
            {
                std::cout << j << ","
                          << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2]
                          << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                          << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << "\n";
            }
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == count*n_out)
            {
                out_file << i*dt << "," << j << ","
                         << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2] << ","
                         << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2] << ","
                         << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << ","
                         << vec_part[j].getRot()[0] << "," << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2] << ","
                         << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
                flag = true;
            }
        }

        if (flag)
        {
            count++;
            flag = false;
        }
    }

    for (int i=0; i<n_part; i++)
    {
        std::cout << i << "," << vec_part[i].getPos()[0] << ","
                              << vec_part[i].getPos()[1] << ","
                              << vec_part[i].getPos()[2] << "\n";
    }

    return 0;
}
