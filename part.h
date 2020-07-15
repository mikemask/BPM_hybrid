class Part
{
public:
    Part() :
        r_{0., 0., 0.},
        r0_{0., 0., 0.},
        q_{1., 0., 0., 0.},
        q_old_{1., 0., 0., 0.},
        U_{0., 0., 0.},
        U_old_{0., 0., 0.},
        om_{0., 0., 0.},
        om_old_{0., 0., 0.},
        f_{0., 0., 0.},
        f_old_{0., 0., 0.},
        t_{0., 0., 0.},
        t_old_{0., 0., 0.},
        m_(0.),
        d_(0.),
        id_(-1)
    {}

    Part(double x, double y, double z, double x0, double y0, double z0, double q0, double q1, double q2, double q3, double u, double v, double w, double omx, double omy, double omz, double fx, double fy, double fz, double tx, double ty, double tz) :
        r_{x, y, z},
        r0_{x0, y0, z0},
        q_{q0, q1, q2, q3},
        U_{u, v, w},
        om_{omx, omy, omz},
        f_{fx, fy, fz},
        t_{tx, ty, tz}
    {}

    Part(double u_old, double v_old, double w_old, double q0_old, double q1_old, double q2_old, double q3_old, double fx_old, double fy_old, double fz_old, double omx_old, double omy_old, double omz_old, double tx_old, double ty_old, double tz_old) :
        U_old_{u_old, v_old, w_old},
        q_old_{q0_old, q1_old, q2_old, q3_old},
        f_old_{fx_old, fy_old, fz_old},
        om_old_{omx_old, omy_old, omz_old},
        t_old_{tx_old, ty_old, tz_old}
    {}

    Part(double m, double d, int id) :
        m_(m),
        d_(d),
        id_(id)
    {}

    ~Part() {}

    double* getPos()
    { return r_; }

    double* getPosIn()
    { return r0_; }

    double* getQuat()
    { return q_; }

    double* getQuatOld()
    { return q_old_; }

    double* getVel()
    { return U_; }

    double* getVelOld()
    { return U_old_; }

    double* getRot()
    { return om_; }

    double* getRotOld()
    { return om_old_; }

    double* getForce()
    { return f_; }

    double* getForceOld()
    { return f_old_; }

    double* getTorque()
    { return t_; }

    double* getTorqueOld()
    { return t_old_; }

    double getMass()
    { return m_; }

    double getDiameter()
    { return d_; }

    double getId()
    { return id_; }


    void setPos(double* r)
    {
        r_[0] = r[0];
        r_[1] = r[1];
        r_[2] = r[2];
    }

    void setPosIn(double* r0)
    {
        r0_[0] = r0[0];
        r0_[1] = r0[1];
        r0_[2] = r0[2];
    }

    void setQuat(double* q)
    {
        q_[0] = q[0];
        q_[1] = q[1];
        q_[2] = q[2];
        q_[3] = q[3];
    }

    void setQuatOld(double* q_old)
    {
        q_old_[0] = q_old[0];
        q_old_[1] = q_old[1];
        q_old_[2] = q_old[2];
        q_old_[3] = q_old[3];
    }

    void setVel(double* U)
    {
        U_[0] = U[0];
        U_[1] = U[1];
        U_[2] = U[2];
    }

    void setVelOld(double* U_old)
    {
        U_old_[0] = U_old[0];
        U_old_[1] = U_old[1];
        U_old_[2] = U_old[2];
    }

    void setRot(double* om)
    {
        om_[0] = om[0];
        om_[1] = om[1];
        om_[2] = om[2];
    }

    void setRotOld(double* om_old)
    {
        om_old_[0] = om_old[0];
        om_old_[1] = om_old[1];
        om_old_[2] = om_old[2];
    }

    void setForce(double* f)
    {
        f_[0] = f[0];
        f_[1] = f[1];
        f_[2] = f[2];
    }

    void setForceOld(double* f_old)
    {
        f_old_[0] = f_old[0];
        f_old_[1] = f_old[1];
        f_old_[2] = f_old[2];
    }

    void setTorque(double* t)
    {
        t_[0] = t[0];
        t_[1] = t[1];
        t_[2] = t[2];
    }

    void setTorqueOld(double* t_old)
    {
        t_old_[0] = t_old[0];
        t_old_[1] = t_old[1];
        t_old_[2] = t_old[2];
    }

    void setMass(double m)
    { m_ = m; }

    void setDiameter(double d)
    { d_ = d; }

    void setId(int id)
    { id_ = id; }


private:
    double r_[3], r0_[3], q_[4], q_old_[4], U_[3], U_old_[3], om_[3], om_old_[3], f_[3], f_old_[3], t_[3], t_old_[3], m_, d_;
    int id_;
};
