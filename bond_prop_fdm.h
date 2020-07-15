class Bond
{
public:
  Bond() :
    sigmac_(0.),
    tauc_(0.),
    rad_(0.),
    kn_(0.),
    ks_(0.),
    kb_(0.),
    kt_(0.),
    partid_{0, 0},
    force_n_{0., 0., 0.},
    force_st_{0., 0., 0.},
    force_sr_{0., 0., 0.},
    torque_n_{0., 0., 0.},
    torque_st_{0., 0., 0.},
    torque_sr_{0., 0., 0.},
    torque_br_{0., 0., 0.},
    torque_tr_{0., 0., 0.}
  {}

  Bond(double sigma, double tau, double radius, double kn, double ks, double area, double kb, double kt) :
    sigmac_(sigma),
    tauc_(tau),
    rad_(radius),
    kn_(kn),
    ks_(ks),
    kb_(kb),
    kt_(kt)
  {}

  Bond(int i, int j) :
    partid_{i, j}
  {}

  Bond(double fnx, double fny, double fnz, double fstx, double fsty, double fstz, double fsrx, double fsry, double fsrz) :
    force_n_{fnx, fny, fnz},
    force_st_{fstx, fsty, fstz},
    force_sr_{fsrx, fsry, fsrz}
  {}

  Bond(double tnx, double tny, double  tnz, double tstx, double tsty, double tstz, double tsrx, double tsry, double tsrz, double tbrx, double tbry, double tbrz, double ttrx, double ttry, double ttrz) :
    torque_n_{tnx, tny, tnz},
    torque_st_{tstx, tsty, tstz},
    torque_sr_{tsrx, tsry, tsrz},
    torque_br_{tbrx, tbry, tbrz},
    torque_tr_{ttrx, ttry, ttrz}
  {}

  ~Bond() {}

  double getSigma()
  { return sigmac_; }

  double getTau()
  { return tauc_; }

  double getRad()
  { return rad_; }

  double getStiffNorm()
  { return kn_; }

  double getStiffShear()
  { return ks_; }

  double getStiffBend()
  { return kb_; }

  double getStiffTor()
  { return kt_; }

  int* getIds()
  { return partid_; }

  double* getForce_n()
  { return force_n_; }

  double* getForce_st()
  { return force_st_; }

  double* getForce_sr()
  { return force_sr_; }

  double* getTorque_st()
  { return torque_st_; }

  double* getTorque_n()
  { return torque_n_; }

  double* getTorque_sr()
  { return torque_sr_; }

  double* getTorque_br()
  { return torque_br_; }

  double* getTorque_tr()
  { return torque_tr_; }


  void setSigma(double sigma)
  { sigmac_ = sigma; }

  void setTau(double tau)
  { tauc_ = tau; }

  void setRad(double radius)
  { rad_ = radius; }

  void setStiffNorm(double kn)
  { kn_ = kn; }

  void setStiffShear(double ks)
  { ks_ = ks; }

  void setStiffBend(double kb)
  { kb_ = kb; }

  void setStiffTor(double kt)
  { kt_ = kt; }

  void setIds(int* ids)
  {
    partid_[0] = ids[0];
    partid_[1] = ids[1];
  }

  void setForce_n(double* fn)
  {
    force_n_[0] = fn[0];
    force_n_[1] = fn[1];
    force_n_[2] = fn[2];
  }

  void setForce_st(double* fst)
  {
    force_st_[0] = fst[0];
    force_st_[1] = fst[1];
    force_st_[2] = fst[2];
  }

  void setForce_sr(double* fsr)
  {
    force_sr_[0] = fsr[0];
    force_sr_[1] = fsr[1];
    force_sr_[2] = fsr[2];
  }

  void setTorque_st(double* tst)
  {
    torque_st_[0] = tst[0];
    torque_st_[1] = tst[1];
    torque_st_[2] = tst[2];
  }

  void setTorque_n(double* tn)
  {
    torque_n_[0] = tn[0];
    torque_n_[1] = tn[1];
    torque_n_[2] = tn[2];
  }

  void setTorque_sr(double* tsr)
  {
    torque_sr_[0] = tsr[0];
    torque_sr_[1] = tsr[1];
    torque_sr_[2] = tsr[2];
  }

  void setTorque_br(double* tbr)
  {
    torque_br_[0] = tbr[0];
    torque_br_[1] = tbr[1];
    torque_br_[2] = tbr[2];
  }

  void setTorque_tr(double* ttr)
  {
    torque_tr_[0] = ttr[0];
    torque_tr_[1] = ttr[1];
    torque_tr_[2] = ttr[2];
  }


private:
  double sigmac_, tauc_, rad_, kn_, ks_, kb_, kt_;
  int partid_[2];
  double force_n_[3], force_st_[3], force_sr_[3], torque_n_[3], torque_st_[3], torque_sr_[3], torque_br_[3], torque_tr_[3];

};




