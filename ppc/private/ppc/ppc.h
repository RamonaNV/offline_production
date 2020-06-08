namespace xppc{
  float xrnd();
  float grnd();

  void start();
  void stop();
  void choose(int);
  void ini();
  void fin();
  void eout();

  void flone(unsigned long long);

  size_t getMaxBunchSize();
  size_t getWorkgroupSize();

  struct mcid:std::pair<int,unsigned long long>{
    int frame;
  };

  struct DOM{
    float r[3];
  };

  struct ikey{
    int str, dom;

    bool isinice() const{
      return str>0 && dom>=1 && dom<=60;
    }

    bool operator< (const ikey & rhs) const {
      return str == rhs.str ? dom < rhs.dom : str < rhs.str;
    }

    bool operator!= (const ikey & rhs) const {
      return str != rhs.str || dom != rhs.dom;
    }
  };

  struct OM:DOM,ikey{};
  extern std::vector<OM> i3oms;
  extern std::map<ikey, float> hvs;
  extern std::map<ikey, std::pair<float, int> > rdes;

  void initialize(float);
  const DOM& flset(int, int);
  void flshift(float [], float [], float * = NULL);

  struct ihit{
    ikey omkey;
    mcid track;
    float time;
    float dir;

    bool operator< (const ihit & rhs) const {
      return
	track!=rhs.track ? track<rhs.track :
	omkey!=rhs.omkey ? omkey<rhs.omkey :
	dir!=rhs.dir ? dir<rhs.dir :
	time<rhs.time;
    }
  };

  void set_res(float);
  void set_res(float, float);

  struct pout{
    float r[4], n[4];
  };

  typedef std::map<ihit, std::vector<pout> > outz;
  extern outz hitz;

  void efin();

  void sett(float, float, float, std::pair<int,unsigned long long>, int);
  template <class T> void addp(float, float, float, float, float, T, float = 1);

  void addp_clst(float, float, float, float, unsigned long long, float, float);
  void addp_mopo(float, float, float, float, float, float, float, float, float);
}
