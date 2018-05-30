class CrystalInfo{
 public:
  // methods
  CrystalInfo();
  ~CrystalInfo();

  // variables
  //DetId id;
  uint32_t rawId;
  int ieta;
  int iphi;
  int ix;
  int iy;
  double energy;
  double time;
  double timeErr;
  int recoFlag;
};

CrystalInfo::CrystalInfo()
{}
CrystalInfo::~CrystalInfo()
{}
