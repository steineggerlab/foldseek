
typedef struct {
  int bins[3];
  float data[8][3];
} nco_struct;

extern nco_struct nco_stat[];
extern nco_struct nco_stat_pro[];
extern float rot_stat_coords[][3];
extern int rot_stat_idx[][6];
