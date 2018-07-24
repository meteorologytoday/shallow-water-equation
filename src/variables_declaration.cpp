float dx, dy, Lx, Ly;

float *Q, *vort, *divg, *geop, *bg_vort, 
      *u2, *u, *u_divg, *u_vort, 
      *v2, *v, *v_divg, *v_vort, 
      *K, *E, 
      *absvort, *absvort_u, *absvort_v, 
      *geop_u, *geop_v, 
       *dvortdx, *dvortdy,
	  *dvortdt, *ddivgdt, *dgeopdt, *workspace;

fftwf_complex *Q_c,
              *vort_c0, *vort_c, *lvort_c, 
              *divg_c0, *divg_c, *ldivg_c, 
              *geop_c0,    *geop_c,    *lgeop_c,    
              *absvort_u_c, *absvort_v_c,
              *geop_u_c, *geop_v_c,
              *E_c,
              *chi_c, *psi_c, 
              *tmp_c, 
              *copy_for_c2r,
              *dvortdt_c, *ddivgdt_c, *dgeopdt_c;
              
fftwf_complex *rk4_vort_c[4], *rk4_divg_c[4], *rk4_geop_c[4];

fftwf_plan p_fwd_vort,    p_bwd_vort,
           p_fwd_divg,    p_bwd_divg,
           p_fwd_geop,    p_bwd_geop,
		   p_bwd_dvortdx, p_bwd_dvortdy,
		   p_bwd_u,       p_bwd_v,
		   p_bwd_u_vort,  p_bwd_v_vort,
		   p_bwd_u_divg,  p_bwd_v_divg,
		   p_fwd_dvortdt, 
		   p_fwd_ddivgdt, 
		   p_fwd_dgeopdt, 
           p_fwd_absvort_u,
           p_fwd_absvort_v,
           p_fwd_geop_u, p_fwd_geop_v,
           p_fwd_E, p_fwd_Q;



fftwf_operation<XPTS,YPTS> fop(LX, LY);

VectorOperation<GRIDS> vop;

char filename[256];


