    printf("Allocating spaces...");

	Q          = malloc_field_re();
	vort       = malloc_field_re();
	divg       = malloc_field_re();
	geop       = malloc_field_re();
	bg_vort    = malloc_field_re();

	u         = malloc_field_re();
	u2        = malloc_field_re();
	u_divg    = malloc_field_re();
	u_vort    = malloc_field_re();

    v         = malloc_field_re();
    v2        = malloc_field_re();
    v_divg    = malloc_field_re();
    v_vort    = malloc_field_re();


    K         = malloc_field_re();
    E         = malloc_field_re();

    absvort   = malloc_field_re();
    absvort_u = malloc_field_re();
    absvort_v = malloc_field_re();
    
    geop_u    = malloc_field_re();
    geop_v    = malloc_field_re();
    
    dvortdx   = malloc_field_re();
    dvortdy   = malloc_field_re();
    
    dvortdt   = malloc_field_re();
    ddivgdt   = malloc_field_re();
    dgeopdt   = malloc_field_re();
    
    workspace = malloc_field_re();


    // complex numbers
    
    vort_c0 = malloc_field_im();
    vort_c  = malloc_field_im(); 
    lvort_c = malloc_field_im();

    divg_c0 = malloc_field_im();
    divg_c  = malloc_field_im(); 
    ldivg_c = malloc_field_im();

    geop_c0 = malloc_field_im();
    geop_c  = malloc_field_im(); 
    lgeop_c = malloc_field_im();
    
    absvort_u_c  = malloc_field_im();
    absvort_v_c  = malloc_field_im();
    geop_u_c     = malloc_field_im();
    geop_v_c     = malloc_field_im();
    E_c          = malloc_field_im();
    Q_c          = malloc_field_im();
    chi_c        = malloc_field_im();
    psi_c        = malloc_field_im();
    tmp_c        = malloc_field_im();
    copy_for_c2r = malloc_field_im();

    dvortdt_c    = malloc_field_im();
    ddivgdt_c    = malloc_field_im();
    dgeopdt_c    = malloc_field_im();
    
    for(int k=0; k < 4; ++k) {
        rk4_vort_c[k] = malloc_field_im();
        rk4_divg_c[k] = malloc_field_im();
        rk4_geop_c[k] = malloc_field_im();
    }

    printf("done.\n");
    printf("Making plans...");
    // plans
    p_fwd_Q          = crt_fwd_plan(Q, Q_c);
    p_fwd_E          = crt_fwd_plan(E, E_c);
    
    p_fwd_vort       = crt_fwd_plan(vort, vort_c);
    p_bwd_vort       = crt_bwd_plan(vort, vort_c);

    p_fwd_divg       = crt_fwd_plan(divg, divg_c);
    p_bwd_divg       = crt_bwd_plan(divg, divg_c);

    p_fwd_geop       = crt_fwd_plan(geop, geop_c);
    p_bwd_geop       = crt_bwd_plan(geop, geop_c);

    p_bwd_dvortdx    = crt_bwd_plan(dvortdx, tmp_c);
    p_bwd_dvortdy    = crt_bwd_plan(dvortdy, tmp_c);
    
    p_bwd_u          = crt_bwd_plan(u     , tmp_c);
    p_bwd_u_vort     = crt_bwd_plan(u_vort, tmp_c);
    p_bwd_u_divg     = crt_bwd_plan(u_divg, tmp_c);

    p_bwd_v          = crt_bwd_plan(v     , tmp_c);
    p_bwd_v_vort     = crt_bwd_plan(v_vort, tmp_c);
    p_bwd_v_divg     = crt_bwd_plan(v_divg, tmp_c);


    p_fwd_dvortdt    = crt_fwd_plan(dvortdt, dvortdt_c);
    p_fwd_ddivgdt    = crt_fwd_plan(ddivgdt, ddivgdt_c);
    p_fwd_dgeopdt    = crt_fwd_plan(dgeopdt, dgeopdt_c);
    
    p_fwd_absvort_u  = crt_fwd_plan(absvort_u, absvort_u_c);
    p_fwd_absvort_v  = crt_fwd_plan(absvort_v, absvort_v_c);

    p_fwd_geop_u     = crt_fwd_plan(geop_u, geop_u_c);
    p_fwd_geop_v     = crt_fwd_plan(geop_v, geop_v_c);


    printf("done.\n");

	// read input
	Lx = LX;
	Ly = LY;
	dx = Lx / XPTS;
	dy = Ly / YPTS;

