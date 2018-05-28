cccccccccccccccc include cparam.h cccccccccccccccccccccccccccccccc
c
c       parameter file to over-ride the parameter settings for
c       Dynamic-Thermodynamic coupling if RG model
c       
c       Developed on:           Feb-2017
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        relax_on_off     = 0.0    ! 0.0 is for temp, salt relaxation off
                                  ! utilized in tracer.polar.F
                                  ! Default is 1.0

        alpha_rgm = 0.0  *      mod(loop,2)  
        alpha_rgt = 1.0  * (1 - mod(loop,2)) 
        beta_mldm = 0.0  *      mod(loop,2)  
        beta_mldt = 1.0  * (1 - mod(loop,2)) 

				  ! coupling parameters
				  ! alpha -> RG gets MYM
				  ! beta  -> MYM transport by RG
				  ! m=momentum, t=tracer
				  ! utilized in couple.F
                                  ! Default is all 1.0 with 
                                  ! mod(loop,2) and (1-mod(loop,2))
				  ! as alternate for m and t.
			!CAUTION  ! If beta_mldm, beta_mldt =0.0 then tracer
                                  ! in RG will not advect within MYM domain


        gama = 1.0     		  ! 1.0 is for entrainment in MYM from
				  ! divergence of RG model on, 0.0 is off.
				  ! utilized in interp_extrap.F
                                  ! Default is 1.0

	diffuse_MY = 1.0e3        ! Diffusion for MYM transport

	diffuse_tr = 1.0e3	  ! Diffusion for tracer transport in RG

	rdrag = 1.0/(100*day2sec)  ! Drag between two layers in days

	relax_South = 1.0	  ! If 1.0 the sothern boundary temp and salt
				  ! relaxes linearly to observations
	widS = 3.0                ! to a width of widS grid thickness

	relax_East = 1.0	  ! If 1.0 the eastern boundary temp and salt
				  ! relaxes linearly to observations
	widE = 3.0                ! to a width of widS grid thickness
