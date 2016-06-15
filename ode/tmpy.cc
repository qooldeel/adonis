#ifdef USE_2_STAGE_2ND_ORDER_ROSENBROCK_METHOD
	//1st linear system
	//possible clipping performed
	Repairer.repair(y_n,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	Repairer.info(i,-1, "1st Rosenbrock sys: ");
	
	g = fun_(y_n);
	Gprime.evaluate_jacobian(y_n);
	Gprime.compute_G_prime(h,1.,1+1./sqrt(2.));
	//ScaleSys.scale(Gprime.value_pointer(),g.begin(),y_n,rtol,ntol);
	//Gprime2.get_matrix() = Gprime.get_matrix(); //assign iteration matrix, before possibly being ordered (needed for scaling)!
	Gprime.reorder_values();
	Gp2 = Gprime.get_matrix();   //use G' for second lin. system 
	LSS.solve(); //g (a.k.a. k1) is overwritten now with solution of 1st sys, so might be G'


	//2nd linear system
	ev2 = y_n + h*g;       //1st order consistend at t = t_n+1 ==> cheap local error estimation ofr step size control [4, p. 1461]
	//possible clipping performed
	Repairer.repair(ev2,excessSpecIndex); //repair me or not (obey bounds on species etc.)
	Repairer.info(i,-1,"2nd Rosenbrock sys: ");
	k2 = fun_(ev2) - 2.*g;
	//ScaleSys.scale(Gprime2.value_pointer(),k2.begin(),y_n,rtol,ntol);
	//Gprime2.reorder_values();
	LSS2.solve();  // solves Gp2 · x = k2 and overwrites k2 with sol x

	//! update solution
	y_n += 1.5*h*g + 0.5*h*k2;

#else //NEWTON METHOD for 1-step methods used
