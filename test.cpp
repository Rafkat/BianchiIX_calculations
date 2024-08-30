//struct rhs_van {
//    Doub eps;
//    rhs_van(Doub epss) : eps(epss) {}
//    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
//        dydx[0] = y[1];
//        dydx[1] = ((1.0 - y[0]*y[0])*y[1]-y[0])/eps;
//    }
//};


int test() {
    const Int nvar {2};
    const Doub atol {1.0e-3};
    const Doub rtol {atol};
    const Doub h1{0.01};
    const Doub hmin {0.0};
    const Doub x1{0.0};
    const Doub x2{-2.0};
    VecDoub ystart(nvar);
    ystart[0]=2.0;
    ystart[1]=0.0;
    Output out(100);
    rhs_van d(1.0e-3);
    Odeint<StepperDopr5<rhs_van> > ode(ystart,x1,x2,atol,rtol,h1,hmin,out,d);
    ode.integrate();

    for (Int i=0;i<out.count;i++)
        cout << out.xsave[i] << " " << out.ysave[0][i]
                << " " << out.ysave[1][i] << endl;

    return 0;
}