from ngsolve import *
from ngsolve.solvers import *

def SolMaxStefan(mesh,order,uin,DD,eps,steps,bndc=[],bnd="",basemeshsize=1,nitsche=1):
    N=uin.dim
    D=mesh.dim
    normal = specialcf.normal(D)

    V = H1(mesh,order=order)
    X = FESpace(([V if i<N else V for i in range(N+N*(D-1))]))

    trial = X.TrialFunction()
    test = X.TestFunction()
    w = CoefficientFunction( tuple([trial[i] for i in range(N)]) )
    v = CoefficientFunction( tuple([ test[i] for i in range(N)]) )
    J = CoefficientFunction( tuple([trial[i+N] for i in range(N*(D-1))]), dims=(N,(D-1)) )
    K = CoefficientFunction( tuple([ test[i+N] for i in range(N*(D-1))]), dims=(N,(D-1)) )

    dtw = CoefficientFunction(tuple([grad(trial[i])[D-1] for i in range(N)]))
    dtv = CoefficientFunction(tuple([grad( test[i])[D-1] for i in range(N)]))
    dxw = [0]*(D-1)
    dnw = [0]*(D-1)
    dnv = [0]*(D-1)
    dxv = [0]*(D-1)
    JN = [0]*(D-1)
    KN = [0]*(D-1)
    for i in range((D-1)):
        dxw[i] = CoefficientFunction(tuple([grad(trial[j])[i] for j in range(N)]))
        dnw[i] = CoefficientFunction(tuple([grad(trial[j]).Trace()[i]*normal[i] for j in range(N)]))
        dnv[i] = CoefficientFunction(tuple([grad(test[j]).Trace()[i]*normal[i] for j in range(N)]))
        dxv[i] = CoefficientFunction(tuple([grad( test[j])[i] for j in range(N)]))
        JN[i] = CoefficientFunction(tuple([J[d,i] for d in range(N)]))
        KN[i] = CoefficientFunction(tuple([K[d,i] for d in range(N)]))

    u = (1/(1+sum(exp(w))))*exp(w)

    mat_M = CoefficientFunction(tuple([u[i]/(DD[i][j])-u[i]/(DD[i][N])-(i==j)*(sum([u[k]/(DD[i][k]) for k in range(N)]) + (1-sum(u))/(DD[i][N])) for i in range(N) for j in range(N)]), dims=(N,N))

    spp = CoefficientFunction(tuple([1/u[i]+1/(1-sum(u)) if i==j else 1/(1-sum(u)) for i in range(N) for j in range(N)]),dims=(N,N))

    uinit = GridFunction(X)
    for i in range(N):
        uinit.components[i].Set(uin[i], BND, definedon=mesh.Boundaries('bottom'))

    for n in range(2*steps):
        print("Step",n+1,"/",2*steps)
        a = BilinearForm(X)
        a += SymbolicBFI(- (-1)**(n) * InnerProduct(u,dtv)
                         + eps*InnerProduct(w,v)
                         )
        for j in range(D-1):
            a += SymbolicBFI(- InnerProduct(JN[j],dxv[j])
                             - InnerProduct(dxw[j],KN[j])
                             + InnerProduct(spp*mat_M*(JN[j]),KN[j])
                            )
            a += SymbolicBFI(eps*InnerProduct(dxw[j],dxv[j]))
        a += SymbolicBFI(eps*InnerProduct(dtw,dtv))

        bottom = 'top' if n%2 else 'bottom';
        top = 'bottom' if n%2 else 'top';
        a += SymbolicBFI( InnerProduct(u,v), BND, definedon=mesh.Boundaries(top))
        a += SymbolicBFI(- InnerProduct(CoefficientFunction(tuple(uinit.components[i] for i in range(N))),v), BND, definedon=mesh.Boundaries(bottom))

        a += SymbolicBFI(nitsche*pow(basemeshsize,-1) * InnerProduct(u-bndc,v), BND, definedon=mesh.Boundaries(bnd))
        for j in range(D-1):
            a += SymbolicBFI(InnerProduct(JN[j]*normal[j],v), BND, definedon=mesh.Boundaries(bnd))
            a += SymbolicBFI(InnerProduct(u-bndc,KN[j]*normal[j]), BND, definedon=mesh.Boundaries(bnd))

        gfu = GridFunction(X)

        with TaskManager():
            Newton(a,gfu,freedofs=X.FreeDofs(),maxit=100,maxerr=1e-11,inverse="umfpack",dampfactor=1,printing=True)

        for i in range(N):
            uinit.components[i].Set((1/(1+sum(exp(gfu.components[j]) for j in range(N))))*exp(gfu.components[i]))

    return uinit
